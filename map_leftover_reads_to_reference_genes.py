#!/usr/bin/env python3

import argparse
import gzip
import glob
import shutil
import subprocess
from pathlib import Path

import yaml
import pandas as pd


FASTQ_SUFFIXES = [
    ".decoy_filtered.PAO1_unmapped.fastq.gz",
    ".decoy_filtered.PAO1_unmapped.fastq",
    ".fastq.gz",
    ".fastq",
    ".fq.gz",
    ".fq",
]


def run(cmd, stdout_path=None, stderr_path=None):
    print("Running:", " ".join(map(str, cmd)))

    stdout = open(stdout_path, "w") if stdout_path else None
    stderr = open(stderr_path, "w") if stderr_path else None

    try:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
    finally:
        if stdout:
            stdout.close()
        if stderr:
            stderr.close()


def check_programs():
    missing = []

    for program in ["bowtie2", "bowtie2-build", "featureCounts"]:
        if shutil.which(program) is None:
            missing.append(program)

    if missing:
        raise RuntimeError(
            "Missing required programs in PATH: "
            + ", ".join(missing)
            + "\nActivate/install the right conda env first."
        )


def sample_name_from_fastq(fastq_path):
    name = Path(fastq_path).name

    for suffix in FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name.replace(suffix, "")

    return Path(fastq_path).stem


def find_sample_ref_match_csv(config):
    configured = (
        config.get("sample_ref_match_csv")
        or config.get("sample_reference_match_csv")
        or config.get("input", {}).get("sample_ref_match_csv")
        or config.get("input", {}).get("sample_reference_match_csv")
    )

    if configured:
        return Path(configured).expanduser()

    default_path = Path("sample_ref_match.csv")
    return default_path if default_path.exists() else None


def load_sample_reference_matches(match_csv):
    """Load sample-to-reference assignments from a two-column CSV file."""
    matches = {}
    df = pd.read_csv(match_csv, header=None, comment="#", names=["sample", "reference"])

    for _, row in df.dropna(how="all").iterrows():
        sample_value = str(row["sample"]).strip()
        reference = str(row["reference"]).strip()

        if not sample_value or not reference:
            continue

        sample_path = Path(sample_value)
        sample_keys = {
            sample_value,
            str(sample_path.expanduser()),
            sample_path.name,
            sample_name_from_fastq(sample_path),
        }

        for key in sample_keys:
            matches[key] = reference

    return matches


def reference_for_sample(fastq, sample, matches):
    if not matches:
        return None

    fastq = Path(fastq)
    for key in [str(fastq), str(fastq.expanduser()), fastq.name, sample]:
        if key in matches:
            return matches[key]

    return None


def open_maybe_gz(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def count_fastq_reads(fastq):
    lines = 0

    with open_maybe_gz(fastq) as f:
        for _ in f:
            lines += 1

    return lines // 4


def resolve_reference_file(reference_dir, filename):
    path = Path(filename).expanduser()

    if path.is_absolute():
        return path

    return Path(reference_dir).expanduser() / filename


def bowtie2_index_exists(index_prefix):
    index_prefix = str(index_prefix)

    small_suffixes = [
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ]

    large_suffixes = [
        ".1.bt2l",
        ".2.bt2l",
        ".3.bt2l",
        ".4.bt2l",
        ".rev.1.bt2l",
        ".rev.2.bt2l",
    ]

    small_exists = all(Path(index_prefix + s).exists() for s in small_suffixes)
    large_exists = all(Path(index_prefix + s).exists() for s in large_suffixes)

    return small_exists or large_exists


def build_bowtie2_index_if_needed(fasta, index_prefix, threads):
    if bowtie2_index_exists(index_prefix):
        print(f"Bowtie2 index exists: {index_prefix}")
        return

    print(f"Bowtie2 index not found. Building: {index_prefix}")

    run([
        "bowtie2-build",
        "--threads", str(threads),
        str(fasta),
        str(index_prefix),
    ])


def map_reads_with_bowtie2(
    fastq,
    index_prefix,
    sam_out,
    log_out,
    threads,
    mode,
    extra_args,
    overwrite
):
    if sam_out.exists() and not overwrite:
        print(f"SAM already exists. Skipping Bowtie2: {sam_out}")
        return

    cmd = [
        "bowtie2",
        mode,
        "-p", str(threads),
        "-x", str(index_prefix),
        "-U", str(fastq),
        "-S", str(sam_out),
    ]

    if extra_args:
        cmd[1:1] = extra_args

    run(cmd, stderr_path=log_out)


def run_featurecounts(
    sam_file,
    gff,
    output_txt,
    log_out,
    threads,
    stranded,
    feature_type,
    gene_attribute,
    allow_multi_overlap,
    overwrite
):
    if output_txt.exists() and not overwrite:
        print(f"featureCounts output exists. Skipping: {output_txt}")
        return

    cmd = [
        "featureCounts",
        "-T", str(threads),
        "-a", str(gff),
        "-s", str(stranded),
        "-t", str(feature_type),
        "-g", str(gene_attribute),
        "-o", str(output_txt),
        str(sam_file),
    ]

    if allow_multi_overlap:
        cmd.insert(1, "-O")

    run(cmd, stderr_path=log_out)


def parse_featurecounts_counts(featurecounts_txt, sample, reference_name, total_input_reads):
    df = pd.read_csv(featurecounts_txt, sep="\t", comment="#")

    annotation_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    count_cols = [c for c in df.columns if c not in annotation_cols]

    if len(count_cols) != 1:
        raise ValueError(f"Expected one count column in {featurecounts_txt}, found: {count_cols}")

    count_col = count_cols[0]

    df = df.rename(columns={
        "Geneid": "gene",
        count_col: "mapped_reads"
    })

    total_gene_counts = int(df["mapped_reads"].sum())

    df.insert(0, "reference", reference_name)
    df.insert(0, "sample", sample)

    df["total_input_reads"] = total_input_reads
    df["total_gene_assigned_counts"] = total_gene_counts

    if total_input_reads > 0:
        df["percent_of_input_reads"] = (df["mapped_reads"] / total_input_reads) * 100
    else:
        df["percent_of_input_reads"] = 0.0

    if total_gene_counts > 0:
        df["percent_of_gene_assigned_counts"] = (df["mapped_reads"] / total_gene_counts) * 100
    else:
        df["percent_of_gene_assigned_counts"] = 0.0

    df = df.sort_values("mapped_reads", ascending=False)

    return df


def parse_featurecounts_summary(summary_file, sample, reference_name, total_input_reads):
    if not Path(summary_file).exists():
        return pd.DataFrame()

    df = pd.read_csv(summary_file, sep="\t")

    count_col = df.columns[-1]

    df = df.rename(columns={
        "Status": "status",
        count_col: "read_count"
    })

    df.insert(0, "reference", reference_name)
    df.insert(0, "sample", sample)
    df["total_input_reads"] = total_input_reads

    if total_input_reads > 0:
        df["percent_of_input_reads"] = (df["read_count"] / total_input_reads) * 100
    else:
        df["percent_of_input_reads"] = 0.0

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Map leftover FASTQ reads to reference FASTA/GFF pairs and count gene-level matches."
    )
    parser.add_argument("--config", required=True, help="YAML config file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    check_programs()

    fastq_glob = config["input"]["fastq_glob"]
    output_dir = Path(config["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    fastqs = sorted(Path(p).expanduser() for p in glob.glob(fastq_glob))

    if not fastqs:
        raise FileNotFoundError(f"No FASTQ files found with glob: {fastq_glob}")

    reference_by_name = {ref["name"]: ref for ref in config["references"]}
    match_csv = find_sample_ref_match_csv(config)
    sample_reference_matches = {}
    if match_csv:
        if not match_csv.exists():
            raise FileNotFoundError(f"sample_ref_match CSV not found: {match_csv}")
        sample_reference_matches = load_sample_reference_matches(match_csv)
        unknown_refs = sorted(set(sample_reference_matches.values()) - set(reference_by_name))
        if unknown_refs:
            raise ValueError(
                "sample_ref_match CSV contains references not defined in config: "
                + ", ".join(unknown_refs)
            )
        unmatched_fastqs = [
            str(fastq) for fastq in fastqs
            if reference_for_sample(fastq, sample_name_from_fastq(fastq), sample_reference_matches) is None
        ]
        if unmatched_fastqs:
            raise ValueError(
                "No sample_ref_match assignment found for FASTQ(s): "
                + ", ".join(unmatched_fastqs)
            )
        print(f"Loaded sample/reference assignments from: {match_csv}")
    else:
        print("No sample_ref_match CSV configured/found; all FASTQs will be mapped to all references.")

    bowtie_cfg = config.get("bowtie2", {})
    bowtie_threads = bowtie_cfg.get("threads", 16)
    bowtie_mode = bowtie_cfg.get("mode", "--end-to-end")
    bowtie_extra_args = bowtie_cfg.get("extra_args", [])

    fc_cfg = config.get("featurecounts", {})
    fc_threads = fc_cfg.get("threads", 16)
    fc_stranded = fc_cfg.get("stranded", 0)
    fc_allow_multi_overlap = fc_cfg.get("allow_multi_overlap", False)

    settings = config.get("settings", {})
    overwrite = settings.get("overwrite_existing", False)
    keep_sam = settings.get("keep_sam", True)
    continue_on_error = settings.get("continue_on_error", False)

    print(f"Found {len(fastqs)} FASTQ files.")

    all_gene_csvs = []
    all_summary_csvs = []
    manifest_rows = []

    for ref_name, ref in reference_by_name.items():
        ref_name = ref["name"]
        ref_fastqs = [
            fastq for fastq in fastqs
            if not sample_reference_matches
            or reference_for_sample(fastq, sample_name_from_fastq(fastq), sample_reference_matches) == ref_name
        ]

        if not ref_fastqs:
            print(f"No FASTQs assigned to reference {ref_name}; skipping.")
            continue

        ref_dir = Path(ref["reference_dir"]).expanduser()

        fasta = resolve_reference_file(ref_dir, ref["fasta"])
        gff = resolve_reference_file(ref_dir, ref["gff"])

        if not fasta.exists():
            raise FileNotFoundError(f"Missing FASTA for {ref_name}: {fasta}")

        if not gff.exists():
            raise FileNotFoundError(f"Missing GFF for {ref_name}: {gff}")

        index_prefix_name = ref.get("index_prefix", fasta.stem)
        index_prefix = resolve_reference_file(ref_dir, index_prefix_name)

        feature_type = ref.get("feature_type", fc_cfg.get("feature_type", "gene"))
        gene_attribute = ref.get("gene_attribute", fc_cfg.get("gene_attribute", "Name"))

        print("========================================")
        print(f"Reference: {ref_name}")
        print(f"FASTA: {fasta}")
        print(f"GFF:   {gff}")
        print(f"Index: {index_prefix}")
        print("========================================")

        build_bowtie2_index_if_needed(
            fasta=fasta,
            index_prefix=index_prefix,
            threads=bowtie_threads
        )

        for fastq in ref_fastqs:
            sample = sample_name_from_fastq(fastq)

            sample_ref_dir = output_dir / ref_name / sample
            sample_ref_dir.mkdir(parents=True, exist_ok=True)

            sam_out = sample_ref_dir / f"{sample}.{ref_name}.sam"
            bowtie_log = sample_ref_dir / f"{sample}.{ref_name}.bowtie2.log"

            featurecounts_txt = sample_ref_dir / f"{sample}.{ref_name}.featureCounts.txt"
            featurecounts_log = sample_ref_dir / f"{sample}.{ref_name}.featureCounts.log"
            featurecounts_summary = Path(str(featurecounts_txt) + ".summary")

            gene_csv = sample_ref_dir / f"{sample}.{ref_name}.gene_counts.csv"
            summary_csv = sample_ref_dir / f"{sample}.{ref_name}.featureCounts_summary.csv"

            try:
                print("========================================")
                print(f"Sample: {sample}")
                print(f"FASTQ:  {fastq}")
                print(f"Ref:    {ref_name}")
                print("========================================")

                total_input_reads = count_fastq_reads(fastq)
                print(f"Total input FASTQ reads: {total_input_reads}")

                map_reads_with_bowtie2(
                    fastq=fastq,
                    index_prefix=index_prefix,
                    sam_out=sam_out,
                    log_out=bowtie_log,
                    threads=bowtie_threads,
                    mode=bowtie_mode,
                    extra_args=bowtie_extra_args,
                    overwrite=overwrite
                )

                run_featurecounts(
                    sam_file=sam_out,
                    gff=gff,
                    output_txt=featurecounts_txt,
                    log_out=featurecounts_log,
                    threads=fc_threads,
                    stranded=fc_stranded,
                    feature_type=feature_type,
                    gene_attribute=gene_attribute,
                    allow_multi_overlap=fc_allow_multi_overlap,
                    overwrite=overwrite
                )

                gene_df = parse_featurecounts_counts(
                    featurecounts_txt=featurecounts_txt,
                    sample=sample,
                    reference_name=ref_name,
                    total_input_reads=total_input_reads
                )

                gene_df.to_csv(gene_csv, index=False)
                all_gene_csvs.append(gene_csv)
                print(f"Wrote: {gene_csv}")

                summary_df = parse_featurecounts_summary(
                    summary_file=featurecounts_summary,
                    sample=sample,
                    reference_name=ref_name,
                    total_input_reads=total_input_reads
                )

                if not summary_df.empty:
                    summary_df.to_csv(summary_csv, index=False)
                    all_summary_csvs.append(summary_csv)
                    print(f"Wrote: {summary_csv}")

                if not keep_sam:
                    sam_out.unlink(missing_ok=True)

                manifest_rows.append({
                    "sample": sample,
                    "reference": ref_name,
                    "fastq": str(fastq),
                    "fasta": str(fasta),
                    "gff": str(gff),
                    "index_prefix": str(index_prefix),
                    "sam": str(sam_out) if keep_sam else "",
                    "bowtie2_log": str(bowtie_log),
                    "featurecounts_txt": str(featurecounts_txt),
                    "featurecounts_summary": str(featurecounts_summary),
                    "gene_counts_csv": str(gene_csv),
                    "status": "done"
                })

            except Exception as e:
                manifest_rows.append({
                    "sample": sample,
                    "reference": ref_name,
                    "fastq": str(fastq),
                    "fasta": str(fasta),
                    "gff": str(gff),
                    "index_prefix": str(index_prefix),
                    "sam": str(sam_out),
                    "bowtie2_log": str(bowtie_log),
                    "featurecounts_txt": str(featurecounts_txt),
                    "featurecounts_summary": str(featurecounts_summary),
                    "gene_counts_csv": str(gene_csv),
                    "status": f"failed: {e}"
                })

                if continue_on_error:
                    print(f"WARNING: failed {sample} vs {ref_name}: {e}")
                    continue

                raise

    manifest = pd.DataFrame(manifest_rows)
    manifest_csv = output_dir / "gene_mapping_manifest.csv"
    manifest.to_csv(manifest_csv, index=False)
    print(f"Wrote: {manifest_csv}")

    if all_gene_csvs:
        combined_gene_df = pd.concat(
            [pd.read_csv(p) for p in all_gene_csvs],
            ignore_index=True
        )

        combined_gene_csv = output_dir / "combined_gene_counts.csv"
        combined_gene_df.to_csv(combined_gene_csv, index=False)
        print(f"Wrote: {combined_gene_csv}")

    if all_summary_csvs:
        combined_summary_df = pd.concat(
            [pd.read_csv(p) for p in all_summary_csvs],
            ignore_index=True
        )

        combined_summary_csv = output_dir / "combined_featureCounts_summary.csv"
        combined_summary_df.to_csv(combined_summary_csv, index=False)
        print(f"Wrote: {combined_summary_csv}")

    print("Done.")


if __name__ == "__main__":
    main()
