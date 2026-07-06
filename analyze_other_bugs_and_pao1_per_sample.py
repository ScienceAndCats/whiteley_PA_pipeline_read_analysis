#!/usr/bin/env python3

import argparse
import gzip
import subprocess
from pathlib import Path
from collections import Counter

import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


def open_maybe_gz(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def find_one_file(folder, patterns, label):
    hits = []
    folder = Path(folder)

    for pattern in patterns:
        hits.extend(folder.glob(pattern))

    if not hits:
        raise FileNotFoundError(f"No {label} found in {folder} using {patterns}")

    if len(hits) > 1:
        print(f"WARNING: multiple {label} files found. Using: {hits[0]}")

    return hits[0]


def get_sample_name_from_decoy_sam(path):
    name = Path(path).name
    suffix = ".mapped_to_other_bugs.sam"

    if name.endswith(suffix):
        return name.replace(suffix, "")

    return Path(path).stem


def parse_decoy_sam(sam_file, skip_secondary=True):
    counts = Counter()

    with open_maybe_gz(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue

            flag = int(fields[1])
            ref_name = fields[2]

            # Skip unmapped reads
            if flag & 4 or ref_name == "*":
                continue

            # Skip secondary/supplementary alignments
            if skip_secondary and (flag & 256 or flag & 2048):
                continue

            counts[ref_name] += 1

    df = pd.DataFrame(
        [{"reference_name": ref, "mapped_reads": count}
         for ref, count in counts.items()]
    )

    if not df.empty:
        df = df.sort_values("mapped_reads", ascending=False)

    return df


def plot_top_bar(df, name_col, count_col, title, xlabel, ylabel, out_png, top_n):
    if df.empty:
        print(f"Skipping empty plot: {out_png}")
        return

    top = df.sort_values(count_col, ascending=False).head(top_n)

    plt.figure(figsize=(10, max(6, len(top) * 0.25)))
    plt.barh(top[name_col][::-1], top[count_col][::-1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"Wrote {out_png}")


def build_bowtie2_index(fasta, index_prefix, threads):
    index_prefix = Path(index_prefix)

    possible_index_files = [
        Path(str(index_prefix) + ".1.bt2"),
        Path(str(index_prefix) + ".1.bt2l")
    ]

    if any(p.exists() for p in possible_index_files):
        print("Bowtie2 index already exists. Skipping index build.")
        return

    run([
        "bowtie2-build",
        "--threads", str(threads),
        str(fasta),
        str(index_prefix)
    ])


def map_fastq_to_pao1(sample, fastq, index_prefix, sample_dir, threads, overwrite=False):
    sam_out = sample_dir / f"{sample}.PAO1.sam"
    log_out = sample_dir / f"{sample}.PAO1.bowtie2.log"

    if sam_out.exists() and not overwrite:
        print(f"PAO1 SAM already exists. Skipping: {sam_out}")
        return sam_out

    run([
        "bowtie2",
        "--end-to-end",
        "-p", str(threads),
        "-x", str(index_prefix),
        "-q",
        "-U", str(fastq),
        "-S", str(sam_out)
    ], stderr_path=log_out)

    return sam_out


def run_featurecounts_single_sample(
    sam_file,
    gff,
    out_txt,
    threads,
    stranded,
    feature_type,
    gene_attribute,
    overwrite=False
):
    if Path(out_txt).exists() and not overwrite:
        print(f"featureCounts output already exists. Skipping: {out_txt}")
        return

    run([
        "featureCounts",
        "-T", str(threads),
        "-a", str(gff),
        "-O",
        "-s", str(stranded),
        "-g", str(gene_attribute),
        "-t", str(feature_type),
        "-o", str(out_txt),
        str(sam_file)
    ])


def parse_featurecounts_single_sample(featurecounts_txt):
    df = pd.read_csv(featurecounts_txt, sep="\t", comment="#")

    annotation_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    count_cols = [c for c in df.columns if c not in annotation_cols]

    if len(count_cols) != 1:
        raise ValueError(f"Expected one sample count column, found: {count_cols}")

    count_col = count_cols[0]

    df = df.rename(columns={
        "Geneid": "gene",
        count_col: "mapped_reads"
    })

    df = df.sort_values("mapped_reads", ascending=False)

    return df


def find_matching_unmapped_fastq(input_dir, sample):
    input_dir = Path(input_dir)

    candidates = [
        input_dir / f"{sample}_unmapped_to_other_bugs.fastq.gz",
        input_dir / f"{sample}_unmapped_to_other_bugs.fastq",
        input_dir / f"{sample}_unmapped_to_other_bugs.fq.gz",
        input_dir / f"{sample}_unmapped_to_other_bugs.fq",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate

    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="YAML config file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    input_dir = Path(config["input_dir"]).expanduser()
    output_dir = Path(config["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    threads = config["settings"].get("threads", 16)
    top_n = config["settings"].get("top_n", 30)
    skip_secondary = config["settings"].get("skip_secondary_alignments", True)
    overwrite = config["settings"].get("overwrite_existing", False)

    decoy_pattern = config["file_patterns"]["decoy_sam"]

    pao1_dir = Path(config["pao1"]["reference_dir"]).expanduser()
    fasta = find_one_file(pao1_dir, config["pao1"]["fasta_extensions"], "PAO1 FASTA")
    gff = find_one_file(pao1_dir, config["pao1"]["gff_extensions"], "PAO1 GFF")

    print(f"Using PAO1 FASTA: {fasta}")
    print(f"Using PAO1 GFF:   {gff}")

    # Build PAO1 Bowtie2 index once, shared by all samples
    index_dir = output_dir / "_shared_PAO1_bowtie2_index"
    index_dir.mkdir(exist_ok=True)
    index_prefix = index_dir / "PAO1"

    build_bowtie2_index(
        fasta=fasta,
        index_prefix=index_prefix,
        threads=threads
    )

    decoy_sams = sorted(input_dir.glob(decoy_pattern))

    if not decoy_sams:
        raise FileNotFoundError(f"No decoy SAM files found in {input_dir} using {decoy_pattern}")

    manifest_rows = []

    for decoy_sam in decoy_sams:
        sample = get_sample_name_from_decoy_sam(decoy_sam)
        sample_dir = output_dir / sample
        sample_dir.mkdir(parents=True, exist_ok=True)

        print("========================================")
        print(f"Processing sample: {sample}")
        print("========================================")

        # ----------------------------
        # 1. Decoy/other-bug counts
        # ----------------------------
        decoy_df = parse_decoy_sam(
            sam_file=decoy_sam,
            skip_secondary=skip_secondary
        )

        decoy_csv = sample_dir / "decoy_mapped_reference_counts.csv"
        decoy_df.to_csv(decoy_csv, index=False)
        print(f"Wrote {decoy_csv}")

        plot_top_bar(
            df=decoy_df,
            name_col="reference_name",
            count_col="mapped_reads",
            title=f"{sample}: top decoy references",
            xlabel="Mapped reads",
            ylabel="Decoy reference",
            out_png=sample_dir / "decoy_top_references_barplot.png",
            top_n=top_n
        )

        # ----------------------------
        # 2. Find matching unmapped FASTQ
        # ----------------------------
        unmapped_fastq = find_matching_unmapped_fastq(input_dir, sample)

        if unmapped_fastq is None:
            print(f"WARNING: no matching unmapped FASTQ found for {sample}. Skipping PAO1 mapping.")

            manifest_rows.append({
                "sample": sample,
                "decoy_sam": str(decoy_sam),
                "unmapped_fastq": "",
                "sample_output_dir": str(sample_dir),
                "status": "decoy_done_no_unmapped_fastq"
            })
            continue

        # ----------------------------
        # 3. Map this sample to PAO1
        # ----------------------------
        pao1_sam = map_fastq_to_pao1(
            sample=sample,
            fastq=unmapped_fastq,
            index_prefix=index_prefix,
            sample_dir=sample_dir,
            threads=threads,
            overwrite=overwrite
        )

        # ----------------------------
        # 4. featureCounts for this sample only
        # ----------------------------
        featurecounts_txt = sample_dir / "PAO1_featureCounts.txt"

        run_featurecounts_single_sample(
            sam_file=pao1_sam,
            gff=gff,
            out_txt=featurecounts_txt,
            threads=threads,
            stranded=config["settings"].get("featurecounts_stranded", 1),
            feature_type=config["settings"].get("featurecounts_feature_type", "CDS"),
            gene_attribute=config["settings"].get("featurecounts_gene_attribute", "locus"),
            overwrite=overwrite
        )

        # ----------------------------
        # 5. PAO1 gene counts CSV + graph
        # ----------------------------
        pao1_df = parse_featurecounts_single_sample(featurecounts_txt)

        pao1_csv = sample_dir / "PAO1_gene_counts.csv"
        pao1_df.to_csv(pao1_csv, index=False)
        print(f"Wrote {pao1_csv}")

        plot_top_bar(
            df=pao1_df[pao1_df["mapped_reads"] > 0],
            name_col="gene",
            count_col="mapped_reads",
            title=f"{sample}: top PAO1 genes",
            xlabel="Mapped reads",
            ylabel="PAO1 gene/locus",
            out_png=sample_dir / "PAO1_top_genes_barplot.png",
            top_n=top_n
        )

        manifest_rows.append({
            "sample": sample,
            "decoy_sam": str(decoy_sam),
            "unmapped_fastq": str(unmapped_fastq),
            "sample_output_dir": str(sample_dir),
            "decoy_counts_csv": str(decoy_csv),
            "pao1_sam": str(pao1_sam),
            "pao1_featurecounts": str(featurecounts_txt),
            "pao1_gene_counts_csv": str(pao1_csv),
            "status": "done"
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest_csv = output_dir / "sample_manifest.csv"
    manifest.to_csv(manifest_csv, index=False)
    print(f"Wrote {manifest_csv}")

    print("Done.")


if __name__ == "__main__":
    main()