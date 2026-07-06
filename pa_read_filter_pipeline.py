#!/usr/bin/env python3
"""Combined PA read-analysis pipeline.

This script replaces the separate PAO1/decoy analysis, PAO1-unmapped read
extraction, and sourmash leftover-classification scripts.  It writes per-sample
outputs plus a combined read-count report and Sankey diagrams.
"""

import argparse
import gzip
import shutil
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml


def run(cmd, stdout_path=None, stderr_path=None):
    print("Running:", " ".join(map(str, cmd)), flush=True)
    stdout = open(stdout_path, "w") if stdout_path else None
    stderr = open(stderr_path, "w") if stderr_path else None
    try:
        subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)
    finally:
        if stdout:
            stdout.close()
        if stderr:
            stderr.close()


def open_text_maybe_gz(path, mode="rt"):
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def count_fastq_reads(path):
    if not path or not Path(path).exists():
        return 0
    with open_text_maybe_gz(path, "rt") as handle:
        return sum(1 for _ in handle) // 4


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
    return name.replace(suffix, "") if name.endswith(suffix) else Path(path).stem


def parse_decoy_sam(sam_file, skip_secondary=True):
    counts = Counter()
    total = mapped = skipped_secondary_count = 0
    with open_text_maybe_gz(sam_file) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            total += 1
            flag = int(fields[1])
            ref_name = fields[2]
            if skip_secondary and (flag & 256 or flag & 2048):
                skipped_secondary_count += 1
                continue
            if flag & 4 or ref_name == "*":
                continue
            mapped += 1
            counts[ref_name] += 1
    df = pd.DataFrame([{"reference_name": r, "mapped_reads": c} for r, c in counts.items()])
    if not df.empty:
        df = df.sort_values("mapped_reads", ascending=False)
    return df, {"decoy_sam_records": total, "decoy_mapped_reads": mapped, "decoy_skipped_secondary": skipped_secondary_count}


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


def build_bowtie2_index(fasta, index_prefix, threads):
    if any(Path(str(index_prefix) + ext).exists() for ext in [".1.bt2", ".1.bt2l"]):
        print("Bowtie2 index already exists. Skipping index build.")
        return
    run(["bowtie2-build", "--threads", str(threads), str(fasta), str(index_prefix)])


def map_fastq_to_pao1(sample, fastq, index_prefix, sample_dir, threads, overwrite=False):
    sam_out = sample_dir / f"{sample}.PAO1.sam"
    log_out = sample_dir / f"{sample}.PAO1.bowtie2.log"
    if sam_out.exists() and not overwrite:
        return sam_out
    run(["bowtie2", "--end-to-end", "-p", str(threads), "-x", str(index_prefix), "-q", "-U", str(fastq), "-S", str(sam_out)], stderr_path=log_out)
    return sam_out


def run_featurecounts_single_sample(sam_file, gff, out_txt, threads, stranded, feature_type, gene_attribute, overwrite=False):
    if Path(out_txt).exists() and not overwrite:
        return
    run(["featureCounts", "-T", str(threads), "-a", str(gff), "-O", "-s", str(stranded), "-g", str(gene_attribute), "-t", str(feature_type), "-o", str(out_txt), str(sam_file)])


def parse_featurecounts_single_sample(featurecounts_txt):
    df = pd.read_csv(featurecounts_txt, sep="\t", comment="#")
    annotation_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    count_cols = [c for c in df.columns if c not in annotation_cols]
    if len(count_cols) != 1:
        raise ValueError(f"Expected one sample count column, found: {count_cols}")
    return df.rename(columns={"Geneid": "gene", count_cols[0]: "mapped_reads"}).sort_values("mapped_reads", ascending=False)


def find_matching_unmapped_fastq(input_dir, sample):
    for suffix in ["_unmapped_to_other_bugs.fastq.gz", "_unmapped_to_other_bugs.fastq", "_unmapped_to_other_bugs.fq.gz", "_unmapped_to_other_bugs.fq"]:
        candidate = Path(input_dir) / f"{sample}{suffix}"
        if candidate.exists():
            return candidate
    return None


def sam_unmapped_to_fastq(pao1_sam, output_fastq, gzip_output=True):
    opener = gzip.open if gzip_output else open
    total = mapped = unmapped = written = skipped = 0
    with open_text_maybe_gz(pao1_sam, "rt") as fin, opener(output_fastq, "wt") as fout:
        for line in fin:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            total += 1
            flag = int(fields[1])
            if flag & 4:
                unmapped += 1
                if fields[9] == "*" or fields[10] == "*":
                    skipped += 1
                    continue
                fout.write(f"@{fields[0]}\n{fields[9]}\n+\n{fields[10]}\n")
                written += 1
            else:
                mapped += 1
    return {"pao1_sam_records": total, "pao1_mapped_reads": mapped, "pao1_unmapped_reads": unmapped, "leftover_fastq_reads": written, "pao1_unmapped_missing_seq_or_qual": skipped}


def sourmash_param_string(ksize, scaled, abundance):
    return f"k={ksize},scaled={scaled}" + (",abund" if abundance else "")


def sample_name_from_fastq(fastq_path):
    name = Path(fastq_path).name
    for suffix in [".decoy_filtered.PAO1_unmapped.fastq.gz", ".decoy_filtered.PAO1_unmapped.fastq", ".fastq.gz", ".fastq", ".fq.gz", ".fq"]:
        if name.endswith(suffix):
            return name.replace(suffix, "")
    return Path(fastq_path).stem


def run_sourmash_for_sample(fastq, database_cfg, overwrite):
    if shutil.which("sourmash") is None:
        raise RuntimeError("Could not find sourmash in PATH.")
    fastq = Path(fastq)
    sample = sample_name_from_fastq(fastq)
    db_name = database_cfg["name"]
    database = Path(database_cfg["database"]).expanduser()
    ksize = database_cfg.get("ksize", 51)
    scaled = database_cfg.get("scaled", 1000)
    abundance = database_cfg.get("abundance", True)
    if not database.exists():
        raise FileNotFoundError(f"Sourmash database not found: {database}")
    sig = fastq.parent / f"{sample}.{db_name}.sourmash.k{ksize}.sig.zip"
    gather_csv = fastq.parent / f"{sample}.{db_name}.sourmash_gather.k{ksize}.csv"
    if not sig.exists() or overwrite:
        run(["sourmash", "sketch", "dna", "-p", sourmash_param_string(ksize, scaled, abundance), "--name", sample, "-o", str(sig), str(fastq)], stderr_path=fastq.parent / f"{sample}.{db_name}.sourmash_sketch.k{ksize}.log")
    if not gather_csv.exists() or overwrite:
        run(["sourmash", "gather", str(sig), str(database), "-k", str(ksize), "-o", str(gather_csv)], stderr_path=fastq.parent / f"{sample}.{db_name}.sourmash_gather.k{ksize}.log")
    matched_fraction = 0.0
    if gather_csv.exists():
        try:
            df = pd.read_csv(gather_csv)
            for col in ["f_unique_to_query", "f_match", "intersect_bp"]:
                if col in df.columns and col != "intersect_bp":
                    matched_fraction = min(1.0, float(df[col].fillna(0).sum()))
                    break
        except pd.errors.EmptyDataError:
            pass
    return {"database": db_name, "sourmash_gather_csv": str(gather_csv), "sourmash_matched_fraction": matched_fraction}


def write_sankey_html(sample, links, out_html):
    labels = []
    label_to_index = {}
    sources = []
    targets = []
    values = []
    for source, target, value in links:
        for label in [source, target]:
            if label not in label_to_index:
                label_to_index[label] = len(labels)
                labels.append(label)
        sources.append(label_to_index[source])
        targets.append(label_to_index[target])
        values.append(max(0, int(value)))
    try:
        import plotly.graph_objects as go
        fig = go.Figure(data=[go.Sankey(node={"label": labels}, link={"source": sources, "target": targets, "value": values})])
        fig.update_layout(title_text=f"{sample}: read filtering flow", font_size=10)
        fig.write_html(out_html)
    except Exception as exc:
        pd.DataFrame(links, columns=["source", "target", "reads"]).to_html(out_html, index=False)
        print(f"WARNING: wrote table fallback for Sankey ({exc}): {out_html}")


def summarize_pao1_sam(pao1_sam):
    total = mapped = unmapped = skipped = 0
    with open_text_maybe_gz(pao1_sam, "rt") as fin:
        for line in fin:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            total += 1
            flag = int(fields[1])
            if flag & 4:
                unmapped += 1
                if fields[9] == "*" or fields[10] == "*":
                    skipped += 1
            else:
                mapped += 1
    return {"pao1_sam_records": total, "pao1_mapped_reads": mapped, "pao1_unmapped_reads": unmapped, "pao1_unmapped_missing_seq_or_qual": skipped}


def process_one_sample(decoy_sam, config, fasta, gff, index_prefix):
    input_dir = Path(config["input_dir"]).expanduser()
    output_dir = Path(config["output_dir"]).expanduser()
    settings = config.get("settings", {})
    sample = get_sample_name_from_decoy_sam(decoy_sam)
    sample_dir = output_dir / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    top_n = settings.get("top_n", 30)
    overwrite = settings.get("overwrite_existing", False)

    decoy_df, stats = parse_decoy_sam(decoy_sam, settings.get("skip_secondary_alignments", True))
    decoy_csv = sample_dir / "decoy_mapped_reference_counts.csv"
    decoy_df.to_csv(decoy_csv, index=False)
    plot_top_bar(decoy_df, "reference_name", "mapped_reads", f"{sample}: top decoy references", "Mapped reads", "Decoy reference", sample_dir / "decoy_top_references_barplot.png", top_n)

    unmapped_fastq = find_matching_unmapped_fastq(input_dir, sample)
    started_reads = stats["decoy_mapped_reads"] + count_fastq_reads(unmapped_fastq)
    stats.update({"sample": sample, "started_reads": started_reads, "input_decoy_sam": str(decoy_sam), "decoy_unmapped_fastq": str(unmapped_fastq or "")})
    if not unmapped_fastq:
        stats["status"] = "decoy_done_no_unmapped_fastq"
        return stats

    pao1_sam = map_fastq_to_pao1(sample, unmapped_fastq, index_prefix, sample_dir, settings.get("threads", 16), overwrite)
    featurecounts_txt = sample_dir / "PAO1_featureCounts.txt"
    run_featurecounts_single_sample(pao1_sam, gff, featurecounts_txt, settings.get("threads", 16), settings.get("featurecounts_stranded", 1), settings.get("featurecounts_feature_type", "CDS"), settings.get("featurecounts_gene_attribute", "locus"), overwrite)
    pao1_df = parse_featurecounts_single_sample(featurecounts_txt)
    pao1_df.to_csv(sample_dir / "PAO1_gene_counts.csv", index=False)
    plot_top_bar(pao1_df[pao1_df["mapped_reads"] > 0], "gene", "mapped_reads", f"{sample}: top PAO1 genes", "Mapped reads", "PAO1 gene/locus", sample_dir / "PAO1_top_genes_barplot.png", top_n)

    output_fastq = sample_dir / f"{sample}{settings.get('output_filename_suffix', '.decoy_filtered.PAO1_unmapped.fastq.gz')}"
    if output_fastq.exists() and not overwrite:
        extract_stats = summarize_pao1_sam(pao1_sam)
        extract_stats["leftover_fastq_reads"] = count_fastq_reads(output_fastq)
    else:
        extract_stats = sam_unmapped_to_fastq(pao1_sam, output_fastq, settings.get("gzip_output", True))
    stats.update(extract_stats)
    stats.update({"pao1_sam": str(pao1_sam), "leftover_fastq": str(output_fastq), "status": "done"})

    remaining = int(stats.get("leftover_fastq_reads") or 0)
    for db_cfg in config.get("sourmash", {}).get("databases", []):
        sm = run_sourmash_for_sample(output_fastq, db_cfg, overwrite)
        est = int(round(remaining * sm["sourmash_matched_fraction"]))
        stats[f"{sm['database']}_sourmash_estimated_matched_reads"] = est
        stats[f"{sm['database']}_sourmash_unmatched_reads"] = max(0, remaining - est)
        stats[f"{sm['database']}_sourmash_gather_csv"] = sm["sourmash_gather_csv"]
        remaining = max(0, remaining - est)
    stats["final_unmatched_unfiltered_reads"] = remaining

    not_decoy = count_fastq_reads(unmapped_fastq)
    not_pao1 = int(stats.get("leftover_fastq_reads") or 0)
    links = [
        ("started", "aligned to decoy/other bugs", int(stats.get("decoy_mapped_reads") or 0)),
        ("started", "not decoy", not_decoy),
        ("not decoy", "aligned to PAO1", int(stats.get("pao1_mapped_reads") or 0)),
        ("not decoy", "not PAO1", not_pao1),
    ]
    previous = "not PAO1"
    for db_cfg in config.get("sourmash", {}).get("databases", []):
        name = db_cfg["name"]
        matched = int(stats.get(f"{name}_sourmash_estimated_matched_reads", 0) or 0)
        unmatched = int(stats.get(f"{name}_sourmash_unmatched_reads", 0) or 0)
        links.append((previous, f"matched {name}", matched))
        links.append((previous, f"not {name}", unmatched))
        previous = f"not {name}"
    write_sankey_html(sample, links, sample_dir / "read_filter_sankey.html")
    return stats


def main():
    parser = argparse.ArgumentParser(description="Combined decoy, PAO1, extraction, sourmash, report, and Sankey pipeline.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)

    output_dir = Path(config["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    settings = config.get("settings", {})
    threads = int(settings.get("threads", 16))
    sample_workers = int(settings.get("max_parallel_samples", 1))
    continue_on_error = bool(settings.get("continue_on_error", False))

    pao1_dir = Path(config["pao1"]["reference_dir"]).expanduser()
    fasta = find_one_file(pao1_dir, config["pao1"]["fasta_extensions"], "PAO1 FASTA")
    gff = find_one_file(pao1_dir, config["pao1"]["gff_extensions"], "PAO1 GFF")
    index_dir = output_dir / "_shared_PAO1_bowtie2_index"
    index_dir.mkdir(exist_ok=True)
    index_prefix = index_dir / "PAO1"
    build_bowtie2_index(fasta, index_prefix, threads)

    decoy_sams = sorted(Path(config["input_dir"]).expanduser().glob(config["file_patterns"]["decoy_sam"]))
    if not decoy_sams:
        raise FileNotFoundError("No decoy SAM files found")

    rows = []
    if sample_workers == 1:
        for sam in decoy_sams:
            rows.append(process_one_sample(sam, config, fasta, gff, index_prefix))
    else:
        with ThreadPoolExecutor(max_workers=sample_workers) as executor:
            futures = {executor.submit(process_one_sample, sam, config, fasta, gff, index_prefix): sam for sam in decoy_sams}
            for future in as_completed(futures):
                try:
                    rows.append(future.result())
                except Exception as exc:
                    if not continue_on_error:
                        raise
                    rows.append({"sample": get_sample_name_from_decoy_sam(futures[future]), "status": "failed", "error": repr(exc)})

    report = pd.DataFrame(rows).sort_values("sample")
    report_csv = output_dir / settings.get("read_report_csv", "read_filter_report.csv")
    report.to_csv(report_csv, index=False)
    report.to_csv(output_dir / "sample_manifest.csv", index=False)
    print(f"Wrote read report: {report_csv}")


if __name__ == "__main__":
    main()
