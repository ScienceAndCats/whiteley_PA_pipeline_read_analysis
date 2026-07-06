#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import yaml
import pandas as pd


def open_text_maybe_gz(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def open_output_fastq(path, gzip_output=True):
    path = Path(path)

    if gzip_output:
        return gzip.open(path, "wt")

    return open(path, "w")


def sam_unmapped_to_fastq(pao1_sam, output_fastq, gzip_output=True):
    """
    Extract reads from a PAO1 SAM file where SAM flag 4 is set.

    These are reads that:
      1. already passed the decoy-filtering step
      2. failed to map to PAO1

    Writes them as FASTQ.
    """
    total_records = 0
    unmapped_records = 0
    written_records = 0
    skipped_missing_seq_or_qual = 0

    with open_text_maybe_gz(pao1_sam, "rt") as fin, open_output_fastq(output_fastq, gzip_output) as fout:
        for line in fin:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 11:
                continue

            total_records += 1

            qname = fields[0]
            flag = int(fields[1])
            seq = fields[9]
            qual = fields[10]

            # SAM flag 4 = read unmapped
            if not (flag & 4):
                continue

            unmapped_records += 1

            if seq == "*" or qual == "*":
                skipped_missing_seq_or_qual += 1
                continue

            fout.write(f"@{qname}\n")
            fout.write(f"{seq}\n")
            fout.write("+\n")
            fout.write(f"{qual}\n")

            written_records += 1

    return {
        "total_sam_records": total_records,
        "pao1_unmapped_sam_records": unmapped_records,
        "fastq_reads_written": written_records,
        "skipped_missing_seq_or_qual": skipped_missing_seq_or_qual,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract decoy-filtered reads that failed to map to PAO1 into per-sample FASTQ files."
    )
    parser.add_argument("--config", required=True, help="YAML config file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    manifest_csv = Path(config["manifest_csv"]).expanduser()
    manifest = pd.read_csv(manifest_csv)

    settings = config.get("settings", {})
    output_suffix = settings.get("output_filename_suffix", ".decoy_filtered.PAO1_unmapped.fastq.gz")
    summary_filename = settings.get("summary_filename", "non_decoy_non_PAO1_fastq_summary.csv")
    overwrite = settings.get("overwrite_existing", False)
    gzip_output = settings.get("gzip_output", True)

    required_cols = ["sample", "pao1_sam", "sample_output_dir"]

    for col in required_cols:
        if col not in manifest.columns:
            raise ValueError(f"Manifest is missing required column: {col}")

    summary_rows = []

    for _, row in manifest.iterrows():
        sample = row["sample"]
        pao1_sam = Path(row["pao1_sam"]).expanduser()
        sample_output_dir = Path(row["sample_output_dir"]).expanduser()

        sample_output_dir.mkdir(parents=True, exist_ok=True)

        output_fastq = sample_output_dir / f"{sample}{output_suffix}"

        print("========================================")
        print(f"Processing sample: {sample}")
        print("========================================")

        if not pao1_sam.exists():
            print(f"WARNING: missing PAO1 SAM: {pao1_sam}")

            summary_rows.append({
                "sample": sample,
                "pao1_sam": str(pao1_sam),
                "output_fastq": str(output_fastq),
                "status": "missing_pao1_sam",
                "total_sam_records": 0,
                "pao1_unmapped_sam_records": 0,
                "fastq_reads_written": 0,
                "skipped_missing_seq_or_qual": 0,
            })
            continue

        if output_fastq.exists() and not overwrite:
            print(f"Output already exists. Skipping: {output_fastq}")

            summary_rows.append({
                "sample": sample,
                "pao1_sam": str(pao1_sam),
                "output_fastq": str(output_fastq),
                "status": "skipped_existing",
                "total_sam_records": "",
                "pao1_unmapped_sam_records": "",
                "fastq_reads_written": "",
                "skipped_missing_seq_or_qual": "",
            })
            continue

        stats = sam_unmapped_to_fastq(
            pao1_sam=pao1_sam,
            output_fastq=output_fastq,
            gzip_output=gzip_output
        )

        print(f"Wrote: {output_fastq}")
        print(f"FASTQ reads written: {stats['fastq_reads_written']}")

        summary_rows.append({
            "sample": sample,
            "pao1_sam": str(pao1_sam),
            "output_fastq": str(output_fastq),
            "status": "done",
            **stats
        })

    summary_df = pd.DataFrame(summary_rows)

    report_root = manifest_csv.parent
    summary_csv = report_root / summary_filename
    summary_df.to_csv(summary_csv, index=False)

    print("========================================")
    print(f"Wrote summary: {summary_csv}")
    print("Done.")


if __name__ == "__main__":
    main()