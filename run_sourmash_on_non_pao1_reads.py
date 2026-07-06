#!/usr/bin/env python3

# Example:
# python run_sourmash_on_non_pao1_reads.py --config sourmash_config.yaml

import argparse
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml
import pandas as pd


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


def check_sourmash_available():
    if shutil.which("sourmash") is None:
        raise RuntimeError(
            "Could not find sourmash in PATH. Install/activate it first, e.g.:\n"
            "conda install -c conda-forge -c bioconda sourmash"
        )


def find_fastqs(input_report_dir, fastq_pattern):
    input_report_dir = Path(input_report_dir)

    # Searches inside per-sample folders
    fastqs = sorted(input_report_dir.glob(f"*/{fastq_pattern}"))

    if not fastqs:
        raise FileNotFoundError(
            f"No FASTQ files found under {input_report_dir} using pattern: */{fastq_pattern}"
        )

    return fastqs


def sample_name_from_fastq(fastq_path):
    name = Path(fastq_path).name

    suffixes = [
        ".decoy_filtered.PAO1_unmapped.fastq.gz",
        ".decoy_filtered.PAO1_unmapped.fastq",
        ".fastq.gz",
        ".fastq",
        ".fq.gz",
        ".fq",
    ]

    for suffix in suffixes:
        if name.endswith(suffix):
            return name.replace(suffix, "")

    return Path(fastq_path).stem


def sourmash_param_string(ksize, scaled, abundance):
    params = f"k={ksize},scaled={scaled}"

    if abundance:
        params += ",abund"

    return params


def run_sourmash_sketch(fastq, sig_out, log_out, sample, ksize, scaled, abundance, overwrite):
    if sig_out.exists() and not overwrite:
        print(f"Signature already exists. Skipping: {sig_out}", flush=True)
        return "skipped_existing"

    params = sourmash_param_string(ksize, scaled, abundance)

    cmd = [
        "sourmash",
        "sketch",
        "dna",
        "-p", params,
        "--name", sample,
        "-o", str(sig_out),
        str(fastq),
    ]

    run(cmd, stderr_path=log_out)
    return "done"


def run_sourmash_gather(sig, database, gather_csv, log_out, ksize, overwrite):
    if gather_csv.exists() and not overwrite:
        print(f"Gather CSV already exists. Skipping: {gather_csv}", flush=True)
        return "skipped_existing"

    cmd = [
        "sourmash",
        "gather",
        str(sig),
        str(database),
        "-k", str(ksize),
        "-o", str(gather_csv),
    ]

    run(cmd, stderr_path=log_out)
    return "done"


def process_one_sample(fastq, database, ksize, scaled, abundance, overwrite):
    fastq = Path(fastq)
    sample = sample_name_from_fastq(fastq)
    sample_dir = fastq.parent

    print("========================================", flush=True)
    print(f"Processing sample: {sample}", flush=True)
    print("========================================", flush=True)

    sig_out = sample_dir / f"{sample}.sourmash.k{ksize}.sig.zip"
    sketch_log = sample_dir / f"{sample}.sourmash_sketch.k{ksize}.log"

    gather_csv = sample_dir / f"{sample}.sourmash_gather.k{ksize}.csv"
    gather_log = sample_dir / f"{sample}.sourmash_gather.k{ksize}.log"

    sketch_status = run_sourmash_sketch(
        fastq=fastq,
        sig_out=sig_out,
        log_out=sketch_log,
        sample=sample,
        ksize=ksize,
        scaled=scaled,
        abundance=abundance,
        overwrite=overwrite,
    )

    gather_status = run_sourmash_gather(
        sig=sig_out,
        database=database,
        gather_csv=gather_csv,
        log_out=gather_log,
        ksize=ksize,
        overwrite=overwrite,
    )

    return {
        "sample": sample,
        "input_fastq": str(fastq),
        "sourmash_signature": str(sig_out),
        "sourmash_gather_csv": str(gather_csv),
        "sourmash_sketch_log": str(sketch_log),
        "sourmash_gather_log": str(gather_log),
        "sketch_status": sketch_status,
        "gather_status": gather_status,
        "status": "done",
        "error": "",
    }


def combine_gather_csvs(gather_csvs, combined_csv):
    rows = []

    for csv_path in gather_csvs:
        csv_path = Path(csv_path)
        sample = csv_path.parent.name

        if not csv_path.exists():
            continue

        try:
            df = pd.read_csv(csv_path)
        except pd.errors.EmptyDataError:
            continue

        if df.empty:
            continue

        df.insert(0, "sample", sample)
        rows.append(df)

    if not rows:
        print("No non-empty sourmash gather CSVs to combine.", flush=True)
        return

    combined = pd.concat(rows, ignore_index=True)
    combined.to_csv(combined_csv, index=False)
    print(f"Wrote combined gather CSV: {combined_csv}", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run sourmash sketch/gather on per-sample decoy-filtered PAO1-unmapped FASTQ files."
    )
    parser.add_argument("--config", required=True, help="YAML config file")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    check_sourmash_available()

    input_report_dir = Path(config["input_report_dir"]).expanduser()
    fastq_pattern = config.get("fastq_pattern", "*.decoy_filtered.PAO1_unmapped.fastq.gz")

    sourmash_cfg = config["sourmash"]
    database = Path(sourmash_cfg["database"]).expanduser()
    ksize = sourmash_cfg.get("ksize", 21)
    scaled = sourmash_cfg.get("scaled", 1000)
    abundance = sourmash_cfg.get("abundance", True)

    settings = config.get("settings", {})
    overwrite = settings.get("overwrite_existing", False)
    combine_outputs = settings.get("combine_gather_csvs", True)
    combined_output_csv = input_report_dir / settings.get(
        "combined_output_csv",
        "sourmash_gather_combined.csv",
    )
    max_parallel_samples = int(settings.get("max_parallel_samples", 1))
    continue_on_error = settings.get("continue_on_error", False)

    if max_parallel_samples < 1:
        raise ValueError("settings.max_parallel_samples must be >= 1")

    if not database.exists():
        raise FileNotFoundError(f"Sourmash database not found: {database}")

    fastqs = find_fastqs(input_report_dir, fastq_pattern)

    print(f"Found {len(fastqs)} FASTQ files.", flush=True)
    print(f"Using sourmash database: {database}", flush=True)
    print(f"Using k={ksize}, scaled={scaled}, abundance={abundance}", flush=True)
    print(f"Running up to {max_parallel_samples} samples in parallel.", flush=True)

    manifest_rows = []

    if max_parallel_samples == 1:
        # Serial mode. Easier to debug.
        for fastq in fastqs:
            try:
                manifest_rows.append(
                    process_one_sample(fastq, database, ksize, scaled, abundance, overwrite)
                )
            except Exception as e:
                if not continue_on_error:
                    raise
                sample = sample_name_from_fastq(fastq)
                manifest_rows.append({
                    "sample": sample,
                    "input_fastq": str(fastq),
                    "sourmash_signature": "",
                    "sourmash_gather_csv": "",
                    "sourmash_sketch_log": "",
                    "sourmash_gather_log": "",
                    "sketch_status": "failed",
                    "gather_status": "failed",
                    "status": "failed",
                    "error": repr(e),
                })
    else:
        with ThreadPoolExecutor(max_workers=max_parallel_samples) as executor:
            futures = {
                executor.submit(
                    process_one_sample,
                    fastq,
                    database,
                    ksize,
                    scaled,
                    abundance,
                    overwrite,
                ): fastq
                for fastq in fastqs
            }

            for future in as_completed(futures):
                fastq = futures[future]
                sample = sample_name_from_fastq(fastq)

                try:
                    manifest_rows.append(future.result())
                except Exception as e:
                    if not continue_on_error:
                        raise

                    print(f"ERROR processing {sample}: {e}", file=sys.stderr, flush=True)
                    manifest_rows.append({
                        "sample": sample,
                        "input_fastq": str(fastq),
                        "sourmash_signature": "",
                        "sourmash_gather_csv": "",
                        "sourmash_sketch_log": "",
                        "sourmash_gather_log": "",
                        "sketch_status": "failed",
                        "gather_status": "failed",
                        "status": "failed",
                        "error": repr(e),
                    })

    manifest = pd.DataFrame(manifest_rows).sort_values("sample")
    manifest_csv = input_report_dir / "sourmash_manifest.csv"
    manifest.to_csv(manifest_csv, index=False)
    print(f"Wrote manifest: {manifest_csv}", flush=True)

    if combine_outputs:
        gather_csvs = [row["sourmash_gather_csv"] for row in manifest_rows if row.get("sourmash_gather_csv")]
        combine_gather_csvs(gather_csvs, combined_output_csv)

    failed = manifest[manifest["status"] == "failed"] if "status" in manifest.columns else pd.DataFrame()
    if not failed.empty:
        print(f"WARNING: {len(failed)} sample(s) failed. See sourmash_manifest.csv.", file=sys.stderr, flush=True)
        if not continue_on_error:
            sys.exit(1)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
