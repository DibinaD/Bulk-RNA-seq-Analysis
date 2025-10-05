import subprocess
import time
import os

# List of SRA IDs to download
sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540","SRR7179541"
]

# Make output folder for FASTQ
os.makedirs("fastq", exist_ok=True)

# Step 1: Prefetch each SRA file
for sra_id in sra_numbers:
    print(f"\n=== Downloading: {sra_id} ===")
    prefetch_cmd = f"prefetch {sra_id}"
    print("Command:", prefetch_cmd)

    start_time = time.time()
    subprocess.call(prefetch_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ Download time for {sra_id}: {elapsed_min:.2f} minutes")

# Step 2: Convert SRA to FASTQ
for sra_id in sra_numbers:
    sra_path = os.path.expanduser(f"~/ncbi/public/sra/{sra_id}.sra")
    print(f"\n=== Generating FASTQ for: {sra_id} ===")

    fasterq_cmd = (
        f"fasterq-dump {sra_path} --outdir fastq --split-files --threads 4"
    )
    print("Command:", fasterq_cmd)

    start_time = time.time()
    subprocess.call(fasterq_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ FASTQ generation time for {sra_id}: {elapsed_min:.2f} minutes")

print("\n🎉 All downloads and FASTQ conversions complete!")

