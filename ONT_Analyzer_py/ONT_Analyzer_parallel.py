import subprocess
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
import Run_Quality_Analyzer

class FastqPipeline:
    """Class to process Oxford Nanopore (ONT) sequencing data efficiently."""

    def __init__(self, base_path, reference_genome, summary_file, num_workers=4):
        """Initialize paths and settings."""
        self.base_path = base_path
        self.results_path = f"{base_path}/results"
        self.input_path = f"{base_path}/data_separated/"
        self.reference_genome = reference_genome
        self.summary_file = summary_file
        self.conda_env = "long_reads"
        self.num_workers = num_workers if num_workers else os.cpu_count()
        self.barcodes = range(1, 7)  # Adjust barcode range if needed

        # Ensure necessary directories exist
        os.makedirs(self.results_path, exist_ok=True)

    def run_command(self, command):
        """Executes a shell command and handles errors."""
        try:
            return subprocess.run(command, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing: {' '.join(command)}\n{e}")
            exit(1)

    def run_quality_check(self):
        """Run PycoQC for quality assessment (single run)."""
        print("Running PycoQC...")
        Run_Quality_Analyzer.Graphs(self.summary_file)
        cmd = ["conda", "run", "-n", self.conda_env, "pycoQC", "--summary_file", self.summary_file,
               "--html_outfile", f"/Users/giuse/pythonProject/Mycodes/MyProjects/Long-Reads_Analyzer/results/Run_quality.html"]
        return self.run_command(cmd)

    def process_barcode(self, barcode, tool):
        """Process a single barcode using a specified tool."""
        match tool:
            case "filtlong":
                input_file = f"{self.base_path}/data_separated/barcode_{barcode}.fastq"
                output_file = f"{self.results_path}/barcode_{barcode}.filtered.fastq"
                cmd = ["conda", "run", "filtlong", "--min_length", "1000",
                       "--keep_percent", "90", "--trim", "-a", self.reference_genome, input_file]

            case "fastqc":
                input_file = f"{self.results_path}/barcode_{barcode}.filtered.fastq"
                cmd = ["conda", "run", "-n", self.conda_env, "fastqc", input_file, "--outdir", f"{self.results_path}/fastqc"]

            #case "genome_assembly":
            #    cmd = ["conda", "run", "-n", self.conda_env, "flye", "--nano-hq",
            #           f"{self.results_path}/barcode_{barcode}.filtered.fastq", "--genome-size", "4.6m",
            #           "--out-dir", f"{self.results_path}/flye_assembly", "--threads", "7"]

            case "mapping":
                input_file = f"{self.results_path}/barcode_{barcode}.filtered.fastq"
                sam_output = f"{self.results_path}/sample_{barcode}.sam"
                cmd = ["minimap2", "-ax", "map-ont", "-t", "4", self.reference_genome, input_file, ">", sam_output]

            case "samtools":
                input_bam = f"{self.results_path}/sample_{barcode}.bam"
                cmd = ["samtools", "view", "-bS", "-@", "4", f"{self.results_path}/sample_{barcode}.sam",
                       "|", "samtools", "sort", "-o", input_bam]
                self.run_command(cmd)
                cmd = ["samtools", "index", input_bam]

            case "picard":
                input_bam = f"{self.results_path}/sample_{barcode}.bam"
                output_bam = f"{self.results_path}/sample_{barcode}.markdup.bam"
                cmd = ["conda", "run", "-n", self.conda_env, "picard", "-Xmx8G", "MarkDuplicates",
                       f"I={input_bam}", f"O={output_bam}", "REMOVE_DUPLICATES=true", "CREATE_INDEX=true"]

            case "variant_calling":
                input_bam = f"{self.results_path}/sample_{barcode}.markdup.bam"
                vcf_output = f"{self.results_path}/sample_{barcode}.vcf"
                cmd = ["conda", "run", "-n", self.conda_env, "sniffles", "--input", input_bam, "--vcf", vcf_output]

            case _:
                raise ValueError(f"Unsupported tool: {tool}")

        print(f"Running {tool} on barcode {barcode}...")
        return self.run_command(cmd)

    def run_parallel(self, tool):
        """Run barcode-dependent steps in parallel."""
        print(f"Starting parallel processing for: {tool}")
        with ProcessPoolExecutor(max_workers=self.num_workers) as executor:
            futures = {executor.submit(self.process_barcode, barcode, tool): barcode for barcode in self.barcodes}

            for future in as_completed(futures):
                barcode = futures[future]
                try:
                    future.result()  # Ensure no errors
                    print(f"{tool} completed for barcode {barcode}")
                except Exception as e:
                    print(f"{tool} failed for barcode {barcode}: {e}")

    def run_pipeline(self):
        """Run the full ONT sequencing pipeline."""
        print("Starting ONT sequencing pipeline...")

        # Step 1: Quality Check (single run)
        #self.run_quality_check()
        #print("Run Quality Check Completed!")

        # Step 2: Parallel Processing
        self.run_parallel("filtlong")
        self.run_parallel("fastqc")
        self.run_parallel("genome_assembly")
        self.run_parallel("mapping")
        self.run_parallel("samtools")
        self.run_parallel("picard")
        self.run_parallel("variant_calling")

        print("Pipeline completed successfully!")


# Run Pipeline
if __name__ == "__main__":
    """ Define Paths - Replace with actual file locations """
    BASE_PATH = "/Users/giuse/pythonProject/Mycodes/MyProjects"
    REFERENCE_GENOME_PATH = f"{BASE_PATH}/Escherichia_coli_reference"
    SUMMARY_FILE_PATH = f"{BASE_PATH}/ONT_simulated_summary.txt"

    # Check required files
    for path in [REFERENCE_GENOME_PATH, SUMMARY_FILE_PATH]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

    # Initialize and run the pipeline
    pipeline = FastqPipeline(BASE_PATH, REFERENCE_GENOME_PATH, SUMMARY_FILE_PATH)
    pipeline.run_pipeline()
