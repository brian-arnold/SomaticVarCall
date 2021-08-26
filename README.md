Simple snakemake workflow to detect somatic variants in tumors. At the moment, this code is optimized to run on Computer Science Department's HPC at Princeton University. To run, create  

1. A samples file (e.g. like ones in `samples_files` directory) that gives the location of BAM files, whether it's tumor/normal, and what patient it belongs to
2. A config file (e.g. like ones in `config_files` directory), and specify this file within the `functions.py` module. This file also specifies the location of supplementary files needed for Mutect2.
3. A conda environment with snakemake installed.

After these are created, you may run the workflow locally (e.g. on turing) using the `run_local.sh` script or on the CS cluster via the SLURM scheduler using the `run_cluster.sh` script. You may have to modify these files based on your setup, for instance the name of your conda environment, and I also include a `source` command in there to enable conda for external jobs. To test the workflow before executing these scripts, type `snakemake -n -p` on the command line after activating your conda environment containing snakemake.


Please make sure that BAM files have been pre-processed for variant calling. This may include e.g. marking duplicates for WGS BAMs.

Not included in this repo are various external resources that facilitate identifying somatic mutations and excluding germline mutations. These can be found in the `resources` subdir.

Many resources for Mutect2 can be found in google cloud storage. See what is available using:

gsutil ls gs://gatk-best-practices/

and copy relevant files into the appropriate directories. Due to their sizes, these files are excluded here but may easily be re-downloaded using the `gsutil cp` command.

A panel of normals is used for filtering out technical artifacts. E.g. in dir `resources/panel_of_normals/hg19`:
- Mutect2-WGS-panel-b37.vcf.gz
- Mutect2-WGS-panel-b37.vcf.gz.tbi
- Mutect2-exome-panel.vcf.gz
- Mutect2-exome-panel.vcf.gz.tbi

A germline resource is used such that allele frequencies of germline mutations in the population are used as a prior for calling variants in the normal sample. E.g. in dir `resources/germline_resource/hg19`:
- af-only-gnomad.raw.sites.vcf.gz
- af-only-gnomad.raw.sites.vcf.gz.tbi

A list of common alleles is used to get pileup information at these sites and is used to calculate contamination. E.g. in dir `resources/common_alleles/hg19`:
- small_exac_common_3.vcf.gz
- small_exac_common_3.vcf.gz.tbi


Strelka needs a reference genome that has the exact same contigs as those seen in the BAM files (including small contigs). Extra reference genomes may be found in the resources directory, and the reference `hs37d5.fa` was downloaded from GCP at gs://genomics-public-data/references.

