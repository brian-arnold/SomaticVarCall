Simple snakemake workflow to detect somatic variants in tumors. At the moment, this code is optimized to run on Computer Science Department's HPC at Princeton University.

Please make sure that BAM files have been pre-processed for variant calling. This may include e.g. marking duplicates for WGS BAMs.

Not included in this repo are various external resources that facilitate identifying somatic mutations and excluding germline mutations. These can be found in the `resources` subdir.

Many resources for Mutect2 can be found in google cloud storage. See what is available using:

gsutil ls gs://gatk-best-practices/

and copy relevant files into the appropriate directories. Due to their sizes, these files are excluded here but may easily be re-downloaded using the `gsutil cp` command.

A panel of normals is used for filtering out technical artifacts. E.g. in dir `resources/panel_of_normals/hg19`:
- Mutect2-WGS-panel-b37.vcf.gz
- Mutect2-exome-panel.vcf.gz
- Mutect2-WGS-panel-b37.vcf.gz.tbi
- Mutect2-exome-panel.vcf.gz.tbi

A germline resource is used such that allele frequencies of germline mutations in the population are used as a prior for calling variants in the normal sample. E.g. in dir `resources/germline_resource/hg19`:
- af-only-gnomad.raw.sites.vcf.gz
- af-only-gnomad.raw.sites.vcf.gz.tbi

A list of common alleles is used to get pileup information at these sites and is used to calculate contamination. E.g. in dir `resources/common_alleles/hg19`:
- small_exac_common_3.vcf.gz
- small_exac_common_3.vcf.gz.tbi

