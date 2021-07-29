# WORKFLOW HAS 2 INDEPENDENT PATHS THAT CONVERGE
# mutect2 -> learnROM + gatherVcfs
# getPileupSummaries -> calculateContamination
# learnROM + gatherVcfs + calculateContamination -> FilterMutectCalls

rule mutect2:
    # scattered by patient (all samples from patient are together) and interval
    input:
        normal = normal_bam,
        tumors = tumor_bams, # plural bc collects all BAMs from patient
        ref = config['reference'],
        pon = config['mutect2']['pon'], # use panel of normals to filter out sequence artifacts
        pon_idx = config['mutect2']['pon'] + ".tbi", 
        gr = config['mutect2']['germline_resource'], # use germline resource as prior prob that normal sample carries an allele
        gr_idx = config['mutect2']['germline_resource'] + ".tbi" # use germline resource as prior prob that normal sample carries an allele
    output: 
        vcf = temp("mutect/{patient}/unfiltered_{interval}.vcf"),
        vcf_stats = temp("mutect/{patient}/unfiltered_{interval}.vcf.stats"),
        f1r2 = temp("mutect/{patient}/f1r2_{interval}.tar.gz")
    params:
        input_string = get_mutect_input,
        interval = "{interval}" # to use wildcards in 'run' statement below, specify them here
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk Mutect2 -R {input.ref} "
        "-O {output.vcf} "
        "-L {params.interval} "
        "--f1r2-tar-gz {output.f1r2} "
        "-pon {input.pon} "
        "-germline-resource {input.gr} "
        "{params.input_string}"

rule learnROM:
    # collect f1r2 files across scattered VCF jobs and learn the read orientation model
    input:
        # NOTE: the double curly brackets around 'patient' prevent the expand function from operating on that variable
        # thus, we expand by interval but not by patient, such that we gather by sample for each list value
        vcf = expand("mutect/{{patient}}/unfiltered_{interval}.vcf", interval=INTERVALS),
        f1r2 = expand("mutect/{{patient}}/f1r2_{interval}.tar.gz", interval=INTERVALS)
    output: 
        rom = "mutect/{patient}/read-orientation-model.tar.gz"
    params:
        input_string = get_ROM_input("{patient}", INTERVALS)
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk LearnReadOrientationModel -O {output.rom} "
        "{params.input_string} "

rule getPileupSummaries:
    # get pileup info at selected positions for EACH tumor sample, scattered by tumor
    input:
        vcf = tumor_bam, # singular bc collects a single BAM from a patient
        common_variants = config['mutect2']['common_alleles'],
        idx = config['mutect2']['common_alleles'] + ".tbi"
    output: 
        temp("mutect/{patient}/{sample}_getpileupsummaries.table")
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk GetPileupSummaries "
        "-I {input.vcf} "
        "-V {input.common_variants} "
        "-L {input.common_variants} "
        "-O {output} "

rule calculateContamination:
    # get pileup info at selected positions for EACH tumor sample, scattered by tumor
    input:
        "mutect/{patient}/{sample}_getpileupsummaries.table"
    output: 
        seg = temp("mutect/{patient}/{sample}_segments.table"),
        contam = temp("mutect/{patient}/{sample}_contamination.table")
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk CalculateContamination "
        "-I {input} "
        "-tumor-segmentation {output.seg} "
        "-O {output.contam} "

rule gatherVCFs:
    input:
        vcf = expand("mutect/{{patient}}/unfiltered_{interval}.vcf", interval=INTERVALS)
    output:
        vcf = temp("mutect/{patient}/unfiltered.vcf")
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    params:
        input_vcf_string = get_gatherVCFs_input("{patient}", INTERVALS), # precede each VCF with "-I"
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk GatherVcfs "
        "{params.input_vcf_string} "
        "-O {output.vcf}\n"

rule gatherStats:
    input:
        vcf_stats = expand("mutect/{{patient}}/unfiltered_{interval}.vcf.stats", interval=INTERVALS)
    output:
        stats = temp("mutect/{patient}/merged.stats")
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    params:
        input_stat_string = get_stats_input("{patient}", INTERVALS) # precede each VCF with "-stats"
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk MergeMutectStats "
        "{params.input_stat_string} "
        "-O {output.stats}"

rule filterMutectCalls:
    # Applies filters to unfiltered VCF, labeling somatic mutations as PASS, 
    # and also labeling germline mutations in the FILTER column
    input:
        # collect segments and contamination tables across all samples for a patient
        seg_tables = get_seg_tables,
        contam_tables = get_contam_tables,
        rom = "mutect/{patient}/read-orientation-model.tar.gz",
        vcf = "mutect/{patient}/unfiltered.vcf",
        ref = config['reference'],
        stats = "mutect/{patient}/merged.stats"
    output: 
        vcf = "mutect/{patient}/filtered.vcf"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    params:
        seg_string = get_seg_input("{patient}"), # precede each table with "--tumor-segmentation"
        contam_string = get_contam_input("{patient}") # precede each table with "--contamination-table"
    conda:
        "../envs/varcall.yml"
    shell:
        "gatk FilterMutectCalls -V {input.vcf} "
        "-R {input.ref} "
        "--stats {input.stats} "
        "{params.seg_string} "
        "{params.contam_string} "
        "--ob-priors {input.rom} "
        "-O {output.vcf}"

rule filterVcfPass:
    # Filter the mutect VCF for PASS, otherwise it still contains germline variants
    # Only take Bi-allelic sites
    input:
        vcf = "mutect/{patient}/filtered.vcf"
    output: 
        vcf = "mutect/{patient}/filtered_PASS_biallelic.vcf"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varcall.yml"
    shell:
        "vcftools --vcf {input.vcf} --remove-filtered-all --recode --stdout "
        "--min-alleles 2 --max-alleles 2 "
        "> {output.vcf}"




