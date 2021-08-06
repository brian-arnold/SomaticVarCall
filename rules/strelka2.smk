# Since strelka2 does not support multi sample calling, the strategy here is to perform
# calling on individual tumor samples, combine the calls, then recall variants using the 
# combined VCF from individual samples.

rule strelka2_config:
    # create strelka configuration script
    input:
        normal = normal_bam,
        tumor = tumor_bam, # plural bc collects all BAMs from patient
        ref = config['reference']
    output: 
        vcf = "strelka/{patient}/{sample}/runWorkflow.py"
    params:
        rundir = "strelka/{patient}/{sample}/"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input.normal} "
        "--tumorBam {input.tumor} "
        "--referenceFasta {input.ref} "
        "--runDir {params.rundir} "

rule strelka2_exec:
    # run each strelka configureation script
    input:
        runpy = "strelka/{patient}/{sample}/runWorkflow.py"
    output: 
        vcf = "strelka/{patient}/{sample}/results/variants/somatic.snvs.vcf.gz"
    params:
        cores = config['strelka2']['cores']
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "{input.runpy} -m local -j {params.cores}"

rule subset_1:
    # Rename VCF files so that they can be combined. By default, all VCFs have NORMAL and TUMOR as sample names
    input:
        vcf = "strelka/{patient}/{sample}/results/variants/somatic.snvs.vcf.gz" 
    output: 
        vcf = "strelka/{patient}/{sample}/results/variants/somatic.snvs.tumor.vcf.gz"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "bcftools view --samples TUMOR "
        "--output-type z "
        "--output {output.vcf} "
        "{input.vcf}"

rule rename_1:
    # Rename VCF files so that they can be combined. By default, all VCFs have NORMAL and TUMOR as sample names
    input:
        vcf = "strelka/{patient}/{sample}/results/variants/somatic.snvs.tumor.vcf.gz"
    output: 
        vcf = "strelka/{patient}/{sample}/results/variants/somatic.snvs.tumor.renamed.vcf.gz",
        tbi = "strelka/{patient}/{sample}/results/variants/somatic.snvs.tumor.renamed.vcf.gz.tbi"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    params:
        sample = "{sample}"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "picard RenameSampleInVcf "
        "INPUT={input.vcf} "
        "OUTPUT={output.vcf} "
        "NEW_SAMPLE_NAME={params.sample}\n"
        "tabix {output.vcf}"

rule combine_calls1:
    # combine single sample calls across samples for each patient, no need to gather by patient at this point
    # just gather by sample
    input:
        vcfs = get_vcfs_per_patient
    output: 
        vcf = "strelka/{patient}/combined_1.vcf.gz",
        idx = "strelka/{patient}/combined_1.vcf.gz.tbi"
    params:
        vcfs = get_vcfs_per_patient_string
    wildcard_constraints:
        patient = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "bcftools merge --output-type z --output {output.vcf} {params.vcfs}\n"
        "tabix {output.vcf}"


rule strelka2_config2:
    # create strelka configuration script
    input:
        normal = normal_bam,
        tumor = tumor_bam, # plural bc collects all BAMs from patient
        ref = config['reference'],
        vcf = "strelka/{patient}/combined_1.vcf.gz",
        idx = "strelka/{patient}/combined_1.vcf.gz.tbi"
    output: 
        vcf = "strelka/{patient}/{sample}/joint/runWorkflow.py"
    params:
        rundir = "strelka/{patient}/{sample}/joint/"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "configureStrelkaSomaticWorkflow.py "
        "--normalBam {input.normal} "
        "--tumorBam {input.tumor} "
        "--referenceFasta {input.ref} "
        "--forcedGT {input.vcf} "
        "--runDir {params.rundir} "

rule strelka2_exec2:
    # run each strelka configureation script
    input:
        runpy = "strelka/{patient}/{sample}/joint/runWorkflow.py"
    output: 
        vcf = "strelka/{patient}/{sample}/joint/results/variants/somatic.snvs.vcf.gz"
    params:
        cores = 8
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "{input.runpy} -m local -j {params.cores}"

rule subset_2:
    # Rename VCF files so that they can be combined. By default, all VCFs have NORMAL and TUMOR as sample names
    input:
        vcf = "strelka/{patient}/{sample}/joint/results/variants/somatic.snvs.vcf.gz" 
    output: 
        vcf = "strelka/{patient}/{sample}/joint/results/variants/somatic.snvs.tumor.vcf.gz"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "bcftools view --samples TUMOR "
        "--output-type z "
        "--output {output.vcf} "
        "{input.vcf}"

rule rename_2:
    # Rename VCF files so that they can be combined. By default, all VCFs have NORMAL and TUMOR as sample names
    input:
        vcf = "strelka/{patient}/{sample}/joint/results/variants/somatic.snvs.tumor.vcf.gz"
    output: 
        vcf = "strelka/{patient}/{sample}/joint/results/variants/somatic.snvs.tumor.renamed.vcf.gz",
        tbi = "strelka/{patient}/{sample}/joint/results/variants/somatic.snvs.tumor.renamed.vcf.gz.tbi"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    params:
        sample = "{sample}"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "picard RenameSampleInVcf "
        "INPUT={input.vcf} "
        "OUTPUT={output.vcf} "
        "NEW_SAMPLE_NAME={params.sample}\n"
        "tabix {output.vcf}"

rule combine_calls2:
    # combine single sample calls across samples for each patient, no need to gather by patient at this point
    # just gather by sample
    input:
        vcfs = get_vcfs_per_patient2
    output: 
        vcf = "strelka/{patient}/multisample_calls.vcf.gz",
        idx = "strelka/{patient}/multisample_calls.vcf.gz.tbi"
    params:
        vcfs = get_vcfs_per_patient_string2
    wildcard_constraints:
        patient = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "bcftools merge --output-type z --output {output.vcf} {params.vcfs}\n"
        "tabix {output.vcf}"

