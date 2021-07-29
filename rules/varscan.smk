#

#rule mpileup_tumors:
#    input:
#        tumor = tumor_bam,
#        ref = config['reference']
#    output: 
#        mpile = "varscan/{patient}/{sample}.mpileup" 
#    resources: 
#        mem_mb = lambda wildcards, attempt: attempt * 3000
#    conda:
#        "../envs/varscan.yml"
#    shell:
#        "samtools mpileup -f {input.ref} -q 1 {input.tumor} > {output}"

#rule mpileup_normals:
#    input:
#        normal = normal_bam,
#        ref = config['reference']
#    output: 
#        mpile = "varscan/{patient}/normal.mpileup" 
#    resources: 
#        mem_mb = lambda wildcards, attempt: attempt * 3000
#    conda:
#        "../envs/varscan.yml"
#    shell:
#        "samtools mpileup -f {input.ref} -q 1 {input.normal} > {output}"

#rule varscan:
#    input:
#        tumor_mpile = "varscan/{patient}/{sample}.mpileup",
#        normal_mpile = "varscan/{patient}/normal.mpileup"
#    output: 
#        vcf = "varscan/{patient}/{sample}.snp"
#    params:
#        name = "varscan/{patient}/{sample}"
#    resources: 
#        mem_mb = lambda wildcards, attempt: attempt * 3000
#    conda:
#        "../envs/varscan.yml"
#    shell:
#        "varscan somatic {input.normal_mpile} {input.tumor_mpile} {params.name}"

rule varscan2:
    input:
        normal = normal_bam,
        tumor = tumor_bam,
        ref = config['reference']
    output: 
        vcf = "varscan/{patient}/{sample}.snp"
    params:
        name = "varscan/{patient}/{sample}",
        n_fifo = "varscan/{patient}/{sample}_normal.fifo",
        t_fifo = "varscan/{patient}/{sample}.fifo"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/varscan.yml"
    shell:
        "mkfifo {params.n_fifo}\n"
        "mkfifo {params.t_fifo}\n"
        "samtools mpileup -B -f {input.ref} -q 15 {input.normal} > {params.n_fifo} & "
        "samtools mpileup -B -f {input.ref} -q 15 {input.tumor} > {params.t_fifo} & "
        "varscan somatic {params.n_fifo} {params.t_fifo} {params.name}\n"

rule bcftools_mpileup_tumors:
    # generates an mpileup using bcftools, used in Decifer input
    input:
        tumor = tumor_bam,
        snvs = "varscan/{patient}/{sample}.snp",
        ref = config['reference']
    output: 
        mpile = "varscan/{patient}/{sample}.mpileup.tsv" 
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/hatchet.yml"
    shell:
        "bcftools mpileup {input.tumor} -f {input.ref} -T <(cut -f1-2 {input.snvs} | grep -v position) -a INFO/AD -Ou | "
        "bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%AD\n' > {output.mpile}"
