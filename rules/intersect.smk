rule intesect_vcfs:
    # create strelka configuration script
    input:
        mutect_vcf = "mutect/{patient}/filtered_PASS_biallelic.vcf.gz",
        strelka_vcf = "strelka/{patient}/multisample_calls.vcf.gz"
    output: 
        vcf = "intersect/{patient}/intersection.vcf.gz"
    params:
        rundir = "intersect/{patient}/"
    wildcard_constraints:
        patient = "\w+",
        sample = "\w+"
    resources: 
        mem_mb = lambda wildcards, attempt: attempt * 3000
    conda:
        "../envs/strelka2.yml"
    shell:
        "bcftools isec -n~11 "
        "--output-type z "
        "-p {params.rundir} "
        "{input.mutect_vcf} {input.strelka_vcf}\n"
        "mv {params.rundir}0000.vcf.gz {output.vcf}\n"

