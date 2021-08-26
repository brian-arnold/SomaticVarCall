import os
from collections import defaultdict
from snakemake import load_configfile

#config = load_configfile('config_files/config_stjude.yaml')
config = load_configfile('config_files/config_Gundem.yaml')

# This file contains (1) functions used by the snakemake rules
# and (2) global variables accessed by these functions, e.g.
# the dictionaries that organize samples by patient

##################
# FUNCTIONS
##################

def get_samples_by_patient():
    # Builds two data structures for grouping samples by patient: 
    # tumors[patient][sample] = bam 
    # normals[patient] = bam 
    tumors = defaultdict(dict)
    normals = defaultdict()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patient = info[0]
                bam = info[2]
                if info[1] == "tumor":
                    sample = get_sample_name( bam )
                    tumors[patient][sample] = bam
                elif info[1] == "normal":
                    normals[patient] = bam
                else:
                    sys.exit("incorrect samples.txt file")
        return (tumors, normals)

def make_all_input_varscan():
    # output is 1 vcf per tumor sample per patient, 
    # with samples grouped in directories according to patient names
    with open(config["samples"], "r") as f:
        outfiles = []
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                if info[1] == "tumor":
                    sample = get_sample_name( info[2] )
                    path = os.path.join( "varscan", info[0], f"{sample}.mpileup.tsv")
                    #path = os.path.join( "varscan", info[0], f"{sample}.snp")
                    outfiles.append( path ) 
        return outfiles

def make_all_input_mutect():
    # output is 1 vcf per patient, 
    patients = set()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patients.add( info[0] )
        #paths =  [  os.path.join( "mutect", p, f"unfiltered_{i}.vcf") for p in patients for i in range(1,24) ]
        #paths =  [  os.path.join( "mutect", p, "read-orientation-model.tar.gz") for p in patients ]
        #paths =  [  os.path.join( "mutect", p, "filtered.vcf") for p in patients ]
        paths =  [  os.path.join( "mutect", p, "filtered_PASS_biallelic.vcf.gz") for p in patients ]
        print(paths)
        return paths

def make_all_input_intersect():
    # output is 1 vcf per patient, 
    patients = set()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patients.add( info[0] )
        #paths =  [  os.path.join( "mutect", p, f"unfiltered_{i}.vcf") for p in patients for i in range(1,24) ]
        #paths =  [  os.path.join( "mutect", p, "read-orientation-model.tar.gz") for p in patients ]
        #paths =  [  os.path.join( "mutect", p, "filtered.vcf") for p in patients ]
        paths =  [  os.path.join(f"intersect/{p}/intersection.vcf.gz") for p in patients ]
        print(paths)
        return paths

def make_all_input_strelka():
    # output is 1 vcf per patient, 
    patients = set()
    with open(config["samples"], "r") as f:
        for l in f:
            if not l.startswith("PATIENT_ID"):
                info = l.strip().split() 
                patients.add( info[0] )
        #paths =  [  os.path.join( "strelka", p, "combined_1.vcf.gz") for p in patients ]
        paths =  [  os.path.join( "strelka", p, "multisample_calls.vcf.gz") for p in patients ]
        print(paths)
        return paths

def get_sample_name(x):
    sample = os.path.basename(x)
    sample = sample.replace(".bam", "") 
    return sample
        
def tumor_bam(wildcards):
    # tumors dict defined below
    tumor_in = tumors[wildcards.patient][wildcards.sample]
    return tumor_in

def tumor_bams(wildcards):
    # tumors dict defined below
    tumors_in = [ tumors[wildcards.patient][s] for s in tumors[wildcards.patient] ]
    return tumors_in

def get_seg_tables(wildcards):
    x = [ f"mutect/{wildcards.patient}/{s}_segments.table" for s in sorted(tumors[wildcards.patient]) ]
    return x

def get_contam_tables(wildcards):
    x = [ f"mutect/{wildcards.patient}/{s}_contamination.table" for s in sorted(tumors[wildcards.patient]) ]
    return x

def normal_bam(wildcards):
    # normals dict defined below
    normal_in = normals[wildcards.patient]
    return normal_in

##########
# functions for prefixing files with argument names
##########

# Mutect

def get_mutect_input(wildcards):
    tumors_in = []
    [ tumors_in.extend( ["-I", tumors[wildcards.patient][s]] ) for s in sorted(tumors[wildcards.patient]) ]
    mutect_input = " ".join(tumors_in)
    return mutect_input 

def get_ROM_input(patient, INTERVALS):
    x = []
    [ x.extend(["-I", f"mutect/{patient}/f1r2_{i}.tar.gz"]) for i in INTERVALS ]
    r = " ".join(x)
    return r

def get_gatherVCFs_input(patient, INTERVALS):
    x = []
    [ x.extend(["-I", f"mutect/{patient}/unfiltered_{i}.vcf"]) for i in INTERVALS ]
    r = " ".join(x)
    return r

def get_stats_input(patient, INTERVALS):
    x = []
    [ x.extend(["-stats", f"mutect/{patient}/unfiltered_{i}.vcf.stats"]) for i in INTERVALS ]
    r = " ".join(x)
    return r

def get_contam_input(patient):
    x = []
    [ x.extend([ "--contamination-table", f"mutect/{patient}/{s}_contamination.table"]) for s in sorted(tumors[patient]) ]
    return x

def get_seg_input(patient):
    x = []
    [ x.extend([ "--tumor-segmentation", f"mutect/{patient}/{s}_segments.table"]) for s in sorted(tumors[patient]) ]
    return x

# Strelka

def get_vcfs_per_patient(wildcards):
   x = [f"strelka/{wildcards.patient}/{s}/results/variants/somatic.snvs.tumor.renamed.vcf.gz" for s in sorted(tumors[wildcards.patient])]
   y = [f"strelka/{wildcards.patient}/{s}/results/variants/somatic.snvs.tumor.renamed.vcf.gz.tbi" for s in sorted(tumors[wildcards.patient])]
   return x+y

def get_vcfs_per_patient_string(wildcards):
    x = [f"strelka/{wildcards.patient}/{s}/results/variants/somatic.snvs.tumor.renamed.vcf.gz" for s in sorted(tumors[wildcards.patient])]
    r = " ".join(x)
    return r

def get_vcfs_per_patient2(wildcards):
   x = [f"strelka/{wildcards.patient}/{s}/joint/results/variants/somatic.snvs.tumor.renamed.vcf.gz" for s in sorted(tumors[wildcards.patient])]
   y = [f"strelka/{wildcards.patient}/{s}/joint/results/variants/somatic.snvs.tumor.renamed.vcf.gz.tbi" for s in sorted(tumors[wildcards.patient])]
   return x+y

def get_vcfs_per_patient_string2(wildcards):
    x = [f"strelka/{wildcards.patient}/{s}/joint/results/variants/somatic.snvs.tumor.renamed.vcf.gz" for s in sorted(tumors[wildcards.patient])]
    r = " ".join(x)
    return r

##################
# GLOBALS
##################

# Data structures for grouping samples by patient:
# tumors[patient][sample] = bam
# normals[patient] = bam
tumors, normals = get_samples_by_patient() 

# intervals for mutect2
INTERVALS = [str(i) for i in range(1,config['autosomes_to_genotype']+1)]
if config['genotype_sexchroms']:
    INTERVALS.extend(["X","Y"])



