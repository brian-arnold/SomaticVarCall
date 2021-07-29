#

from functions import *
import random
configfile: "config.yaml"

INTERVALS = [str(i) for i in range(1,4)]

rule all:
    input:
        # custom func instead of expand with wildcards 
        # bc patients may have variable/unequal samples
        #make_all_input_varscan()
        make_all_input_mutect()

include: "rules/mutect2.smk"
include: "rules/varscan.smk"

