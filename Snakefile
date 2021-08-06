#

from functions import *
import random

rule all:
    input:
        # custom func instead of expand with wildcards 
        # bc patients may have variable/unequal samples
        #make_all_input_strelka()
        #make_all_input_mutect()
        make_all_input_intersect()

include: "rules/mutect2.smk"
include: "rules/strelka2.smk"
include: "rules/intersect.smk"

