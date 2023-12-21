include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "tables/calls_update.tsv.gz"


##### Modules #####

include: "rules/subset.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
