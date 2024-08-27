#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_mzmls = "/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/generate_library/input"
params.input_csvs = "/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/generate_library/input"
params.data_collector = "Minions"
params.ms2_explanation_cutoff = 0.60

TOOL_FOLDER = "$baseDir/bin"

process createLibrary {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path mzmls
    path csvs

    output:
    path 'all_library.tsv'
    path '*.tsv', optional: true
    path '*.pdf', optional: true
    path '*.png', optional: true

    script:
    """
    python $TOOL_FOLDER/main_batch.py \
    --mzml_files ${mzmls} \
    --csv_files ${csvs} \
    --data_collector ${params.data_collector} \
    --ms2_explanation_cutoff ${params.ms2_explanation_cutoff} --plot
    """
}


workflow {
    mzmls_ch = Channel.fromPath(params.input_mzmls + "/*.mzML")
    csvs_ch = Channel.fromPath(params.input_csvs + "/*.csv")

    createLibrary(mzmls_ch.collect(), csvs_ch.collect())
}
