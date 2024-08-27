#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_mzmls = "/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/generate_library/input/*.mzML"
params.input_csvs = "/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/generate_library/input/*.csv"
params.data_collector = "Minions"

TOOL_FOLDER = "$baseDir/bin"

process createLibrary {
    publishDir "./nf_output", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file input_mzmls
    file input_csvs

    output:
    file 'all_library.tsv'
    file '*.tsv' optional true
    file '*.pdf' optional true
    file '*.svg' optional true

    """
    echo "${input_mzmls.join('\n')}" > temp_mzml_list.txt
    echo "${input_csvs.join('\n')}" > temp_csv_list.txt
    python $TOOL_FOLDER/main_batch.py --mzml_list temp_mzml_list.txt --csv_list temp_csv_list.txt \
    --data_collector ${params.data_collector}
    """
}


workflow {
    mzmls_ch = Channel.fromPath(params.input_mzmls)
    csvs_ch = Channel.fromPath(params.input_csvs)

    createLibrary(mzmls_ch, csvs_ch)
}
