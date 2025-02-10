#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_mzmls = "/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/generate_library/input"
params.input_csvs = "/Users/shipei/Documents/projects/reverse_metabolomics/reverse_metabolomics/generate_library/input"
params.data_collector = "Minions Hello"
params.pi_name = "Pieter Dorrestein"
params.mass_detect_intensity_tolerance = 5e4
params.min_feature_height = 1.5e5
params.mz_tol_ppm = 5
params.ms2_explanation_cutoff = 0.60
params.adduct_type_mode = "simple"
params.core_adduct_filter = "simple"
params.component_precursor_check = "1"

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
    def escapedDataCollector = params.data_collector.replaceAll(/"/, '\\\\"')
    def escapedPiName = params.pi_name.replaceAll(/"/, '\\\\"')
    """
    python $TOOL_FOLDER/main_batch.py \
    --mzml_files ${mzmls} \
    --csv_files ${csvs} \
    --data_collector "${escapedDataCollector}" \
    --pi_name "${escapedPiName}" \
    --mass_detect_int_tol ${params.mass_detect_intensity_tolerance} \
    --min_feature_height ${params.min_feature_height} \
    --mz_tol_ppm ${params.mz_tol_ppm} \
    --ms2_explanation_cutoff ${params.ms2_explanation_cutoff} \
    --core_adduct_filter "${params.core_adduct_filter}" \
    --adduct_type_mode ${params.adduct_type_mode} \
    --plot \
    --component_precursor_check ${params.component_precursor_check} \
    --preprocessed_pkl "$TOOL_FOLDER/cmpd_name_to_mass.pkl"
    """
}

workflow {
    mzmls_ch = Channel.fromPath(params.input_mzmls + "/*.mzML")
    csvs_ch = Channel.fromPath(params.input_csvs + "/*.csv")

    createLibrary(mzmls_ch.collect(), csvs_ch.collect())
}