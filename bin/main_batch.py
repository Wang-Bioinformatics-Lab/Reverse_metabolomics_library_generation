import argparse
import os
import shutil
import pandas as pd

from cmpd import prepare_cmpd_df
from create_library import create_library, append_file_summary, plot_ms2_annotation_distribution
from feature_extraction import feature_extraction_single, plot_all_ms2, plot_all_eic, plot_mz_rt


def main_batch(mzml_files, csv_files,
               data_collector='Minions',
               pi_name='Pieter Dorrestein',
               mass_detect_int_tol=5e4,
               min_feature_height=1.5e5,
               mz_tol_ppm=5,
               ms2_explanation_cutoff=0.60,
               core_adduct_filter='full',
               adduct_type_mode='full',
               plot=False,
               component_precursor_check=True,
               preprocessed_pkl_path=None):
    """
    Process a batch of mzML files and csv files.
    """

    os.makedirs('details', exist_ok=True)
    os.makedirs('tmp_data', exist_ok=True)

    # Load the csv file
    print('Calculating compound masses...')
    all_cmpd_df = pd.DataFrame()
    for csv in csv_files:
        print('Loading', csv)
        cmpd_df = prepare_cmpd_df(csv, adduct_type_mode)
        all_cmpd_df = pd.concat([all_cmpd_df, cmpd_df])

    # Get unique mzmls
    unique_mzmls = all_cmpd_df['unique_sample_id'].unique()

    # split all_cmpd_df into multiple dataframes by unique_sample_id and save them
    for mzml in unique_mzmls:
        cmpd_df = all_cmpd_df[all_cmpd_df['unique_sample_id'] == mzml]
        mzml_basename = mzml.split('.mz')[0]
        cmpd_df.to_csv(f'tmp_data/{mzml_basename}_cmpd_df.csv', index=False)

    del all_cmpd_df, unique_mzmls  # free up memory

    # scan number
    scans_no = 1

    # file summary rows
    all_file_summary_rows = []

    # Load the mzml file
    for mzml in mzml_files:
        print('Processing', mzml)

        mzml_name = os.path.basename(mzml)
        mzml_basename = mzml_name.split('.mz')[0]
        cmpd_df_path = f'tmp_data/{mzml_basename}_cmpd_df.csv'
        if not os.path.exists(cmpd_df_path):
            print(f'Compound list for {mzml_basename} not found. Skipping...')
            continue

        # Extract features from mzML file
        try:
            feature_df, ion_mode = feature_extraction_single(file_path=mzml,
                                                             mass_detect_int_tol=mass_detect_int_tol,
                                                             min_feature_height=min_feature_height,
                                                             save=False, out_dir='details')
        except:
            print(f'Error extracting features from {mzml}. Skipping...')
            continue

        # Load the compound list for the mzml file
        cmpd_df = pd.read_csv(cmpd_df_path, low_memory=False)
        # Filter compound list by ion mode
        cmpd_df = cmpd_df[cmpd_df['ion_mode'] == ion_mode].reset_index(drop=True)

        # Filter library
        print('Creating MS/MS library...')
        df, library_df, scans_no = create_library(cmpd_df, feature_df, ion_mode,
                                                  data_collector, pi_name, mzml_basename,
                                                  mz_tol_ppm=mz_tol_ppm,
                                                  filter_library=True,
                                                  ms2_explanation_cutoff=ms2_explanation_cutoff,
                                                  core_adduct_filter=core_adduct_filter,
                                                  component_precursor_check=component_precursor_check,
                                                  preprocessed_pkl_path=preprocessed_pkl_path,
                                                  scans_start=scans_no)

        # Append file summary rows
        all_file_summary_rows = append_file_summary(all_file_summary_rows, mzml_name, cmpd_df, feature_df, df)

        if library_df is not None:
            library_df.to_csv(f'tmp_data/{mzml_basename}_library.tsv', sep='\t', index=False, na_rep='N/A')

        # Save
        if df is not None:
            df.to_csv(f'details/{mzml_basename}_metadata.tsv', sep='\t', index=False)

        if plot and df is not None and feature_df is not None:
            try:
                # Plot all MS2 spectra
                print('Plotting all spectra...')
                plot_all_ms2(df, mzml, out_dir='details')

                # Plot all EICs
                print('Plotting all EICs...')
                plot_all_eic(df, mzml, out_dir='details')

                # Plot mz-rt scatter plot
                print('Plotting mz-rt scatter plot...')
                plot_mz_rt(feature_df, df, mzml_basename, out_dir='details')
            except Exception as e:
                print(e)

    # Merge all library files
    # find all library files (tmp_data/*_library.tsv)
    library_files = [f for f in os.listdir('tmp_data') if f.endswith('_library.tsv')]
    if library_files:
        all_library_df = pd.DataFrame()
        for library_file in library_files:
            library_df = pd.read_csv(f'tmp_data/{library_file}', sep='\t', low_memory=False)
            all_library_df = pd.concat([all_library_df, library_df], ignore_index=True)

        # Save the merged library
        all_library_df.to_csv('all_library.tsv', sep='\t', index=False, na_rep='N/A')

    # clean up the entire tmp_data folder
    if os.path.exists('tmp_data'):
        shutil.rmtree('tmp_data')

    # Save the file summary
    file_summary_df = pd.DataFrame(all_file_summary_rows)
    file_summary_df.to_csv('file_summary.tsv', sep='\t', index=False)

    # # Plot MS2 annotation distribution
    # plot_ms2_annotation_distribution(file_summary_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a batch of mzML files and csv files.')
    parser.add_argument('--mzml_files', nargs='+', help='List of mzML files')
    parser.add_argument('--csv_files', nargs='+', help='List of CSV files')
    parser.add_argument('--data_collector', type=str, default='Minions', help='Data collector.')
    parser.add_argument('--pi_name', type=str, default='Pieter Dorrestein', help='PI name.')
    parser.add_argument('--mass_detect_int_tol', type=float, default=5e4, help='Mass detection intensity tolerance.')
    parser.add_argument('--min_feature_height', type=float, default=1.5e5, help='Minimum feature height.')
    parser.add_argument('--mz_tol_ppm', type=float, default=10, help='m/z tolerance in ppm.')
    parser.add_argument('--ms2_explanation_cutoff', type=float, default=0.60, help='MS2 explanation cutoff.')
    parser.add_argument('--core_adduct_filter', type=str, default='simple',
                        help='Core adduct filter. Available options: none, full, simple.')
    parser.add_argument('--adduct_type_mode', type=str, default='full', help='Adduct type mode.')
    parser.add_argument('--plot', action='store_true', help='Plot the results.')
    parser.add_argument('--component_precursor_check', type=str, default='0', help='Combinatorial synthesis.')
    parser.add_argument('--preprocessed_pkl', type=str, default=None, help='Preprocessed pkl file.')
    args = parser.parse_args()

    main_batch(args.mzml_files, args.csv_files,
               data_collector=args.data_collector,
               pi_name=args.pi_name,
               mz_tol_ppm=args.mz_tol_ppm,
               mass_detect_int_tol=args.mass_detect_int_tol,
               min_feature_height=args.min_feature_height,
               ms2_explanation_cutoff=args.ms2_explanation_cutoff,
               core_adduct_filter=args.core_adduct_filter,
               adduct_type_mode=args.adduct_type_mode,
               plot=False,
               component_precursor_check=True if args.component_precursor_check == '1' else False,
               preprocessed_pkl_path=args.preprocessed_pkl)

    #############################################################################################################

    # main_batch([
    #     '../test_data/AP_176.mzML',
    #     # '../test_data/AP_177.mzML',
    #     '../test_data/AP_178.mzML'
    # ],
    #            ['../test_data/AP_176_178.csv'],
    #            adduct_type_mode='full',
    #            component_precursor_check=False,
    #            ms2_explanation_cutoff=0.60,
    #            core_adduct_filter='full',
    #            preprocessed_pkl_path='cmpd_name_to_mass.pkl')

    # main_batch([
    #     # '../test_data/AP_175.mzML',
    #             '../test_data/AP_176.mzML', '../test_data/AP_177.mzML',
    #             '../test_data/AP_178.mzML', '../test_data/AP_179.mzML', '../test_data/AP_180.mzML',
    #             '../test_data/AP_181.mzML', '../test_data/AP_182.mzML', '../test_data/AP_183.mzML',
    #             '../test_data/AP_184.mzML', '../test_data/AP_185.mzML', '../test_data/AP_186.mzML',
    #             '../test_data/AP_187.mzML', '../test_data/AP_188.mzML', '../test_data/AP_189.mzML',
    #             '../test_data/AP_190.mzML', '../test_data/AP_191.mzML', '../test_data/AP_192.mzML',
    #             '../test_data/AP_193.mzML', '../test_data/AP_194.mzML', '../test_data/AP_195.mzML',
    #             '../test_data/AP_196.mzML', '../test_data/AP_197.mzML', '../test_data/AP_198.mzML',
    #             '../test_data/AP_199.mzML', '../test_data/AP_200.mzML', '../test_data/AP_201.mzML',
    #             '../test_data/AP_202.mzML', '../test_data/AP_203.mzML', '../test_data/AP_204.mzML',
    #             '../test_data/AP_205.mzML', '../test_data/AP_206.mzML', '../test_data/AP_207.mzML',
    #             '../test_data/AP_208.mzML', '../test_data/AP_209.mzML', '../test_data/AP_210.mzML',
    #             '../test_data/AP_211.mzML', '../test_data/AP_212.mzML', '../test_data/AP_213.mzML',
    #             '../test_data/AP_214.mzML', '../test_data/AP_215.mzML', '../test_data/AP_216.mzML',
    #             '../test_data/AP_217.mzML', '../test_data/AP_218.mzML', '../test_data/AP_219.mzML',
    #             '../test_data/AP_220.mzML', '../test_data/AP_221.mzML', '../test_data/AP_222.mzML',
    #             '../test_data/AP_223.mzML', '../test_data/AP_224.mzML', '../test_data/AP_225.mzML',
    #             '../test_data/AP_226.mzML', '../test_data/AP_227.mzML', '../test_data/AP_228.mzML',
    #             '../test_data/AP_229.mzML', '../test_data/AP_230.mzML'],
    #            ['../test_data/AP_175_230_AA_1_modified_output_file.csv'],
    #            adduct_type_mode='full',
    #            component_precursor_check=True,
    #            ms2_explanation_cutoff=0.60,
    #            core_adduct_filter='full',
    #            preprocessed_pkl_path='cmpd_name_to_mass.pkl')

    # main_batch(['../test/VD_52.mzML'], ['../test/3_OH_VD_KV_saturated.csv'],
    #            adduct_type_mode='full',
    #            core_adduct_filter='none',
    #            component_precursor_check=False,
    #            preprocessed_pkl_path='cmpd_name_to_mass.pkl')
    # main_batch(['../test/P1_A1_510.mzML'],
    #            ['../test/PCP.csv'],
    #            adduct_type_mode='full',
    #            core_adduct_filter='none', component_precursor_check=False)
    # main_batch(['../test/reframe_drugs_pos_P1_A10_id.mzML'],
    #            ['../test/20241017_reframe_metadata_pos_gnps2_workflow.csv'],
    #            adduct_type_mode='full',
    #            core_adduct_filter='none')
