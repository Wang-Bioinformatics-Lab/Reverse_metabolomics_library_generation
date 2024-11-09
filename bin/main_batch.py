import os
import argparse
import pandas as pd

from feature_extraction import feature_extraction_single, plot_all_ms2, plot_all_eic, plot_mz_rt
from cmpd import prepare_cmpd_df
from create_library import create_library, append_file_summary, plot_ms2_annotation_distribution


def main_batch(mzml_files, csv_files,
               data_collector='Minions',
               pi_name='Pieter Dorrestein',
               mass_detect_int_tol=5e4,
               min_feature_height=1.5e5,
               mz_tol_ppm=10,
               ms2_explanation_cutoff=0.60,
               core_adduct_filter='full',
               adduct_type_mode='full',
               plot=False,
               write_individual_mgf=False):
    """
    Process a batch of mzML files and csv files.
    """

    # Load the csv file
    print('Calculating compound masses...')
    all_cmpd_df = pd.DataFrame()
    for csv in csv_files:
        print('Loading', csv)
        cmpd_df = prepare_cmpd_df(csv, adduct_type_mode)
        all_cmpd_df = pd.concat([all_cmpd_df, cmpd_df])

    # Get unique mzmls
    unique_mzmls = all_cmpd_df['unique_sample_id'].unique()

    # Initialize the library dataframe, for uploading to GNPS
    all_library_df = pd.DataFrame()

    # file summary rows
    all_file_summary_rows = []

    # Load the mzml file
    for mzml in mzml_files:

        print('Processing', mzml)
        mzml_name = os.path.basename(mzml)
        mzml_basename = mzml_name.split('.mz')[0]

        if mzml_name not in unique_mzmls:
            print(f'{mzml_name} not in the csv file. Skipping...')

        # Extract features from mzML file
        try:
            feature_df, ion_mode = feature_extraction_single(file_path=mzml,
                                                             mass_detect_int_tol=mass_detect_int_tol,
                                                             min_feature_height=min_feature_height,
                                                             save=False)
        except:
            print(f'Error extracting features from {mzml}. Skipping...')
            continue

        # Load the compound list for the mzml file
        cmpd_df = all_cmpd_df[(all_cmpd_df['unique_sample_id'] == mzml_name) &
                              (all_cmpd_df['ion_mode'] == ion_mode)].copy().reset_index(drop=True)

        # Filter library
        print('Creating MS/MS library...')
        df, library_df = create_library(cmpd_df, feature_df, ion_mode,
                                        data_collector, pi_name, mzml_basename,
                                        mz_tol_ppm=mz_tol_ppm,
                                        filter_library=True,
                                        ms2_explanation_cutoff=ms2_explanation_cutoff,
                                        core_adduct_filter=core_adduct_filter,
                                        metadata_dir=None,
                                        write_individual_mgf=write_individual_mgf)

        # Append file summary rows
        all_file_summary_rows = append_file_summary(all_file_summary_rows, mzml_name, cmpd_df, feature_df, df)

        if library_df is not None:
            all_library_df = pd.concat([all_library_df, library_df])

        # Save
        if df is not None:
            df.to_csv(f'{mzml_basename}_metadata.tsv', sep='\t', index=False)

        if plot and df is not None and feature_df is not None:
            try:
                # Plot all MS2 spectra
                print('Plotting all spectra...')
                plot_all_ms2(df, mzml)

                # Plot all EICs
                print('Plotting all EICs...')
                plot_all_eic(df, mzml)

                # Plot mz-rt scatter plot
                print('Plotting mz-rt scatter plot...')
                plot_mz_rt(feature_df, df, mzml_basename)
            except Exception as e:
                print(e)

    # Save the library dataframe
    all_library_df.to_csv('all_library.tsv', sep='\t', index=False, na_rep='N/A')

    # Save the file summary
    file_summary_df = pd.DataFrame(all_file_summary_rows)
    file_summary_df.to_csv('file_summary.tsv', sep='\t', index=False)

    # Plot MS2 annotation distribution
    plot_ms2_annotation_distribution(file_summary_df)


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
    parser.add_argument('--core_adduct_filter', type=str, default='full', help='Core adduct filter. Available options: none, full, simple.')
    parser.add_argument('--adduct_type_mode', type=str, default='full', help='Adduct type mode.')
    parser.add_argument('--plot', action='store_true', help='Plot the results.')
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
               plot=args.plot)

    ##############################################################################################################
    # main_batch(['../test/VD_52.mzML'], ['../test/3_OH_VD_KV_saturated.csv'],
    #            ms2_explanation_cutoff=0.0,
    #            adduct_type_mode='full',
    #            core_adduct_filter='none')
    # main_batch(['../test/P1_A1_510.mzML'], ['../test/PCP.csv'],
    #            ms2_explanation_cutoff=0.0,
    #            adduct_type_mode='full',
    #            core_adduct_filter='none')
    # main_batch(['../test/reframe_drugs_pos_P1_A10_id.mzML'],
    #            ['../test/20241017_reframe_metadata_pos_gnps2_workflow.csv'],
    #            ms2_explanation_cutoff=0.0,
    #            adduct_type_mode='full',
    #            core_adduct_filter='none')
