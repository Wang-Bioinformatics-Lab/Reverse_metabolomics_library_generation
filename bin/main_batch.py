import os
import argparse
import pandas as pd

from feature_extraction import feature_extraction_single, plot_all_ms2, plot_all_eic, plot_mz_rt
from cmpd import prepare_cmpd_df
from create_library import create_library


def main_batch(mzml_files, csv_files,
               data_collector='Minions',
               ms2_explanation_cutoff=0.60,
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
        cmpd_df = prepare_cmpd_df(csv)
        all_cmpd_df = pd.concat([all_cmpd_df, cmpd_df])

    # Get unique mzmls
    unique_mzmls = all_cmpd_df['unique_sample_id'].unique()

    # Load the mzml file
    all_library_df = pd.DataFrame()
    for mzml in mzml_files:

        print('Processing', mzml)
        mzml_name = os.path.basename(mzml).split('.mz')[0]

        if mzml_name + '.mzML' not in unique_mzmls:
            print(f'{mzml_name}.mzML not in the csv file. Skipping...')

        # Extract features from mzML file
        feature_df, ion_mode, intensity_threshold = feature_extraction_single(file_path=mzml, save=False)

        # Load the compound list for the mzml file
        cmpd_df = all_cmpd_df[(all_cmpd_df['unique_sample_id'] == mzml_name + '.mzML') &
                              (all_cmpd_df['ion_mode'] == ion_mode)].reset_index(drop=True).copy()

        # Filter library
        print('Creating MS/MS library...')
        df, library_df = create_library(cmpd_df, feature_df, ion_mode, intensity_threshold,
                                        data_collector, mzml_name,
                                        filter_library=True, ms2_explanation_cutoff=ms2_explanation_cutoff,
                                        metadata_dir=None,
                                        write_individual_mgf=write_individual_mgf)

        all_library_df = pd.concat([all_library_df, library_df])

        # Save the filtered library
        df.to_csv(f'{mzml_name}_metadata.tsv', sep='\t', index=False)

        if plot:
            # Plot all MS2 spectra
            print('Plotting all MS2 spectra...')
            plot_all_ms2(df, mzml)

            # Plot all EICs
            print('Plotting all EICs...')
            plot_all_eic(df, mzml)

            # Plot mz-rt scatter plot
            print('Plotting mz-rt scatter plot...')
            plot_mz_rt(feature_df, df, mzml_name)

    all_library_df.to_csv('all_library.tsv', sep='\t', index=False, na_rep='N/A')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a batch of mzML files and csv files.')
    parser.add_argument('--mzml_files', nargs='+', help='List of mzML files')
    parser.add_argument('--csv_files', nargs='+', help='List of CSV files')
    parser.add_argument('--data_collector', type=str, default='Minions', help='Data collector.')
    parser.add_argument('--ms2_explanation_cutoff', type=float, default=0.60, help='MS2 explanation cutoff.')
    parser.add_argument('--plot', action='store_true', help='Plot the results.')
    args = parser.parse_args()

    main_batch(args.mzml_files, args.csv_files,
               data_collector=args.data_collector,
               ms2_explanation_cutoff=args.ms2_explanation_cutoff,
               plot=args.plot)
