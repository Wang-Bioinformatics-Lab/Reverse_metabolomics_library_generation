import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math


def append_file_summary(all_file_summary_rows, mzml_name, cmpd_df, feature_df, metadata_df):
    total_features = None
    total_features_with_ms2 = None

    unique_target_cmpds = None

    unique_cmpds_with_ms1 = None
    unique_cmpds_with_ms1_M_H = None

    unique_cmpds_with_ms1_ms2 = None
    unique_cmpds_with_ms1_ms2_M_H = None

    unique_cmpds_with_ms1_ms2_selected = None
    unique_cmpds_with_ms1_ms2_selected_M_H = None

    annotation_distribution = None

    if cmpd_df is not None and not cmpd_df.empty:
        unique_target_cmpds = cmpd_df['inchi'].nunique()

    if feature_df is not None and not feature_df.empty:
        total_features = len(feature_df)
        total_features_with_ms2 = len(feature_df[feature_df['best_MS2_scan_idx'].notnull()])

    if metadata_df is not None and not metadata_df.empty:

        try:
            unique_cmpds_with_ms1 = int(metadata_df['inchi'].nunique())
            unique_cmpds_with_ms1_M_H = int(metadata_df[metadata_df['t_adduct'].isin(['[M+H]+', '[M-H]-'])]['inchi'].nunique())

            df = metadata_df[metadata_df['best_MS2_scan_idx'].notnull()].reset_index(drop=True)
            unique_cmpds_with_ms1_ms2 = int(df['inchi'].nunique())
            unique_cmpds_with_ms1_ms2_M_H = int(df[df['t_adduct'].isin(['[M+H]+', '[M-H]-'])]['inchi'].nunique())

            # only keep selected
            df = df[df['selected']].reset_index(drop=True)
            unique_cmpds_with_ms1_ms2_selected = int(df['inchi'].nunique())
            unique_cmpds_with_ms1_ms2_selected_M_H = int(df[df['t_adduct'].isin(['[M+H]+', '[M-H]-'])]['inchi'].nunique())

            # count # of annotations for each MS2 (ie, how many MS2 are annotated with one compound, how many are annotated with two compounds, etc.)
            annotation_counts = df.groupby('best_MS2_scan_idx').size().value_counts().sort_index()

            # Calculate percentages
            total_spectra = annotation_counts.sum()
            annotation_distribution = {
                n_compounds: [count, round(count / total_spectra * 100, 2)]
                for n_compounds, count in annotation_counts.items()
            }

        except:
            pass

    all_file_summary_rows.append({
        'mzml_name': mzml_name,
        'total_features': total_features,
        'total_features_with_ms2': total_features_with_ms2,
        'unique_target_cmpds': unique_target_cmpds,
        'unique_cmpds_with_ms1': unique_cmpds_with_ms1,
        'unique_cmpds_with_ms1_M_H': unique_cmpds_with_ms1_M_H,
        'unique_cmpds_with_ms1_ms2': unique_cmpds_with_ms1_ms2,
        'unique_cmpds_with_ms1_ms2_M_H': unique_cmpds_with_ms1_ms2_M_H,
        'unique_cmpds_with_ms1_ms2_selected': unique_cmpds_with_ms1_ms2_selected,
        'unique_cmpds_with_ms1_ms2_selected_M_H': unique_cmpds_with_ms1_ms2_selected_M_H,
        'ms2_annotation_number_distribution': annotation_distribution
    })

    return all_file_summary_rows


def plot_ms2_annotation_distribution(file_summary_df):
    if file_summary_df is None or file_summary_df.empty:
        return

    try:
        # Get max annotation count across all files
        max_annotation = 1  # default to at least 1
        for dist in file_summary_df['ms2_annotation_number_distribution'].dropna():
            if dist:  # check if not None
                max_annotation = max(max_annotation, max(dist.keys()))

        # Calculate number of pages needed (4 plots per page)
        n_files = len(file_summary_df)
        plots_per_page = 4
        n_pages = math.ceil(n_files / plots_per_page)

        # Create PDF
        with PdfPages('ms2_annotation_distributions.pdf') as pdf:
            for page in range(n_pages):
                # Create figure with 4 subplots
                fig, axes = plt.subplots(plots_per_page, 1, figsize=(8, 16))
                # fig.suptitle('MS2 Annotation Distribution', fontsize=16, y=0.95)

                # Process each file for this page
                start_idx = page * plots_per_page
                end_idx = min(start_idx + plots_per_page, n_files)

                for i, (idx, row) in enumerate(file_summary_df.iloc[start_idx:end_idx].iterrows()):
                    ax = axes[i]

                    # Get distribution dictionary
                    dist = row['ms2_annotation_number_distribution']

                    if dist is not None:
                        # Create x-axis values (1 to max_annotation)
                        x = list(range(1, max_annotation + 1))

                        # Create y-axis values (counts)
                        y = [dist.get(j, [0, 0])[0] for j in x]  # get count or 0 if missing

                        # Create bar plot
                        bars = ax.bar(x, y)

                        # Customize plot
                        ax.set_xlabel('Number of annotations per MS2')
                        ax.set_ylabel('Count of MS2 spectra')
                        ax.set_title(f'File: {row["mzml_name"]}')
                        ax.grid(True, alpha=0.3)

                        # Set x-axis ticks to integers only
                        ax.set_xticks(x)

                        # Add value labels on top of bars
                        for j, bar in enumerate(bars):
                            height = bar.get_height()
                            if height > 0:  # Only add label if there are spectra
                                # Get percentage from distribution dictionary
                                percent = dist.get(j + 1, [0, 0])[1]  # j+1 because x starts from 1
                                ax.text(bar.get_x() + bar.get_width() / 2., height,
                                        f'{int(height)} ({percent}%)',
                                        ha='center', va='bottom')
                    else:
                        ax.text(0.5, 0.5, 'No annotation data available',
                                ha='center', va='center',
                                transform=ax.transAxes)

                    # Remove empty subplot if we don't have data for all 4 spots
                    if i >= end_idx - start_idx:
                        ax.remove()

                # Adjust layout and save page
                plt.tight_layout()
                pdf.savefig()
                plt.close()

        print('Saved ms2_annotation_distributions.pdf')

    except:
        pass
