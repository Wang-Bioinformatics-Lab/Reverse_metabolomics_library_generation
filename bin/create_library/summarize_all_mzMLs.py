import pandas as pd


def append_file_summary(all_file_summary_rows, mzml_name, cmpd_df, feature_df, metadata_df):
    total_features = None
    total_features_with_ms2 = None

    unique_cmpds_in_target_list = None

    unique_ms2_with_mz_matched_to_target_list = None
    unique_ms2_annotated = None
    unique_cmpds_annotated = None
    unique_cmpd_adduct_pairs_annotated = None

    try:
        if cmpd_df is not None and not cmpd_df.empty:
            unique_cmpds_in_target_list = cmpd_df['inchi'].nunique()

        if feature_df is not None and not feature_df.empty:
            total_features = len(feature_df)
            total_features_with_ms2 = len(feature_df[feature_df['best_MS2_scan_idx'].notnull()])

        if metadata_df is not None and not metadata_df.empty:
            unique_ms2 = metadata_df['best_MS2_scan_idx'].unique()
            unique_ms2 = unique_ms2[~pd.isnull(unique_ms2)]
            unique_ms2_with_mz_matched_to_target_list = len(unique_ms2)

            # only keep selected
            metadata_df = metadata_df[metadata_df['selected']].reset_index(drop=True)
            unique_ms2 = metadata_df['best_MS2_scan_idx'].unique()
            unique_ms2 = unique_ms2[~pd.isnull(unique_ms2)]
            unique_ms2_annotated = len(unique_ms2)

            unique_cmpds_annotated = metadata_df['inchi'].nunique()

            unique_cmpd_adduct_pairs_annotated = metadata_df['inchi_adduct'].nunique()
    except:
        pass

    all_file_summary_rows.append({
        'mzml_name': mzml_name,
        'total_features': total_features,
        'total_features_with_ms2': total_features_with_ms2,
        'unique_cmpds_in_target_list': unique_cmpds_in_target_list,
        'unique_ms2_with_mz_matched_to_target_list': unique_ms2_with_mz_matched_to_target_list,
        'unique_ms2_annotated': unique_ms2_annotated,
        'unique_cmpds_annotated': unique_cmpds_annotated,
        'unique_cmpd_adduct_pairs_annotated': unique_cmpd_adduct_pairs_annotated
    })

    return all_file_summary_rows
