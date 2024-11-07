import pandas as pd


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

    try:
        if cmpd_df is not None and not cmpd_df.empty:
            unique_target_cmpds = cmpd_df['inchi'].nunique()

        if feature_df is not None and not feature_df.empty:
            total_features = len(feature_df)
            total_features_with_ms2 = len(feature_df[feature_df['best_MS2_scan_idx'].notnull()])

        if metadata_df is not None and not metadata_df.empty:

            unique_cmpds_with_ms1 = metadata_df['inchi'].nunique()
            unique_cmpds_with_ms1_M_H = metadata_df[metadata_df['t_adduct'].isin(['[M+H]+', '[M-H]-'])]['inchi'].nunique()

            df = metadata_df[metadata_df['best_MS2_scan_idx'].notnull()].reset_index(drop=True)
            unique_cmpds_with_ms1_ms2 = df['inchi'].nunique()
            unique_cmpds_with_ms1_ms2_M_H = df[df['t_adduct'].isin(['[M+H]+', '[M-H]-'])]['inchi'].nunique()

            # only keep selected
            df = df[df['selected']].reset_index(drop=True)
            unique_cmpds_with_ms1_ms2_selected = df['inchi'].nunique()
            unique_cmpds_with_ms1_ms2_selected_M_H = df[df['t_adduct'].isin(['[M+H]+', '[M-H]-'])]['inchi'].nunique()
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
        'unique_cmpds_with_ms1_ms2_selected_M_H': unique_cmpds_with_ms1_ms2_selected_M_H
    })

    return all_file_summary_rows
