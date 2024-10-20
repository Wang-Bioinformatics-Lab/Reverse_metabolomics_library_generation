# from .basic_filter import remove_smiles_with_empty_valid_ms2, remove_doubly_charged_ions, remove_isotopes
from .ms2_explanation_filter import filter_by_ms2_explanation
from .core_adduct_filter import filter_by_core_adduct


def filter_feature_df(feature_df, intensity_threshold):
    """
    Filter the feature DataFrame based on intensity threshold
    """

    # filter by ROI length
    feature_df = feature_df[feature_df['length'] >= 4].reset_index(drop=True)

    # filter by intensity
    feature_df = feature_df[feature_df['peak_height'] >= intensity_threshold].reset_index(drop=True)

    return feature_df


def filter_df(df,
              ion_mode,
              ms2_explanation_cutoff=0.60,
              core_adduct_filter=True):
    """
    Filter the merged DataFrame
    """

    # # if matched to a doubly charged ion, remove the match
    # df = df.apply(remove_doubly_charged_ions, axis=1)

    # # if matched to an isotope, remove the match
    # df = df.apply(remove_isotopes, axis=1)

    # filter by MS2 explained intensity
    if ms2_explanation_cutoff > 0.0:
        print('Calculating MS2 explanation intensity...')
        df = df.apply(lambda row: filter_by_ms2_explanation(row, explanation_cutoff=ms2_explanation_cutoff), axis=1)

    # # remove smiles with 0 valid MS2
    # df = remove_smiles_with_empty_valid_ms2(df)

    # filter by core adducts
    if core_adduct_filter:
        print('Filtering by core adducts...')
        df = filter_by_core_adduct(df, ion_mode, core_adduct_ls=None, rt_tol=0.05)

    return df


