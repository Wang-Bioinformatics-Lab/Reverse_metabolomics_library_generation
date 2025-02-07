# from .basic_filter import remove_smiles_with_empty_valid_ms2, remove_doubly_charged_ions, remove_isotopes
from .ms2_explanation_filter import filter_by_ms2_explanation
from .core_adduct_filter import filter_by_core_adduct
from .component_precursor_filter import filter_by_component_precursor


def filter_df(df,
              ion_mode,
              ms2_explanation_cutoff=0.60,
              core_adduct_filter='simple',
              component_precursor_check=True, preprocessed_pkl_path=None):
    """
    Filter the merged DataFrame

    core_adduct_filter: str, one of ['full', 'simple', 'none']
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

    # filter by component precursor
    if component_precursor_check:
        print('Filtering by component precursor...')
        df = filter_by_component_precursor(df, ion_mode, preprocessed_pkl_path=preprocessed_pkl_path)

    # filter by core adducts
    if core_adduct_filter != 'none':
        print('Filtering by core adducts...')
        df = filter_by_core_adduct(df, ion_mode, core_adduct_filter_mode=core_adduct_filter, rt_tol=0.05)

    return df


