import pandas as pd
from .utils import neutralize_formula, calc_exact_mass, smiles_to_formula_and_inchi


adduct_pos_simple = [
    {'name': '[M+H]+', 'm': 1, 'charge': 1, 'mass': 1.00727645223},
    {'name': '[M+Na]+', 'm': 1, 'charge': 1, 'mass': 22.989220702},
    {'name': '[M+NH4]+', 'm': 1, 'charge': 1, 'mass': 18.03382555335},
    # {'name': '[M+K]+', 'm': 1, 'charge': 1, 'mass': 38.9631579064}
]

adduct_pos_full = [
    {'name': '[M+H]+', 'm': 1, 'charge': 1, 'mass': 1.00727645223},
    {'name': '[M+Na]+', 'm': 1, 'charge': 1, 'mass': 22.989220702},
    {'name': '[M+K]+', 'm': 1, 'charge': 1, 'mass': 38.9631579064},
    {'name': '[M+NH4]+', 'm': 1, 'charge': 1, 'mass': 18.03382555335},
    # {'name': '[M-H2O+NH4]+', 'm': 1, 'charge': 1, 'mass': 0.023260869},
    {'name': '[M+H-H2O]+', 'm': 1, 'charge': 1, 'mass': -17.0032882318},
    {'name': '[M+H-2H2O]+', 'm': 1, 'charge': 1, 'mass': -35.01385291583},
    # {'name': '[M+H-3H2O]+', 'm': 1, 'charge': 1, 'mass': -53.02441759986}
]

adduct_neg_simple = [
    {'name': '[M-H]-', 'm': 1, 'charge': 1, 'mass': -1.00727645223}
]

adduct_neg_full = [
    {'name': '[M-H]-', 'm': 1, 'charge': 1, 'mass': -1.00727645223},
    {'name': '[M+Cl]-', 'm': 1, 'charge': 1, 'mass': 34.968304102},
    # {'name': '[M+Br]-', 'm': 1, 'charge': 1, 'mass': 78.91778902},
    {'name': '[M+FA]-', 'm': 1, 'charge': 1, 'mass': 44.99710569137},
    {'name': '[M+Ac]-', 'm': 1, 'charge': 1, 'mass': 59.01275575583},
    {'name': '[M-H-H2O]-', 'm': 1, 'charge': 1, 'mass': -19.01784113626}
]


def prepare_cmpd_df(cmpd_df, adduct_type_mode):

    """
    Calculate the exact mass for each compound in the compound list
    """

    if adduct_type_mode == 'simple':
        adduct_pos = adduct_pos_simple
        adduct_neg = adduct_neg_simple
    else:
        adduct_pos = adduct_pos_full
        adduct_neg = adduct_neg_full

    # drop rows w/o SMILES
    cmpd_df = cmpd_df.dropna(subset=['unique_sample_id', 'SMILES']).reset_index(drop=True)
    # unique_sample_id should contain '.mz'
    cmpd_df = cmpd_df[cmpd_df['unique_sample_id'].str.contains('.mz')].reset_index(drop=True)

    # Process unique SMILES to get formula and InChI
    unique_smiles = cmpd_df['SMILES'].unique()
    smiles_to_formula_inchi = {}

    for smiles in unique_smiles:
        formula, inchi = smiles_to_formula_and_inchi(smiles)
        smiles_to_formula_inchi[smiles] = (formula, inchi)

    # Apply to original dataframe
    cmpd_df['formula'] = cmpd_df['SMILES'].map(lambda x: smiles_to_formula_inchi[x][0])
    cmpd_df['inchi'] = cmpd_df['SMILES'].map(lambda x: smiles_to_formula_inchi[x][1])
    del unique_smiles, smiles_to_formula_inchi

    # dereplicate by unique_sample_id, inchi (same reaction, same compound)
    cmpd_df = cmpd_df.drop_duplicates(subset=['unique_sample_id', 'inchi'], keep='first').reset_index(drop=True)

    # Process unique formulas to get neutralized formulas
    unique_formulas = cmpd_df['formula'].unique()
    formula_to_neutralized = {}

    for formula in unique_formulas:
        neutralized = neutralize_formula(formula)
        formula_to_neutralized[formula] = neutralized

    # Apply to original dataframe
    cmpd_df['neutralized_formula'] = cmpd_df['formula'].map(lambda x: formula_to_neutralized[x])
    del unique_formulas, formula_to_neutralized

    # calculate the exact mass
    cmpd_df['exact_mass'] = cmpd_df['neutralized_formula'].apply(calc_exact_mass)

    # create a list to store the new rows
    new_rows = []

    # Iterate through each compound and adduct
    for _, compound in cmpd_df.iterrows():
        for adduct in adduct_pos:
            mz = compound['exact_mass'] + adduct['mass']
            new_row = compound.to_dict()
            new_row['t_mz'] = mz
            new_row['t_adduct'] = adduct['name']
            new_row['ion_mode'] = 'positive'
            new_rows.append(new_row)

        for adduct in adduct_neg:
            mz = compound['exact_mass'] + adduct['mass']
            new_row = compound.to_dict()
            new_row['t_mz'] = mz
            new_row['t_adduct'] = adduct['name']
            new_row['ion_mode'] = 'negative'
            new_rows.append(new_row)

    # Create a new DataFrame from the list of new rows
    df = pd.DataFrame(new_rows)

    return df
