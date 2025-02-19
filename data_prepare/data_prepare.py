"""
preprocess all reactants tsv file
"""

from _utils import neutralize_formula, calc_exact_mass, smiles_to_formula, chloride_to_acid
import pandas as pd
import pickle


def main():
    df = pd.read_csv('all_reactants.tsv', sep='\t', low_memory=False)

    # drop rows w/o SMILES
    df = df.dropna(subset=['SMILES']).reset_index(drop=True)

    # convert compound_name to lower case
    df['compound_name'] = df['compound_name'].str.lower()

    # dereplicate by compound name
    df = df.drop_duplicates(subset=['compound_name']).reset_index(drop=True)

    # Add formula and InChI information
    df['corrected_SMILES'] = df['SMILES'].apply(chloride_to_acid)

    # Add formula and InChI information
    df['formula'] = df['corrected_SMILES'].apply(smiles_to_formula)

    # neutralize the formula, deal with the charge (e.g., C5H5N+)
    df['neutralized_formula'] = df['formula'].apply(neutralize_formula)

    # calculate the exact mass
    df['exact_mass'] = df['neutralized_formula'].apply(calc_exact_mass)

    print(f'Number of unique reactants: {len(df)}')

    # save
    df.to_csv('all_reactants_preprocessed.tsv', sep='\t', index=False)

    # remove rows with missing exact_mass
    df = df.dropna(subset=['exact_mass']).reset_index(drop=True)

    # save dict from compound name to exact_mass
    compound_name_to_mass = dict(zip(df['compound_name'], df['exact_mass']))

    with open('../bin/cmpd_name_to_mass.pkl', 'wb') as f:
        pickle.dump(compound_name_to_mass, f)


if __name__ == '__main__':
    main()
