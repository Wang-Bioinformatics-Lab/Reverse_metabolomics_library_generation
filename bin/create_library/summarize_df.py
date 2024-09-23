import pandas as pd


def summarize_df(df):
    """
    Name: XXXX (known structural isomers: 2; isobaric peaks in run: 4)
    include isomers' inchis in a separate column
    """

    df = df[df['selected']].reset_index(drop=True)
    df['SMILES_adduct'] = df.apply(lambda x: str(x['SMILES']) + ' ' + str(x['t_adduct']), axis=1)
    df['inchi_adduct'] = df.apply(lambda x: str(x['inchi']) + ' ' + str(x['t_adduct']), axis=1)

    df['name'] = df['compound_name']

    # Count the number of isomers (how many unique SMILES_adduct are using the same MS2)
    df['isomer_count'] = 1
    df['isomer_inchis'] = None
    # Count the number of isobaric peaks in the run (how many features are associating the same SMILES_adduct)
    df['isobaric_peak_count'] = 1

    for i, row in df.iterrows():
        # the mask for the same MS2 scan
        mask = df['best_MS2_scan_idx'] == row['best_MS2_scan_idx']
        row['isomer_count'] = df.loc[mask, 'SMILES_adduct'].nunique()

        # for INCHI_AUX
        isomer_inchi_ls = df.loc[mask, 'inchi_adduct'].unique().tolist()
        # remove the current inchi_adduct
        if row['inchi_adduct'] in isomer_inchi_ls:
            isomer_inchi_ls.remove(row['inchi_adduct'])
        row['isomer_inchis'] = ';'.join(isomer_inchi_ls)

        mask2 = df['SMILES_adduct'] == row['SMILES_adduct']
        row['isobaric_peak_count'] = len(df.loc[mask2])

        row['name'] = f"{row['compound_name']} (known structural isomers: {row['isomer_count'] - 1}; isobaric peaks in run: {row['isobaric_peak_count'] - 1})"
        df.loc[i] = row

    return df

