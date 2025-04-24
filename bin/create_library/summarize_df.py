import pandas as pd


def summarize_df(df):
    """
    Name: XXXX (known structural isomers: 2; isobaric peaks in run: 4)
    include isomers' inchis in a separate column
    """

    df['inchi_adduct'] = df.apply(lambda x: str(x['inchi']) + ' ' + str(x['t_adduct']), axis=1)

    df['name'] = df['compound_name']

    # Count the number of isomers (how many unique inchi_adduct are using the same MS2)
    df['isomer_count'] = None
    df['isomer_inchis'] = None
    # Count the number of isobaric peaks in the run (how many features are associating the same inchi_adduct)
    df['isobaric_peak_count'] = None

    for i, row in df.iterrows():

        if not row['selected']:
            continue

        if pd.isnull(row['best_MS2_scan_idx']):
            continue

        # the mask for the same MS2 scan (excluding current row)
        mask = (pd.notnull(df['best_MS2_scan_idx'])) & \
               (df['best_MS2_scan_idx'] == row['best_MS2_scan_idx']) & \
               (df['selected']) & \
               (df.index != i)  # This excludes the current row

        isomer_inchi_ls = df.loc[mask, 'inchi_adduct'].unique().tolist()
        df.at[i, 'isomer_count'] = len(isomer_inchi_ls)
        df.at[i, 'isomer_inchis'] = ';'.join(isomer_inchi_ls)

        mask2 = df['inchi_adduct'] == row['inchi_adduct']
        isobaric_peak_count = len(df.loc[mask2]) - 1
        df.at[i, 'isobaric_peak_count'] = isobaric_peak_count

        new_name = f"{row['compound_name']} (known structural isomers: {len(isomer_inchi_ls)}; isobaric peaks in run: {isobaric_peak_count})"
        df.at[i, 'name'] = new_name

    ##############################
    # # one MS2 should be associated with one compound at most
    # # dereplicate by best_MS2_scan_idx if it is not null
    # df1 = df[pd.isnull(df['best_MS2_scan_idx'])].reset_index(drop=True)
    # df2 = df[pd.notnull(df['best_MS2_scan_idx'])].reset_index(drop=True)
    #
    # # Define the custom order for adducts
    # adduct_order = ['[M+H]+', '[M-H]-', '[M+H-H2O]+']  # The adducts not in this list will be sorted after these
    # # Convert t_adduct to Categorical with custom ordering
    # df2['t_adduct'] = pd.Categorical(df2['t_adduct'],
    #                                  categories=adduct_order + [x for x in df2['t_adduct'].unique() if
    #                                                             x not in adduct_order],
    #                                  ordered=True)
    #
    # # Sort by selected and then t_adduct
    # df2 = df2.sort_values(by=['selected', 't_adduct'], ascending=[False, True]).reset_index(drop=True)
    # df2 = df2.drop_duplicates(subset='best_MS2_scan_idx').reset_index(drop=True)
    # df = pd.concat([df1, df2], axis=0).reset_index(drop=True)

    return df

