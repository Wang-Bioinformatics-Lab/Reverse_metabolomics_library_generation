import pickle
import numpy as np


def filter_by_component_precursor(df, ion_mode, preprocessed_pkl_path=None):
    """
    Filter by component precursor
    """

    with open(preprocessed_pkl_path, 'rb') as f:
        cmpd_name_to_mass = pickle.load(f)

    for i, row in df.iterrows():
        if not row['selected']:
            continue

        if not '_' in row['compound_name']:
            continue

        cmpd_name_ls = row['compound_name'].split('_')

        mass_ls = []
        for cmpd_name in cmpd_name_ls:
            if cmpd_name in cmpd_name_to_mass:
                this_mass = cmpd_name_to_mass[cmpd_name]
                mass_ls.append(this_mass)

        if len(mass_ls) == 0:
            continue

        # if more than 2 components, add the mass of all combinations
        if len(mass_ls) > 2:
            for m in range(len(mass_ls)):
                for n in range(m + 1, len(mass_ls)):
                    mass_ls.append(mass_ls[m] + mass_ls[n] - 18.010565)

        if ion_mode == 'pos':
            frag_mz_ls = [mass + 1.007276 for mass in mass_ls]
            if row['t_adduct'] == '[M+Na]+':
                frag_mz_ls.extend([mass + 22.98922 for mass in mass_ls])
        else:
            frag_mz_ls = [mass - 1.007276 for mass in mass_ls]

        ms2_mzs = row['MS2'][:, 0]
        if not np.any(np.abs(np.subtract.outer(frag_mz_ls, ms2_mzs)) < 0.02):
            df.at[i, 'selected'] = False

            if df.at[i, 'discard_reason'] == '':
                df.at[i, 'discard_reason'] = 'No any component precursor'
            else:
                df.at[i, 'discard_reason'] += '; No any component precursor'

    return df
