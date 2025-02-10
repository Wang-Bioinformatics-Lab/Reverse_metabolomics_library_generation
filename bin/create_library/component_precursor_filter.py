import pickle
import numpy as np


def filter_by_component_precursor(df, ion_mode, preprocessed_pkl_path=None):
    """
    Filter by component precursor
    """

    with open(preprocessed_pkl_path, 'rb') as f:
        cmpd_name_to_mass = pickle.load(f)

    for i, row in df.iterrows():

        if not '_' in row['compound_name']:
            continue

        if not (row['selected'] and row['MS2'] is not None):
            continue

        # split the compound name by '_'
        cmpd_name_ls = row['compound_name'].split('_')

        mass_ls = []
        valid_compounds = True  # Flag to track if all compounds are valid

        for cmpd_name in cmpd_name_ls:
            if cmpd_name not in cmpd_name_to_mass:
                valid_compounds = False
                break  # Exit the inner loop

            this_mass = cmpd_name_to_mass[cmpd_name]
            mass_ls.append(this_mass)
            mass_ls.append(this_mass - 18.010565)

        if not valid_compounds:
            df.at[i, 'selected'] = False

            if df.at[i, 'discard_reason'] == '':
                df.at[i, 'discard_reason'] = 'Component not found in DB'
            else:
                df.at[i, 'discard_reason'] += '; Component not found in DB'
            continue  # Skip to next row in outer loop

        if len(mass_ls) == 0:
            continue

        # if more than 2 components, add the mass of all combinations
        if len(mass_ls) > 2:
            for m in range(len(mass_ls)):
                for n in range(m + 1, len(mass_ls)):
                    mass_ls.append(mass_ls[m] + mass_ls[n] - 18.010565)
                    mass_ls.append(mass_ls[m] + mass_ls[n] - 18.010565 * 2)

        if ion_mode == 'pos':
            frag_mz_ls = [mass + 1.007276 for mass in mass_ls]
            if row['t_adduct'] == '[M+Na]+':
                frag_mz_ls.extend([mass + 22.98922 for mass in mass_ls])
            if row['t_adduct'] == '[M+K]+':
                frag_mz_ls.extend([mass + 38.96316 for mass in mass_ls])
        else:
            frag_mz_ls = [mass - 1.007276 for mass in mass_ls]

        ms2_peaks = row['MS2']
        # minimum 1% intensity
        ms2_peaks = ms2_peaks[ms2_peaks[:, 1] > 0.01 * np.max(ms2_peaks[:, 1])]
        ms2_mzs = ms2_peaks[:, 0]

        if not np.any(np.abs(np.subtract.outer(frag_mz_ls, ms2_mzs)) < 0.02):
            df.at[i, 'selected'] = False

            if df.at[i, 'discard_reason'] == '':
                df.at[i, 'discard_reason'] = 'No any component precursor'
            else:
                df.at[i, 'discard_reason'] += '; No any component precursor'

    return df
