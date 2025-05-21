import pandas as pd
import os


def write_mgf(df, file_name, scans_start=0):
    """
    Write the filtered library to a file.
    """

    df = df[(df['selected']) & (~pd.isnull(df['best_MS2_scan_idx']))].reset_index(drop=True)

    out_path = 'all_ms2.mgf'

    with open(out_path, 'a', encoding='utf-8') as f:
        for _, row in df.iterrows():

            ms2_scan = int(row['best_MS2_scan_idx'])

            f.write('BEGIN IONS\n')
            f.write(f'NAME={row["name"]}\n')
            f.write(f'PEPMASS={row["t_mz"]}\n')
            # f.write(f'MSLEVEL=2\n')
            f.write(f'TITLE={file_name}:scan:{ms2_scan}\n')
            f.write(f'SMILES={row["SMILES"]}\n')
            f.write(f'INCHI={row["inchi"]}\n')
            # inchi_aux = row['isomer_inchis'] if row['isomer_inchis'] != '' else 'N/A'
            # f.write(f'INCHIAUX={inchi_aux}\n')
            f.write(f'ADDUCT={row["t_adduct"]}\n')
            f.write(f'SCANS={scans_start}\n')

            mzs = row['MS2'][:, 0]
            intensities = row['MS2'][:, 1]
            for mz, intensity in zip(mzs, intensities):
                mz = round(mz, 5)
                intensity = round(intensity, 4)
                f.write(f'{mz} {intensity}\n')
            f.write('END IONS\n\n')
            scans_start += 1

    return scans_start
