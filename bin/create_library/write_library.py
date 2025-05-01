import pandas as pd


def write_library(df, data_collector, pi_name, file_name, ion_mode, scans_start=0):
    """
    Write the filtered library to a file.
    """

    df = df[(df['selected']) & (~pd.isnull(df['best_MS2_scan_idx']))].reset_index(drop=True)

    charge = 1 if ion_mode == 'positive' else -1
    _ion_mode = 'Positive' if ion_mode == 'positive' else 'Negative'

    rows = []
    scan = scans_start
    for _, row in df.iterrows():
        ms2_scan = int(row['best_MS2_scan_idx'])
        rows.append({
            'FILENAME': 'all_ms2.mgf',
            'SEQ': '*..*',
            'COMPOUND_NAME': row['name'],
            'MOLECULEMASS': row['MS2_precursor_mz'],
            'INSTRUMENT': 'Orbitrap',
            'IONSOURCE': 'LC-ESI',
            'EXTRACTSCAN': scan,
            'SMILES': row['SMILES'],
            'INCHI': row['inchi'],
            'INCHIAUX': row['isomer_inchis'],
            'CHARGE': charge,
            'IONMODE': _ion_mode,
            'PUBMED': f'{file_name}:scan:{ms2_scan}',
            'ACQUISITION': 'Crude',
            'EXACTMASS': row['exact_mass'],
            'DATACOLLECTOR': data_collector,
            'ADDUCT': row['t_adduct'],
            'CASNUMBER': str(row['CASNUMBER']) if 'CASNUMBER' in df.columns else None,
            'PI': pi_name,
            'LIBQUALITY': 1,
            'GENUS': None,
            'SPECIES': None,
            'INTEREST': None,
            'STRAIN': None
        })
        scan += 1

    return pd.DataFrame(rows)
