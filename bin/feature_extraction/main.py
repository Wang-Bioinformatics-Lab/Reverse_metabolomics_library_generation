import os
from .config import Params, find_ms_info
from .raw_data_utils import MSData


def feature_extraction_single(file_path, mass_detect_int_tol=None,
                              min_feature_height=1.5e5,
                              peak_cor_rt_tol=0.025, min_ppc=0.8,
                              ms1_tol=0.005, ms2_tol=0.015,
                              save=True, out_dir=None):
    """
    feature detection from a single file
    """

    ms_type, ion_mode, centroided = find_ms_info(file_path)

    # init a new config object
    config = init_config(ms_type, ion_mode=ion_mode,
                         mz_tol_ms1=ms1_tol, mz_tol_ms2=ms2_tol,
                         peak_cor_rt_tol=peak_cor_rt_tol, min_ppc=min_ppc,
                         mass_detect_int_tol=mass_detect_int_tol)

    # create a MSData object
    d = MSData()

    # read raw data
    d.read_raw_data(file_path, config)

    # detect region of interests (ROIs)
    d.find_rois()

    # cut ROIs
    d.cut_rois()

    # label short ROIs, find the best MS2, and sort ROIs by m/z
    d.summarize_roi()

    # output single file to a tsv file, in the same directory as the raw file
    print('Generating feature table...')
    if out_dir is None:
        out_dir = os.path.dirname(file_path)

    df = d.output_single_file(save, out_dir)

    # filter by ROI length
    df = df[df['length'] >= 4].reset_index(drop=True)

    # filter by intensity
    df = df[df['peak_height'] >= min_feature_height].reset_index(drop=True)

    return df, ion_mode


def init_config(ms_type, ion_mode="positive",
                mz_tol_ms1=0.005, mz_tol_ms2=0.015,
                peak_cor_rt_tol=0.025, min_ppc=0.8,
                mass_detect_int_tol=5e4):
    # init
    config = Params()

    if mass_detect_int_tol is None:
        mass_detect_int_tol = 5e4 if ms_type == "orbitrap" else 1e3

    ##########################
    # MS data acquisition
    config.ion_mode = ion_mode  # Ionization mode, "positive" or "negative"
    config.mz_tol_ms1 = mz_tol_ms1  # m/z tolerance for MS1
    config.mz_tol_ms2 = mz_tol_ms2  # m/z tolerance for MS2
    config.int_tol = mass_detect_int_tol  # Mass detection intensity tolerance
    config.ppr = min_ppc  # Peak-peak correlation threshold for feature grouping
    config.peak_cor_rt_tol = peak_cor_rt_tol  # RT tolerance for peak correlation

    return config
