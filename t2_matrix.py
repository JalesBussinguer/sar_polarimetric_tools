import rasterio as rio
import numpy as np

def t2_matrix(svv_re, svv_im, svh_re, svh_im):

    """Calculate the coherency matrix components

    Args:
        svv_re: VV band (real)
        svv_im: VV band (imag)
        svh_re: VH band (real)
        svh_im: VH band (imag)

    Returns:
        T11, T12_real, T12_imag, T22 [Array]: _description_
    """

    # k1 = Svv
    # k2 = 2Svh

    k1_re = svv_re
    k1_im = svv_im
    k2_re = 2 * svh_re
    k2_im = 2 * svh_im

    t11 = k1_re * k1_re + k1_im * k1_im

    t12_real = k1_re * k2_re + k1_im * k2_im
    t12_imag = k1_im * k2_re - k1_re * k2_im

    t22 = k2_re * k2_re + k2_im * k2_im

    return t11, t12_real, t12_imag, t22

if __name__ == "__main__":


    with rio.open('path_to_SLC_complex_data') as src:
        # Read the C2 matrix components as numpy arrays
        svv_re = src.read(1)
        svv_im = src.read(2)
        svh_re = src.read(3)
        svh_im = src.read(4)

    t11, t12_real, t12_imag, t22 = t2_matrix(svv_re, svv_im, svh_re, svh_im)

    bands_dict = {'t11': t11, 't12_real': t12_real, 't2_imag': t12_imag, 't22': t22}

    profile = src.meta.copy()

    with rio.open('outpath/T2.tif', 'w', **profile) as dest:
        for id, band_name in enumerate(bands_dict.keys(), start=1):
            dest.write(bands_dict[band_name], id)
