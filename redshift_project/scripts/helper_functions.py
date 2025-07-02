from astropy.table import Table
import numpy as np
    

pretty_plot_params = {'font.size': 18, 
                     'text.usetex': True, 
                     'font.family': 'STIXGeneral',
                     'xtick.top': True,
                     'ytick.right': True,
                     'xtick.major.size': 10.,
                     'xtick.minor.size': 5.,
                     'xtick.major.width': 1.5,
                     'xtick.minor.width': 1.,
                     'ytick.major.size': 10.,
                     'ytick.minor.size': 5.,
                     'ytick.major.width': 1.5,
                     'ytick.minor.width': 1.,
                     'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'xtick.minor.visible': True,
                     'ytick.minor.visible': True}



def get_continuum_mask(wave, wave_range, buffer=20):
    """
    Function to get a boolean mask for the continuum regions based on the given wavelength range.
    
    Parameters:
    wave (array-like): Array of wavelength values.
    wave_range (list of tuples): List of wavelength ranges to exclude from the continuum.
    buffer (float): Buffer region around each wavelength range to exclude.
    
    Returns:
    np.ndarray: Boolean mask where True indicates continuum regions.
    """

    mask = (wave >= wave_range[0]) & (wave < wave_range[1])
    mask = ~mask & ((wave >= wave_range[0] - buffer) & 
                    (wave < wave_range[1] + buffer))
    return mask

def subtract_continuum(wave, flux, wave_range=None, buffer=20):
    """
    Function to subtract the continuum from a spectrum.
    
    Parameters:
    wave (array-like): Array of wavelength values.
    flux (array-like): Array of flux values.
    wave_range (list of tuples): List of wavelength ranges to exclude from the continuum.
    buffer (float): Buffer region around each wavelength range to exclude.
    
    Returns:
    tuple: Continuum-subtracted flux and the fitted continuum.
    """

    if wave_range is None:
        wave_range = [4980, 5020]
    
    mask = np.ones_like(wave, dtype=bool)
    for wr in wave_range:
        mask &= get_continuum_mask(wave, wr, buffer=buffer)
    
    # Fit a 3rd degree polynomial to the continuum regions
    p_coeff = np.polyfit(wave[mask], flux[mask], deg=3)
    continuum = np.polyval(p_coeff, wave)
    
    # Subtract the fitted continuum from the original flux
    flux_cont_subtracted = flux - continuum
    
    return flux_cont_subtracted, continuum