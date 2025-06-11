from astropy.table import Table

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

def read_spectrum(filename):
    """
    Reads a FITS file from SDSS and returns the spectrum data as an Astropy Table.
    
    Parameters:
    filename (str): The path to the FITS file.
    
    Returns:
    astropy.table.Table: The spectrum data.
    """
    
    t = Table.read(filename, hdu='COADD')

    wavelength = 10**t['loglam']
    flux = t['flux']

    return wavelength, flux
