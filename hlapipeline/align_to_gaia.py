import os

from .utils import astrometric_utils

from drizzlepac import tweakreg

def align(expnames, **kwargs):
    """ Align exposures to a GAIA-based source catalog

    Parameters
    ===========
    expnames : str or list of strings
        Filename of exposure or list of filenames to be aligned to GAIA catalog

    """
    shift_name = kwargs.get('shift_name','shifts_gaia.txt')
    updatehdr = kwargs.get('updatehdr', False)
    wcsname = kwargs.get('wcsname', "TWEAK_GAIA")
    searchrad = kwargs.get('searchrad', 15.0)
    searchunits = kwargs.get('searchunits', 'arcseconds')
    ref_cat_file = kwargs.get('output', None)

    gaia_catalog = astrometric_utils.create_astrometric_catalog(expnames, **kwargs)
    gaia_catalog.write(ref_cat_file, format='ascii.no_header')

    if len(gaia_catalog) > 6:
        fitgeometry = 'general'
    else:
        fitgeometry = 'shift'

    tweakreg.TweakReg(expnames,
                      shiftfile=True,
                      interactive=False,
                      clean=True,
                      see2dplot=False,
                      updatehdr=updatehdr,
                      wcsname=wcsname,
                      outshifts=shift_name,
                      outwcs = shift_name.replace('.txt','_wcs.fits'),
                      refcat=ref_cat_file,
                      searchrad=searchrad,
                      searchunits=searchunits,
                      fitgeometry=fitgeometry)
    return gaia_catalog
