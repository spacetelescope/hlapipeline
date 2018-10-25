import os
import hlapipeline.utils.astrometric_utils as hlautils

from drizzlepac import tweakreg

def align(expnames, **kwargs):
    """ Align exposures to a GAIA-based source catalog

    Parameters
    ===========
    expnames : str or list of strings
        Filename of exposure or list of filenames to be aligned to GAIA catalog

    Returns
    =======
    gaia_catalog : Table
        Astropy Table object containing gaia catalog retrieved for exposures

    shift_name : string
        Filename of shift file written out by `tweakreg.Tweakreg`

    """
    setname = os.path.basename(expnames[0]).split('_')[0]
    shift_name = kwargs.get('shift_name',None)
    if shift_name is None:
        shift_name = '{}_shifts_gaia.txt'.format(setname)
    ref_cat_file = kwargs.get('ref_cat_file', '{}_gaia_ref.cat'.format(setname))
    # Set default values for specific Tweakreg parameters which are more
    # appropriate for most of our use cases
    updatehdr = kwargs.get('updatehdr', False)
    wcsname = kwargs.get('wcsname', "TWEAK_GAIA")
    searchrad = kwargs.get('searchrad', 15.0)
    searchunits = kwargs.get('searchunits', 'arcseconds')
    threshold = kwargs.get('threshold', 1000.0)
    conv_width = kwargs.get('conv_width', 3.5)

    gaia_catalog = hlautils.create_astrometric_catalog(expnames, **kwargs)
    gaia_catalog.write(ref_cat_file, format='ascii.no_header', overwrite=True)
    print("FINISHED writing out gaia_catalog: {}".format(os.path.abspath(ref_cat_file)))

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
                      threshold=threshold,
                      conv_width=conv_width,
                      outshifts=shift_name,
                      outwcs = shift_name.replace('.txt','_wcs.fits'),
                      refcat=ref_cat_file,
                      searchrad=searchrad,
                      searchunits=searchunits,
                      fitgeometry=fitgeometry)
    return gaia_catalog, shift_name
