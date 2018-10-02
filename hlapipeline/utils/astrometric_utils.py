"""Utilities to support creation of astrometrically accurate reference catalogs

The function, create_astrometric_catalog, allows the user to query an
astrometric catalog online to generate a catalog of astrometric sources that
should fall within the field-of-view of all the input images.

This module relies on the definition of an environment variable to specify
the URL of the astrometric catalog to use for generating this
reference catalog.

    ASTROMETRIC_CATALOG_URL  -- URL of web service that can be queried to
                                obtain listing of astrometric sources,
                                sky coordinates, and magnitudes.

"""
import os
import io
from io import BytesIO

import csv
import requests
from lxml import etree

import numpy as np

import stwcs
from stwcs.distortion import utils
from stsci.tools import fileutil as fu
from stsci.tools import parseinput

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as pf
from astropy.io import ascii
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import photutils
from photutils import detect_sources, source_properties
from photutils import Background2D, MedianBackground
from scipy.spatial import distance_matrix

import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

ASTROMETRIC_CAT_ENVVAR = "ASTROMETRIC_CATALOG_URL"
DEF_CAT_URL = 'http://gsss.stsci.edu/webservices'

if ASTROMETRIC_CAT_ENVVAR in os.environ:
    SERVICELOCATION = os.environ[ASTROMETRIC_CAT_ENVVAR]
else:
    SERVICELOCATION = DEF_CAT_URL


def create_astrometric_catalog(inputs, **pars):
    """Create an astrometric catalog that covers the inputs' field-of-view.

    Parameters
    ===========
    input : str
        Filenames of images to be aligned to astrometric catalog

    catalog : str, optional
        Name of catalog to extract astrometric positions for sources in the
        input images' field-of-view. Default: GSC241. Options available are
        documented on the catalog web page.

    output : str, optional
        Filename to give to the astrometric catalog read in from the master
        catalog web service.  Default: ref.cat

    gaia_only : bool, optional
        Specify whether or not to only use sources from GAIA in output catalog
        Default: False

    note ::
        This function will point to astrometric catalog web service defined
        through the use of the ASTROMETRIC_CATALOG_URL environment variable.
    """
    # interpret input parameters
    catalog = pars.get("catalog", 'GSC241')
    output = pars.get("output", 'ref.cat')
    gaia_only = pars.get("gaia_only", False)

    inputs, _ = parseinput.parseinput(inputs)
    # start by creating a composite field-of-view for all inputs
    wcslist = []
    for img in inputs:
        nsci = fu.countExtn(img)
        for num in range(nsci):
            extname = '{}[sci,{}]'.format(img, num+1)
            wcslist.append(stwcs.wcsutil.HSTWCS(extname))

    # This default output WCS will have the same plate-scale and orientation
    # as the first chip in the list, which for WFPC2 data means the PC.
    # Fortunately, for alignment, this doesn't matter since no resampling of
    # data will be performed
    outwcs = utils.output_wcs(wcslist)
    radius = compute_radius(outwcs)
    ra, dec = outwcs.wcs.crval

    # perform query for this field-of-view
    ref_dict = get_catalog(ra, dec, sr=radius, catalog=catalog)
    row = "{:14.8f}  {:14.8f}  {:8.4f} {} {}\n"
    num_sources = 0
    with open(output, 'w') as refcat:
        refcat.write("#ra  dec   mag   objID   GaiaID\n")
        for source in ref_dict:
            if 'GAIAsourceID' in source:
                g = source['GAIAsourceID']
                if gaia_only and g.strip() is '':
                    continue
            else:
                g = -1  # indicator for no source ID extracted
            r = float(source['ra'])
            d = float(source['dec'])
            m = float(source['mag'])
            o = source['objID']
            refcat.write(row.format(r, d, m, o, g))
            num_sources += 1

    print("Created catalog '{}' with {} sources".format(output, num_sources))

def get_catalog(ra, dec, sr=0.1, fmt='CSV', catalog='GSC241'):
    """ Extract catalog from VO web service.

    Parameters
    ----------
    ra : float
        Right Ascension (RA) of center of field-of-view (in decimal degrees)

    dec : float
        Declination (Dec) of center of field-of-view (in decimal degrees)

    sr : float, optional
        Search radius (in decimal degrees) from field-of-view center to use
        for sources from catalog.  Default: 0.1 degrees

    fmt : str, optional
        Format of output catalog to be returned.  Options are determined by
        web-service, and currently include (Default: CSV):
            VOTABLE(default) | HTML | KML | CSV | TSV | JSON | TEXT

    catalog : str, optional
        Name of catalog to query, as defined by web-service.  Default: 'GSC241'

    Returns
    -------
    csv : obj
        CSV object of returned sources with all columns as provided by catalog

    """
    serviceType = 'vo/CatalogSearch.aspx'
    spec_str = 'RA={}&DEC={}&SR={}&FORMAT={}&CAT={}&MINDET=5'
    headers = {'Content-Type': 'text/csv'}

    spec = spec_str.format(ra, dec, sr, fmt, catalog)
    serviceUrl = '{}/{}?{}'.format(SERVICELOCATION, serviceType,spec)
    rawcat = requests.get(serviceUrl, headers=headers)
    r_contents = rawcat.content.decode()  # convert from bytes to a String
    rstr = r_contents.split('\r\n')
    # remove initial line describing the number of sources returned
    # CRITICAL to proper interpretation of CSV data
    del rstr[0]
    r_csv = csv.DictReader(rstr)

    return r_csv


def compute_radius(wcs):
    """Compute the radius from the center to the furthest edge of the WCS."""

    ra,dec = wcs.wcs.crval
    img_center = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    wcs_foot = wcs.calc_footprint()
    img_corners = SkyCoord(ra=wcs_foot[:,0]*u.degree,
                           dec=wcs_foot[:,1]*u.degree)
    radius = img_center.separation(img_corners).max().value

    return radius

def find_gsc_offset(image, input_catalog='GSC1', output_catalog='GAIA'):
    """Find the GSC to GAIA offset based on guide star coordinates

    Parameters
    ----------
    image : str
        filename of image to be processed

    Returns
    -------
    delta_ra,delta_dec : tuple of floats
        Offset in decimal degrees of image based on correction to guide star
        coordinates relative to GAIA
    """
    serviceType = "GSCConvert/GSCconvert.aspx"
    spec_str = "TRANSFORM={}-{}&IPPPSSOOT={}"

    if 'rootname' in pf.getheader(image):
        ippssoot = pf.getval(image, 'rootname').upper()
    else:
        ippssoot = fu.buildNewRootname(image).upper()

    spec = spec_str.format(input_catalog, output_catalog, ippssoot)
    serviceUrl = "{}/{}?{}".format(SERVICELOCATION, serviceType,spec)
    rawcat = requests.get(serviceUrl)
    if not rawcat.ok:
        print("Problem accessing service with:\n{{}".format(serviceUrl))
        raise ValueError

    delta_ra = delta_dec = None
    tree = BytesIO(rawcat.content)
    for _,element in etree.iterparse(tree):
        if element.tag == 'deltaRA':
            delta_ra = float(element.text)
        elif element.tag == 'deltaDEC':
            delta_dec = float(element.text)

    return delta_ra,delta_dec


def extract_sources(img, fwhm=2.0, threshold=None, source_box=7,
                    sharp=None, output=None, plot=False, vmax=None):
    """Use photutils to find sources in image based on segmentation."""
    if threshold is None:
        bkg_estimator = MedianBackground()
        bkg = Background2D(img, (50, 50), filter_size=(3, 3),
                           bkg_estimator=bkg_estimator)
        threshold = bkg.background + (3. * bkg.background_rms)
    sigma = fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=source_box, y_size=source_box)
    kernel.normalize()
    segm = detect_sources(img, threshold, npixels=source_box,
                          filter_kernel=kernel)
    cat = source_properties(img, segm)
    print("Total Number of detected sources: {}".format(len(cat)))
    if sharp:
        # Remove sources that do not fall within specified sharpness limit
        newcat = photutils.segmentation.properties.SourceCatalog([])
        for obj in cat:
            src_peak = (obj.max_value - (obj.source_sum/obj.area.value))
            kern_peak = (kernel*obj.source_sum).array.max()
            sharpness = src_peak/kern_peak
            if sharpness < sharp[0] or sharpness > sharp[1]:
                newcat._data.append(obj)
    else:
        newcat = cat
    print("Final Number of selected sources: {}".format(len(newcat)))
    if output:
        tbl = newcat.to_table()
        tbl['xcentroid'].info.format = '.10f'  # optional format
        tbl['ycentroid'].info.format = '.10f'
        tbl['cxy'].info.format = '.10f'
        tbl['cyy'].info.format = '.10f'
        tbl.write(output, format='ascii.ecsv')

    if plot:
        norm = None
        if vmax is None:
            norm = ImageNormalize(stretch=SqrtStretch())
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
        ax1.imshow(img, origin='lower', cmap='Greys_r', norm=norm, vmax=vmax)
        ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))

    return newcat, segm


def find_isolated_source(catalog, columns=None):
    """Find the source in the catalog which is most isolated from all others.

    PARAMETERS
    -----------
    catalog : obj
        list of source positions with (at a minimum) ra,dec or x,y positions
        Catalog could come from photutils.source_properties, for example.

    columns : list, optional
        List of column names for source positions as provided by the catalog.
        If specified (not None), convert input catalog using these column names.
        Default: None.

    Returns
    --------
    index : int
        Index of source from catalog list for source furthest from all others

    nearest_dist : float
        Distance from isolated source to it's nearest neighbor in the catalog

    .. note ::
        This function assumes that the catalog has already been limited to
        only sources which overlap the actual image area.

    """
    if columns:
        cat1 = np.column_stack((catalog[columns[0]], catalog[columns[1]]))
    else:
        cat1 = catalog
    cat_dist = distance_matrix(cat1, cat1)
    y_max = cat_dist.sum(axis=1)
    iso_indx = y_max.argsort()[-1]

    # get index to next nearest neighbor to this isolated source
    next_nearest = cat_dist[iso_indx].argsort()[1]
    # now remember the distance to that source from the isolated source
    iso_dist = cat_dist[iso_indx][next_nearest]

    return iso_indx, iso_dist


def within_footprint(img, wcs, x, y):
    """Determine whether input x,y fall in the science area of the image.

    Parameters
    -----------
    wcs : obj
        HSTWCS or WCS object with naxis terms defined

    img : ndarray
        ndarray of image where non-science areas are marked with value of NaN

    x,y : arrays
        arrays of x,y positions for sources to be checked

    Returns
    -------
    x,y : arrays
        New arrays which have been trimmed of all sources that fall outside
        the science areas of the image

    """
    # start with limits of WCS shape
    maskx = np.bitwise_or(x<0, x>wcs.naxis1)
    masky = np.bitwise_or(y<0, y>wcs.naxis2)
    mask = ~np.bitwise_or(maskx,masky)
    x = x[mask]
    y = y[mask]

    # Now, confirm that these points fall within actual science area of WCS
    nanmask = np.isnan(img[x.astype(np.int32),y.astype(np.int32)])
    x = x[~nanmask]
    y = y[~nanmask]
    return x,y


def find_best_offset(image, wcs, reference, refnames=['ra', 'dec'],
                     match_tolerance=5.):
    """Iteratively look for the best cross-match between the catalog and ref.

    Parameters
    ----------
        image : ndarray
            Image (as ndarray) of image to extract sources for matching to
            the external astrometric catalog

        wcs : object
            WCS object which provides translation from sky coordinates to
            image coordinates for input image

        reference : str or object
            Reference catalog, either as a filename or ``astropy.Table``
            containing astrometrically accurate sky coordinates for astrometric
            standard sources

        refnames : list
            List of table column names for sky coordinates of astrometric
            standard sources from reference catalog

        match_tolerance : float
            Tolerance (in pixels) for recognizing that a source position matches
            an astrometric catalog position.  Larger values allow for lower
            accuracy source positions to be compared to astrometric catalog
            Default: 5 pixels
    Returns
    -------
        best_offset : tuple
            Offset in input image pixels between image source positions and
            astrometric catalog positions that results in largest number of
            matches of astrometric sources with image sources
    """
    # read in reference catalog
    if isinstance(reference, str):
        refcat = ascii.read(reference)
    else:
        refcat = reference
    ref_ra = refcat[refnames[0]]
    ref_dec = refcat[refnames[1]]
    reftab = np.column_stack((ref_ra, ref_dec))
    xref, yref = wcs.all_world2pix(ref_ra, ref_dec, 1)

    # look for the most isolated reference source to serve as a
    # zero-point/anchor
    xref, yref = within_footprint(image, wcs, xref, yref)
    xyarr = np.column_stack((xref, yref))
    # assumption: catalog is column_stack((xc,yc))
    ref_iso, iso_dist = find_isolated_source(xyarr)
    print("Found isolated source at position : {},{}".format(xyarr[ref_iso]))
    print("  with separation from neighbor of: {}".format(iso_dist))

    # compute match limit based on distance to neighbor
    #  WARNING ::
    #  hard-coded limit: 55% from isolated source to neighbor
    #  This allows for a little overlap in areas searched for matches
    match_limit = iso_dist*0.55

    # find sources in image
    seg_cat, segmap = extract_sources(image)
    seg_tab = seg_cat.to_table()
    seg_xy = np.column_stack((seg_tab['xcentroid'], seg_tab['ycentroid']))
    seg_xy = seg_xy[~np.isnan(seg_xy[:, 0])]

    # look for nearest matches to isolated reference source
    distxy = distance_matrix(xyarr[ref_iso], seg_xy)
    maskdist = (distxy < match_limit)[0]
    # identify x,y positions of sources nearest isolated reference
    # close_matches = seg_xy[maskdist]

    # Now, start iterating through all combinations
    delta_refs = seg_xy[maskdist].value - xyarr[ref_iso]
    num_matches = []
    for delta in delta_refs:
        mdist = distance_matrix(xyarr+delta, seg_xy)
        num_matches.append(len(np.where(mdist < match_tolerance)[0]))
    max_matches = max(num_matches)
    best_offset = delta_refs[num_matches.index(max_matches)]
    if max_matches < int(len(xyarr)*0.1):
        best_offset = None
    return best_offset
