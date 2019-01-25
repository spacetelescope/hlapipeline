#!/usr/bin/env python

"""This script is a modernized replacement of tweakreg.

"""

import datetime
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
from astropy import units as u
from collections import OrderedDict
from drizzlepac import updatehdr
import glob
import math
import numpy as np
import os
import pdb
from stsci.tools import fileutil
from stwcs.wcsutil import HSTWCS
import sys
import tweakwcs
try:
    from hlapipeline.utils import astrometric_utils as amutils
    from hlapipeline.utils import astroquery_utils as aqutils
    from hlapipeline.utils import filter
    from hlapipeline.utils import get_git_rev_info
except:
    from utils import astrometric_utils as amutils
    from utils import astroquery_utils as aqutils
    from utils import filter
    from utils import get_git_rev_info

MIN_CATALOG_THRESHOLD = 3
MIN_OBSERVABLE_THRESHOLD = 10
MIN_CROSS_MATCHES = 3
MIN_FIT_MATCHES = 6
MAX_FIT_RMS = 10 # RMS now in mas, 1.0
MAX_FIT_LIMIT = 1000 # Maximum RMS that a result is useful
MAX_SOURCES_PER_CHIP = 250  # Maximum number of sources per chip to include in source catalog

# Module-level dictionary contains instrument/detector-specific parameters used later on in the script.
detector_specific_params = {"acs":
                                {"hrc":
                                     {"fwhmpsf": 0.073,
                                      "classify": True,
                                      "threshold": None},
                                 "sbc":
                                     {"fwhmpsf": 0.065,
                                      "classify": False,
                                      "threshold": 2.0},
                                 "wfc":
                                     {"fwhmpsf": 0.13, #0.076,
                                      "classify": True,
                                      "threshold": -1.1}},
                            "wfc3":
                                {"ir":
                                     {"fwhmpsf": 0.14,
                                      "classify": False,
                                      "threshold": None},
                                 "uvis":
                                     {"fwhmpsf": 0.076,
                                      "classify": True,
                                      "threshold": None}}} # fwhmpsf in units of arcsec


# ----------------------------------------------------------------------------------------------------------------------


def check_and_get_data(input_list,**pars):
    """Verify that all specified files are present. If not, retrieve them from MAST.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    input_file_list : list
        list of full filenames

    """
    totalInputList=[]
    for input_item in input_list:
        filelist = aqutils.retrieve_observation(input_item,**pars)
        if len(filelist) == 0:
            # look for local copy of the file
            fitsfilenames = sorted(glob.glob("{}_fl?.fits".format(input_item)))
            if len(fitsfilenames) > 0:
                imghdu = fits.open(fitsfilenames[0])
                imgprimaryheader = imghdu[0].header
                try:
                    asnid = imgprimaryheader['ASN_ID'].strip().lower()
                    if asnid == 'NONE':
                        asnid = None
                except KeyError:
                    asnid = None
                if asnid:
                    filelist = aqutils.retrieve_observation(asnid,**pars)
        if len(filelist) > 0:
            totalInputList += filelist

    if len(filelist) > 0: totalInputList = sorted(
        list(set(totalInputList)))  # remove duplicate list elements, sort resulting list of unique elements
    print("TOTAL INPUT LIST: ",totalInputList)
    # TODO: add trap to deal with non-existent (incorrect) rootnames
    # TODO: Address issue about how the code will retrieve association information if there isn't a local file to get 'ASN_ID' header info
    return(totalInputList)


# ----------------------------------------------------------------------------------------------------------------------

def convert_string_tf_to_boolean(invalue):
    """Converts string 'True' or 'False' value to Boolean True or Boolean False.

    :param invalue: string
        input true/false value

    :return: Boolean
        converted True/False value
    """
    outvalue = False
    if invalue == 'True':
        outvalue = True
    return(outvalue)


# ----------------------------------------------------------------------------------------------------------------------


def perform_align(input_list, archive=False, clobber=False, update_hdr_wcs=False, print_fit_parameters=True,
                    print_git_info=False):
    """Main calling function.

    Parameters
    ----------
    input_list : list
        List of one or more IPPSSOOTs (rootnames) to align.

    archive : Boolean
        Retain copies of the downloaded files in the astroquery created sub-directories?

    clobber : Boolean
        Download and overwrite existing local copies of input files?

    update_hdr_wcs : Boolean
        Write newly computed WCS information to image image headers?

    print_fit_parameters : Boolean
        Specify whether or not to print out FIT results for each chip.

    print_git_info : Boolean
        Display git repository information?

    Returns
    -------
    int value 0 if successful, int value 1 if unsuccessful

    """

    # Define astrometric catalog list in priority order
    catalogList = ['GAIADR2', 'GSC241']

    # 0: print git info
    if print_git_info:
        print("-------------------- STEP 0: Display Git revision info  --------------------")
        full_path = os.path.dirname(__file__)+"/utils"
        repo_path=None
        if "hlapipeline/hlapipeline" in full_path:
            repo_path = full_path.split("hlapipeline/hlapipeline")[0]+"hlapipeline"
        elif "hlapipeline" in full_path:
            repo_path = full_path.split("hlapipeline")[0]+"hlapipeline"
        else:
            pass
        if not os.path.exists(repo_path): repo_path = None # protect against non-existent paths
        if repo_path:
            get_git_rev_info.print_rev_id(repo_path) # Display git repository information
        else:
            print("WARNING: Unable to display Git repository revision information.")

    # 1: Interpret input data and optional parameters
    print("-------------------- STEP 1: Get data --------------------")
    zeroDT = startingDT = datetime.datetime.now()
    print(str(startingDT))
    imglist = check_and_get_data(input_list, archive=archive, clobber=clobber)
    print("\nSUCCESS")

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 1]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 2: Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
    print("-------------------- STEP 2: Filter data --------------------")
    filteredTable = filter.analyze_data(imglist)

    # Check the table to determine if there is any viable data to be aligned.  The
    # 'doProcess' column (bool) indicates the image/file should or should not be used
    # for alignment purposes.  For filtered data, 'doProcess=0' and 'status=9999' in the table
    # (the status value by default), so there is no need to update the filteredTable here.
    if filteredTable['doProcess'].sum() == 0:
        print("No viable images in filtered table - no processing done.\n")
        return(filteredTable)

    # Get the list of all "good" files to use for the alignment
    processList = filteredTable['imageName'][np.where(filteredTable['doProcess'])]
    processList = list(processList) #Convert processList from numpy list to regular python list
    print("\nSUCCESS")

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 2]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 3: Build WCS for full set of input observations
    print("-------------------- STEP 3: Build WCS --------------------")
    refwcs = amutils.build_reference_wcs(processList)
    print("\nSUCCESS")


    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 3]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 4: Extract catalog of observable sources from each input image
    print("-------------------- STEP 4: Source finding --------------------")
    extracted_sources = generate_source_catalogs(processList,
                                                 centering_mode='starfind',
                                                 nlargest=MAX_SOURCES_PER_CHIP)

    for imgname in extracted_sources.keys():
        table=extracted_sources[imgname]["catalog_table"]
        # The catalog of observable sources must have at least MIN_OBSERVABLE_THRESHOLD entries to be useful
        total_num_sources = 0
        for chipnum in table.keys():
            total_num_sources += len(table[chipnum])

        # Update filtered table with number of found sources
        index = np.where(filteredTable['imageName']==imgname)[0][0]
        filteredTable[index]['foundSources'] = total_num_sources

        if total_num_sources < MIN_OBSERVABLE_THRESHOLD:
            print("Not enough sources ({}) found in image {}".format(total_num_sources,imgname))
            filteredTable[index]['status'] = 1
            return(filteredTable)
    print("\nSUCCESS")
    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 4]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 5: Retrieve list of astrometric sources from database

    # Convert input images to tweakwcs-compatible NDData objects and
    # attach source catalogs to them.
    imglist = []
    for group_id, image in enumerate(processList):
        img = amutils.build_nddata(image, group_id,
                                   extracted_sources[image]['catalog_table'])
        # add the name of the image to the imglist object
        for im in img:
            im.meta['name'] = image
        imglist.extend(img)

    best_fit_rms = -99999.0
    fit_algorithm_list= [match_2dhist_fit,match_default_fit]
    for catalogIndex in range(0, len(catalogList)): #loop over astrometric catalog
        print("-------------------- STEP 5: Detect astrometric sources --------------------")
        print("Astrometric Catalog: ",catalogList[catalogIndex])
        reference_catalog = generate_astrometric_catalog(processList, catalog=catalogList[catalogIndex])

        currentDT = datetime.datetime.now()
        deltaDT = (currentDT - startingDT).total_seconds()
        print('Processing time of [STEP 5]: {} sec'.format(deltaDT))
        startingDT = currentDT

        if len(reference_catalog) < MIN_CATALOG_THRESHOLD:
            print("Not enough sources found in catalog " + catalogList[catalogIndex])
            if catalogIndex < len(catalogList) -1:
                print("Try again with other catalog")
            else:
                print("ERROR! No astrometric sources found in any catalog. Exiting...") #bail out if not enough sources can be found any of the astrometric catalogs
                filteredTable['status'][:] = 1
                return (filteredTable)
        else:
            print("-------------------- STEP 5b: Cross matching and fitting --------------------")
            for algorithm_name in fit_algorithm_list: #loop over fit algorithm type
                print("------------------ Catalog {} matched using {} ------------------ ".format(catalogList[catalogIndex],algorithm_name.__name__))

                #execute the correct fitting/matching algorithm
                fit_rms, fit_num = algorithm_name(imglist, reference_catalog, print_fit_parameters=print_fit_parameters)

                # update the best fit
                if best_fit_rms >= 0.:
                    if fit_rms < best_fit_rms:
                        best_fit_rms = fit_rms
                        best_fit_num = fit_num
                        for item in imglist:
                            item.best_meta = item.meta.copy()
                else:
                    if fit_rms < MAX_FIT_LIMIT:
                        best_fit_rms = fit_rms
                        best_fit_num = fit_num
                        for item in imglist:
                            item.best_meta = item.meta.copy()
                #imglist_temp = imglist.copy() # preserve best fit solution so that it can be inserted into a reinitialized imglist next time through.


    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 5b]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 6: Populate the filteredTable
    print("-------------------- STEP 6: Collect up information and populate the filtered table --------------------")
    if best_fit_rms > 0 and best_fit_rms < MAX_FIT_RMS:
        print("The fitting process was successful with a best fit total rms of {} mas".format(best_fit_rms))
    else:
        print("The fitting process was unsuccessful with a best fit total rms of {} mas".format(best_fit_rms))

    if best_fit_rms > 0 and best_fit_rms < MAX_FIT_LIMIT:
        # update to the meta information with the lowest rms if it is reasonable
        for item in imglist:
            item.meta = item.best_meta.copy()
        filteredTable['status'][:] = 0

    info_keys = OrderedDict(imglist[0].meta['tweakwcs_info']).keys()
    # Update filtered table with number of matched sources and other information
    for item in imglist:
        imgname = item.meta['name']
        index = np.where(filteredTable['imageName'] == imgname)[0][0]

        if item.meta['tweakwcs_info']['status'].startswith("FAILED") != True:
            for tweakwcs_info_key in info_keys:
                if not tweakwcs_info_key.startswith("matched"):
                    if tweakwcs_info_key.lower() == 'rms':
                        filteredTable[index]['rms_x'] = item.meta['tweakwcs_info'][tweakwcs_info_key][0]
                        filteredTable[index]['rms_y'] = item.meta['tweakwcs_info'][tweakwcs_info_key][1]

            filteredTable[index]['catalog'] = item.meta['tweakwcs_info']['catalog']
            filteredTable[index]['catalogSources'] = len(reference_catalog)
            filteredTable[index]['matchSources'] = item.meta['tweakwcs_info']['nmatches']
            filteredTable[index]['rms_ra'] = item.meta['tweakwcs_info']['RMS_RA'].value
            filteredTable[index]['rms_dec'] = item.meta['tweakwcs_info']['RMS_DEC'].value
            filteredTable[index]['fit_rms'] = item.meta['tweakwcs_info']['FIT_RMS']
            filteredTable[index]['total_rms'] = item.meta['tweakwcs_info']['TOTAL_RMS']
            #filteredTable.pprint(max_width=-1)

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 6]: {} sec'.format(deltaDT))
    startingDT = currentDT
    # 7: Write new fit solution to input image headers
    print("-------------------- STEP 7: Update image headers with new WCS information --------------------")
    if best_fit_rms > 0 and update_hdr_wcs:
        update_image_wcs_info(imglist, processList)
        print("\nSUCCESS")
    else:
        print("\n STEP SKIPPED")

    currentDT = datetime.datetime.now()
    deltaDT = (currentDT - startingDT).total_seconds()
    print('Processing time of [STEP 7]: {} sec'.format(deltaDT))
    print('TOTAL Processing time of {} sec'.format((currentDT- zeroDT).total_seconds()))
    return (filteredTable)

def match_default_fit(imglist, reference_catalog, print_fit_parameters=True):
    """Perform cross-matching and final fit using 2dHistogram matching

    Parameters
    ----------
    imglist : list
        List of input image NDData objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    print_fit_parameters : bool
        Specify whether or not to print out FIT results for each chip


    Returns
    --------
    fit_rms : float
        Visit level RMS for the FIT

    fit_num : int
        Number of sources used to generate visit level FIT and `fit_rms`

    """
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=250, separation=0.1,
                             tolerance=100, use2dhist=False)
    # Align images and correct WCS
    tweakwcs.tweak_image_wcs(imglist, reference_catalog, match=match)
    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    # determine the quality of the fit
    fit_rms, fit_num  = determine_fit_quality(imglist, print_fit_parameters=print_fit_parameters)

    return fit_rms, fit_num


def match_2dhist_fit(imglist, reference_catalog, print_fit_parameters=True):
    """Perform cross-matching and final fit using 2dHistogram matching

    Parameters
    ----------
    imglist : list
        List of input image NDData objects with metadata and source catalogs

    reference_catalog : Table
        Astropy Table of reference sources for this field

    print_fit_parameters : bool
        Specify whether or not to print out FIT results for each chip


    Returns
    --------
    fit_rms : float
        Visit level RMS for the FIT

    fit_num : int
        Number of sources used to generate visit level FIT and `fit_rms`

    """
    print("-------------------- STEP 5b: (match_2dhist_fit) Cross matching and fitting --------------------")
    # Specify matching algorithm to use
    match = tweakwcs.TPMatch(searchrad=75, separation=0.1,
                             tolerance=2.0, use2dhist=True)
    # Align images and correct WCS
    tweakwcs.tweak_image_wcs(imglist, reference_catalog, match=match)
    # Interpret RMS values from tweakwcs
    interpret_fit_rms(imglist, reference_catalog)

    # determine the quality of the fit
    fit_rms, fit_num  = determine_fit_quality(imglist, print_fit_parameters=print_fit_parameters)

    return fit_rms, fit_num

def determine_fit_quality(imglist, print_fit_parameters=True):
    """Determine the quality of the fit to the data

    Parameters
    ----------
    imglist : list
        output of interpret_fits. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
        input image.  This list should have been  updated, in-place, with the new RMS values;
        specifically,

            * 'FIT_RMS': RMS of the separations between fitted image positions and reference positions
            * 'TOTAL_RMS': mean of the FIT_RMS values for all observations
            * 'NUM_FITS': number of images/group_id's with successful fits included in the TOTAL_RMS

        These entries are added to the 'tweakwcs_info' dictionary.

    print_fit_parameters : bool
        Specify whether or not to print out FIT results for each chip

    Returns
    -------
    max_rms_val : float
        The best Totol rms dteremined from all of the images

    num_xmatches: int
        The number of stars used in matching the data

    """
    tweakwcs_info_keys = OrderedDict(imglist[0].meta['tweakwcs_info']).keys()
    max_rms_val = 1e9
    num_xmatches = 0

    for item in imglist:
        image_name = item.meta['name']
        #Handle fitting failures (no matches found)
        if item.meta['tweakwcs_info']['status'].startswith("FAILED") == True:
                print("No cross matches found in any catalog for {} - no processing done.".format(image_name))
                continue
        fit_rms_val = item.meta['tweakwcs_info']['FIT_RMS']
        max_rms_val = item.meta['tweakwcs_info']['TOTAL_RMS']
        num_xmatches = item.meta['tweakwcs_info']['nmatches']
        if num_xmatches < MIN_CROSS_MATCHES:
            if catalogIndex < numCatalogs-1:
                print("Not enough cross matches found between astrometric catalog and sources found in {}".format(image_name))
                continue
        print('RESULTS FOR {} Chip {}: FIT_RMS = {} mas, TOTAL_RMS = {} mas, NUM =  {}'.format(image_name, item.meta['chip'], fit_rms_val, max_rms_val, num_xmatches))
        # print fit params to screen
        if print_fit_parameters:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIT PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("image: {}".format(image_name))
            print("chip: {}".format(item.meta['chip']))
            print("group_id: {}".format(item.meta['group_id']))
            for tweakwcs_info_key in tweakwcs_info_keys:
                if not tweakwcs_info_key.startswith("matched"):
                    print("{} : {}".format(tweakwcs_info_key,item.meta['tweakwcs_info'][tweakwcs_info_key]))
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    if max_rms_val > MAX_FIT_RMS:
        print("Total fit RMS value = {} mas greater than the maximum threshold value {}.".format(max_rms_val, MAX_FIT_RMS))
        print("Try again with the next catalog")
    else:
        print("Fit calculations successful.")

    return max_rms_val, num_xmatches
# ----------------------------------------------------------------------------------------------------------------------


def generate_astrometric_catalog(imglist, **pars):
    """Generates a catalog of all sources from an existing astrometric catalog are in or near the FOVs of the images in
        the input list.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for catalog generation.

    Returns
    =======
    ref_table : object
        Astropy Table object of the catalog

    """
    # generate catalog
    out_catalog = amutils.create_astrometric_catalog(imglist,**pars)

    # if the catalog has contents, write the catalog to ascii text file
    if len(out_catalog) > 0:
        catalog_filename = "refcatalog.cat"
        out_catalog.write(catalog_filename, format="ascii.fast_commented_header")
        print("Wrote reference catalog {}.".format(catalog_filename))

    return(out_catalog)


# ----------------------------------------------------------------------------------------------------------------------


def generate_source_catalogs(imglist, **pars):
    """Generates a dictionary of source catalogs keyed by image name.

    Parameters
    ----------
    imglist : list
        List of one or more calibrated fits images that will be used for source detection.

    Returns
    -------
    sourcecatalogdict : dictionary
        a dictionary (keyed by image name) of two element dictionaries which in tern contain 1) a dictionary of the
        detector-specific processing parameters and 2) an astropy table of position and photometry information of all
        detected sources
    """
    output = pars.get('output', False)
    sourcecatalogdict = {}
    for imgname in imglist:
        print("Image name: ", imgname)

        sourcecatalogdict[imgname] = {}

        # open image
        imghdu = fits.open(imgname)
        imgprimaryheader = imghdu[0].header
        instrument = imgprimaryheader['INSTRUME'].lower()
        detector = imgprimaryheader['DETECTOR'].lower()

        # get instrument/detector-specific image alignment parameters
        if instrument in detector_specific_params.keys():
            if detector in detector_specific_params[instrument].keys():
                detector_pars = detector_specific_params[instrument][detector]
                # to allow generate_source_catalog to get detector specific parameters
                detector_pars.update(pars)
                sourcecatalogdict[imgname]["params"] = detector_pars
            else:
                sys.exit("ERROR! Unrecognized detector '{}'. Exiting...".format(detector))
        else:
            sys.exit("ERROR! Unrecognized instrument '{}'. Exiting...".format(instrument))

        # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
        imgwcs = HSTWCS(imghdu, 1)
        fwhmpsf_pix = sourcecatalogdict[imgname]["params"]['fwhmpsf']/imgwcs.pscale #Convert fwhmpsf from arsec to pixels

        sourcecatalogdict[imgname]["catalog_table"] = amutils.generate_source_catalog(imghdu, fwhm=fwhmpsf_pix, **detector_pars)

        # write out coord lists to files for diagnostic purposes. Protip: To display the sources in these files in DS9,
        # set the "Coordinate System" option to "Physical" when loading the region file.
        imgroot = os.path.basename(imgname).split('_')[0]
        numSci = amutils.countExtn(imghdu)
        # Allow user to decide when and how to write out catalogs to files
        if output:
            for chip in range(1,numSci+1):
                regfilename = "{}_sci{}_src.reg".format(imgroot, chip)
                out_table = Table(sourcecatalogdict[imgname]["catalog_table"][chip])
                out_table.write(regfilename, include_names=["xcentroid", "ycentroid"], format="ascii.fast_commented_header")
                print("Wrote region file {}\n".format(regfilename))
        imghdu.close()
    return(sourcecatalogdict)


# ----------------------------------------------------------------------------------------------------------------------


def update_image_wcs_info(tweakwcs_output,imagelist):
    """Write newly computed WCS information to image headers

    Parameters
    ----------
    tweakwcs_output : list
        output of tweakwcs. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
        input image.

    imagelist : list
        list of valid processed images to be updated

    Returns
    -------
    Nothing!
    """
    imgctr = 0
    for item in tweakwcs_output:
        #print('YYYYYYYYY',item.wcs.pscale)
        if item.meta['chip'] == 1:  # to get the image name straight regardless of the number of chips
            image_name = imagelist[imgctr]
            if imgctr > 0: #close previously opened image
                print("CLOSE {}".format(hdulist[0].header['FILENAME'])) #TODO: Remove before deployment
                hdulist.flush()
                hdulist.close()
            hdulist = fits.open(image_name, mode='update')
            sciExtDict = {}
            for sciExtCtr in range(1, amutils.countExtn(hdulist) + 1): #establish correct mapping to the science extensions
                sciExtDict["{}".format(sciExtCtr)] = fileutil.findExtname(hdulist,'sci',extver=sciExtCtr)
            imgctr += 1
        updatehdr.update_wcs(hdulist, sciExtDict["{}".format(item.meta['chip'])], item.wcs, wcsname='TWEAKDEV', reusename=True, verbose=True) #TODO: May want to settle on a better name for 'wcsname'
        print()
    print("CLOSE {}".format(hdulist[0].header['FILENAME'])) #TODO: Remove before deployment
    hdulist.flush() #close last image
    hdulist.close()


# ======================================================================================================================

def interpret_fit_rms(tweakwcs_output, reference_catalog):
    """Interpret the FIT information to convert RMS to physical units

    Parameters
    ----------
    tweakwcs_output : list
        output of tweakwcs. Contains sourcelist tables, newly computed WCS info, etc. for every chip of every valid
        input image.  This list gets updated, in-place, with the new RMS values;
        specifically,

            * 'FIT_RMS': RMS of the separations between fitted image positions and reference positions
            * 'TOTAL_RMS': mean of the FIT_RMS values for all observations
            * 'NUM_FITS': number of images/group_id's with successful fits included in the TOTAL_RMS

        These entries are added to the 'tweakwcs_info' dictionary.

    reference_catalog : astropy.Table
        Table of reference source positions used for the fit

    Returns
    -------
    Nothing
    """
    # Start by collecting information by group_id
    group_ids = [info.meta['group_id'] for info in tweakwcs_output]
    group_dict = {'avg_RMS':None}
    obs_rms = []
    for group_id in group_ids:
        for item in tweakwcs_output:
            if item.meta['tweakwcs_info']['status'].startswith('FAILED'):
                continue
            if item.meta['group_id'] == group_id and \
               group_id not in group_dict:
                    group_dict[group_id] = {'ref_idx':None, 'FIT_RMS':None}
                    tinfo = item.meta['tweakwcs_info']
                    ref_idx = tinfo['fit_ref_idx']
                    group_dict[group_id]['ref_idx'] = ref_idx
                    ref_RA = reference_catalog[ref_idx]['RA']
                    ref_DEC = reference_catalog[ref_idx]['DEC']
                    img_coords = SkyCoord(tinfo['fit_RA'], tinfo['fit_DEC'],
                                          unit='deg',frame='icrs')
                    ref_coords = SkyCoord(ref_RA, ref_DEC, unit='deg',frame='icrs')
                    dra, ddec = img_coords.spherical_offsets_to(ref_coords)
                    ra_rms = np.std(dra.to(u.mas))
                    dec_rms = np.std(ddec.to(u.mas))
                    fit_rms = np.std(Angle(img_coords.separation(ref_coords), unit=u.mas)).value
                    group_dict[group_id]['FIT_RMS'] = fit_rms
                    group_dict[group_id]['RMS_RA'] = ra_rms
                    group_dict[group_id]['RMS_DEC'] = dec_rms
                    obs_rms.append(fit_rms)
    # Compute RMS for entire ASN/observation set
    total_rms = np.mean(obs_rms)
    #total_rms = np.sqrt(np.sum(np.array(obs_rms)**2))

    # Now, append computed results to tweakwcs_output
    for item in tweakwcs_output:
        group_id = item.meta['group_id']
        if group_id in group_dict:
            fit_rms = group_dict[group_id]['FIT_RMS']
            ra_rms = group_dict[group_id]['RMS_RA']
            dec_rms = group_dict[group_id]['RMS_DEC']
        else:
            fit_rms = None
            ra_rms = None
            dec_rms = None

        item.meta['tweakwcs_info']['FIT_RMS'] = fit_rms
        item.meta['tweakwcs_info']['TOTAL_RMS'] = total_rms
        item.meta['tweakwcs_info']['NUM_FITS'] = len(group_ids)
        item.meta['tweakwcs_info']['RMS_RA'] = ra_rms
        item.meta['tweakwcs_info']['RMS_DEC'] = dec_rms
        item.meta['tweakwcs_info']['catalog'] = reference_catalog.meta['catalog']


# ======================================================================================================================


if __name__ == '__main__':
    import argparse
    PARSER = argparse.ArgumentParser(description='Align images')
    PARSER.add_argument('raw_input_list', nargs='+', help='A space-separated list of fits files to align, or a simple '
                    'text file containing a list of fits files to align, one per line')

    PARSER.add_argument( '-a', '--archive', required=False,choices=['True','False'],default='False',help='Retain '
                    'copies of the downloaded files in the astroquery created sub-directories? Unless explicitly set, '
                    'the default is "False".')

    PARSER.add_argument( '-c', '--clobber', required=False,choices=['True','False'],default='False',help='Download and '
                    'overwrite existing local copies of input files? Unless explicitly set, the default is "False".')

    PARSER.add_argument( '-u', '--update_hdr_wcs', required=False,choices=['True','False'],default='False',help='Write '
                    'newly computed WCS information to image image headers? Unless explicitly set, the default is '
                    '"False".')
    ARGS = PARSER.parse_args()

    # Build list of input images
    input_list = []
    for item in ARGS.raw_input_list:
        if os.path.exists(item):
            with open(item, 'r') as infile:
                fileLines = infile.readlines()
            for fileLine in fileLines:
                input_list.append(fileLine.strip())
        else:
            input_list.append(item)

    archive = convert_string_tf_to_boolean(ARGS.archive)

    clobber = convert_string_tf_to_boolean(ARGS.clobber)

    update_hdr_wcs = convert_string_tf_to_boolean(ARGS.update_hdr_wcs)
    # Get to it!
    return_value = perform_align(input_list,archive,clobber,update_hdr_wcs)
