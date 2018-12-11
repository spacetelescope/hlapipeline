#!/usr/bin/env python

"""This script is a modernized replacement of tweakreg.

"""

from astropy.io import fits
from astropy.table import Table
from collections import OrderedDict
import glob
import numpy as np
import os
import pdb
from stwcs.wcsutil import HSTWCS
import sys
import tweakwcs
from utils import astrometric_utils as amutils
from utils import astroquery_utils as aqutils
from utils import filter

MIN_CATALOG_THRESHOLD = 3
MIN_OBSERVABLE_THRESHOLD = 10
MIN_CROSS_MATCHES = 3
MIN_FIT_MATCHES = 6
MAX_FIT_RMS = 1.0

# Module-level dictionary contains instrument/detector-specific parameters used later on in the script.
detector_specific_params = {"acs":
                                {"hrc":
                                     {"fwhmpsf": 0.073,
                                      "classify": True,
                                      "threshold": None},
                                 "sbc":
                                     {"fwhmpsf": 0.065,
                                      "classify": False,
                                      "threshold": 10},
                                 "wfc":
                                     {"fwhmpsf": 0.076,
                                      "classify": True,
                                      "threshold": None}},  # TODO: Verify ACS fwhmpsf values
                            "wfc3":
                                {"ir":
                                     {"fwhmpsf": 0.14,
                                      "classify": False,
                                      "threshold": None},
                                 "uvis":
                                     {"fwhmpsf": 0.076,
                                      "classify": True,
                                      "threshold": None}}}
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
        if input_item.endswith("0"): #asn table
            totalInputList += aqutils.retrieve_observation(input_item,**pars)

        else: #single file rootname.
            fitsfilename = glob.glob("{}_flc.fits".format(input_item))
            if not fitsfilename:
                fitsfilename = glob.glob("{}_flt.fits".format(input_item))
            fitsfilename = fitsfilename[0]

            if not os.path.exists(fitsfilename):
                imghdu = fits.open(fitsfilename)
                imgprimaryheader = imghdu[0].header
                try:
                    asnid = imgprimaryheader['ASN_ID'].strip().lower()
                except:
                    asnid = 'NONE'
                if asnid[0] in ['i','j']:
                    totalInputList += aqutils.retrieve_observation(asnid,**pars)
                else:
                    totalInputList += aqutils.retrieve_observation(input_item, **pars) #try with ippssoot instead

            else: totalInputList.append(fitsfilename)
    print("TOTAL INPUT LIST: ",totalInputList)
    # TODO: add trap to deal with non-existent (incorrect) rootnames
    # TODO: Address issue about how the code will retrieve association information if there isn't a local file to get 'ASN_ID' header info
    return(totalInputList)


def perform_align(input_list):
    """Main calling function.

    Parameters
    ----------
    input_list : list
        List of one or more IPPSSOOTs (rootnames) to align.

    Returns
    -------
    int value 0 if successful, int value 1 if unsuccessful

    """

    # Define astrometric catalog list in priority order
    catalogList = ['GAIADR2', 'GSC241']
    numCatalogs = len(catalogList)

    # 1: Interpret input data and optional parameters
    print("-------------------- STEP 1: Get data --------------------")
    imglist = check_and_get_data(input_list)
    print("\nSUCCESS")

    # 2: Apply filter to input observations to insure that they meet minimum criteria for being able to be aligned
    print("-------------------- STEP 2: Filter data --------------------")
    filteredTable = filter.analyze_data(imglist)

    # Check the table to determine if there is any viable data to be aligned.  The
    # 'doProcess' column (bool) indicates the image/file should or should not be used
    # for alignment purposes.
    if filteredTable['doProcess'].sum() == 0:
        print("No viable images in filtered table - no processing done.\n")
        return(1)

    # Get the list of all "good" files to use for the alignment
    processList = filteredTable['imageName'][np.where(filteredTable['doProcess'])]
    processList = list(processList) #Convert processList from numpy list to regular python list
    print("\nSUCCESS")

    # 3: Build WCS for full set of input observations
    print("-------------------- STEP 3: Build WCS --------------------")
    refwcs = amutils.build_reference_wcs(processList)
    print("\nSUCCESS")

    # 4: Retrieve list of astrometric sources from database
    # While loop to accommodate using multiple catalogs
    doneFitting = False
    catalogIndex = 0
    while not doneFitting:
        print("-------------------- STEP 4: Detect astrometric sources --------------------")
        print("Astrometric Catalog: ",catalogList[catalogIndex])
        reference_catalog = generate_astrometric_catalog(processList, catalog=catalogList[catalogIndex])
        # The table must have at least MIN_CATALOG_THRESHOLD entries to be useful
        if len(reference_catalog) >= MIN_CATALOG_THRESHOLD:
            print("\nSUCCESS")
        # 5: Extract catalog of observable sources from each input image
            print("-------------------- STEP 5: Source finding --------------------")
            if catalogIndex == 0:
                extracted_sources = generate_source_catalogs(processList)
                for imgname in extracted_sources.keys():
                    table=extracted_sources[imgname]["catalog_table"]
                    # The catalog of observable sources must have at least MIN_OBSERVABLE_THRESHOLD entries to be useful
                    total_num_sources = 0
                    for chipnum in table.keys():
                        total_num_sources += len(table[chipnum])
                    if total_num_sources < MIN_OBSERVABLE_THRESHOLD:
                        print("Not enough sources ({}) found in image {}".format(total_num_sources,imgname))
                        return(1)
            # Convert input images to tweakwcs-compatible NDData objects and
            # attach source catalogs to them.
            imglist = []
            for group_id, image in enumerate(processList):
                imglist.extend(amutils.build_nddata(image, group_id,
                                                    extracted_sources[image]['catalog_table']))
            print("\nSUCCESS")

        # 6: Cross-match source catalog with astrometric reference source catalog, Perform fit between source catalog and reference catalog
            print("-------------------- STEP 6: Cross matching and fitting --------------------")
            # Specify matching algorithm to use
            match = tweakwcs.TPMatch(searchrad=250, separation=0.1,
                                     tolerance=250, use2dhist=False)
            # Align images and correct WCS
            tweakwcs.tweak_image_wcs(imglist, reference_catalog, match=match)

            tweakwcs_info_keys = OrderedDict(imglist[0].meta['tweakwcs_info']).keys()
            imgctr=0
            for item in imglist:
                max_rms_val = max(item.meta['tweakwcs_info']['rms'])
                num_xmatches = item.meta['tweakwcs_info']['nmatches']
                # print fit params to screen
                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIT PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                if item.meta['chip'] == 1:
                    image_name = processList[imgctr]
                    imgctr += 1
                print("image: {}".format(image_name))
                print("chip: {}".format(item.meta['chip']))
                print("group_id: {}".format(item.meta['group_id']))
                for tweakwcs_info_key in tweakwcs_info_keys:
                    if not tweakwcs_info_key.startswith("matched"):
                        print("{} : {}".format(tweakwcs_info_key,item.meta['tweakwcs_info'][tweakwcs_info_key]))
                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

                if num_xmatches < MIN_CROSS_MATCHES:
                    if catalogIndex < numCatalogs-1:
                        print("Not enough cross matches found between astrometric catalog and sources found in images")
                        print("Try again with the next catalog")
                        catalogIndex += 1
                        break
                    else:
                        print("Not enough cross matches found in any catalog - no processing done.")
                        return(1)
                elif max_rms_val > MAX_FIT_RMS:
                    if catalogIndex < numCatalogs-1:
                        print("Fit RMS value(s) X_rms= {}, Y_rms = {} greater than the maximum threshold value {}.".format(item.meta['tweakwcs_info']['rms'][0], item.meta['tweakwcs_info']['rms'][1],MAX_FIT_RMS))
                        print("Try again with the next catalog")
                        catalogIndex += 1
                        break
                    else:
                        print("Fit RMS values too large using any catalog - no processing done.")
                        return(1)
                else:
                    print("Fit calculations successful.")
        else:
            if catalogIndex < numCatalogs-1:
                print("Not enough sources found in catalog" + catalogList[catalogIndex])
                print("Try again with the next catalog")
                catalogIndex += 1
            else:
                print("Not enough sources found in any catalog - no processing done.")
                return(1)
    print("\nSUCCESS")
    
    # 7: Write new fit solution to input image headers

    return (0)


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

    # write catalog to ascii text file
    catalog_fileanme = "refcatalog.cat"
    out_catalog.write(catalog_fileanme, format="ascii.fast_commented_header")

    print("Wrote reference catalog {}.".format(catalog_fileanme))
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
                sourcecatalogdict[imgname]["params"] = detector_pars
                # to allow generate_source_catalog to get detector specific parameters
                pars.update(detector_pars)
            else:
                sys.exit("ERROR! Unrecognized detector '{}'. Exiting...".format(detector))
        else:
            sys.exit("ERROR! Unrecognized instrument '{}'. Exiting...".format(instrument))

        # Identify sources in image, convert coords from chip x, y form to reference WCS sky RA, Dec form.
        imgwcs = HSTWCS(imghdu, 1)
        fwhmpsf_pix = sourcecatalogdict[imgname]["params"]['fwhmpsf']/imgwcs.pscale
        sourcecatalogdict[imgname]["catalog_table"] = amutils.generate_source_catalog(imghdu, fwhm=fwhmpsf_pix, **pars)

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
# ======================================================================================================================


if __name__ == '__main__':
    import argparse
    PARSER = argparse.ArgumentParser(description='Align images')
    PARSER.add_argument('raw_input_list', nargs='+', help='A space-separated list of fits files to align, or a simple text '
                                                     'file containing a list of fits files to align, one per line')
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


    # Get to it!
    return_value = perform_align(input_list)

