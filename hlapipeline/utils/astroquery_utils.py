"""Wrappers for astroquery-related functionality"""
import shutil
import os

from astroquery.mast import Observations


def retrieve_observation(obsid, suffix=['FLC']):
    """Simple interface for retrieving an observation from the MAST archive

    If the input obsid is for an association, it will request all members with
    the specified suffixes.

    Parameters
    -----------
    obsid : string
        ID for observation to be retrieved from the MAST archive.  Only the
        IPPSSOOT (rootname) of exposure or ASN needs to be provided; eg., ib6v06060.

    suffix : list
        List containing suffixes of files which should be requested from MAST.

    path : string
        Directory to use for writing out downloaded files.  If `None` (default),
        the current working directory will be used.

    """
    local_files = []
    obsTable = Observations.query_criteria(obs_id=obsid)
    # Catch the case where no files are found for download
    if len(obsTable) == 0:
        print("WARNING: Query for {} returned NO RESULTS!".format(obsid))
        return local_files

    dpobs = Observations.get_product_list(obsTable)
    dataProductsByID = Observations.filter_products(dpobs,
                                              productSubGroupDescription=suffix,
                                              extension='fits',
                                              mrp_only=False)
    manifest = Observations.download_products(dataProductsByID, mrp_only=False)

    download_dir = None
    for file in manifest['Local Path']:
        # Identify what sub-directory was created by astroquery for the download
        if download_dir is None:
            file_path = file.split(os.sep)
            file_path.remove('.')
            download_dir = file_path[0]
        # Move downloaded file to current directory
        local_file = os.path.abspath(os.path.basename(file))
        shutil.move(file, local_file)
        # Record what files were downloaded and their current location
        local_files.append(local_file)
    # Remove astroquery created sub-directories
    shutil.rmtree(download_dir)

    return local_files
