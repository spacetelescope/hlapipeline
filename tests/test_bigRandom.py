import sys
import traceback
import os
import glob
import pytest
import numpy as np
from astropy.table import Table, vstack
from astropy.io import ascii

from base_test import BaseHLATest
from hlapipeline import alignimages
import hlapipeline.utils.catalog_utils as catutils

@pytest.mark.bigdata
class TestRandomAlignMosaic(BaseHLATest):
    """ Tests which validate whether mosaics can be aligned to an astrometric standard.

        Characeteristics of these tests:
          * A single astrometric catalog was obtained with both GAIA and non-GAIA
            (PanSTARRS?) sources for the entire combined field-of-view using the GSSS
            server.
              * These tests assume combined input exposures for each test have enough
                astrometric sources from the external catalog to perform a high quality
                fit.
          * This test only determines the fit between sources
            extracted from the images by Tweakreg and the source positions included in
            the astrometric catalog.
          * The WCS information for the input exposures do not get updated in this test.
          * No mosaic gets generated.

        Success Criteria:
          * Success criteria hard-coded for this test represents 10mas RMS for the
            WFC3 images based on the fit of source positions to the astrometric catalog
            source positions.
              * RMS values are extracted from optional shiftfile output from `tweakreg`
              * Number of stars used for the fit and other information is not available
                with the current version of `tweakreg`.
    """

    ref_loc = ['truth']

    #@pytest.mark.xfail
    #@pytest.mark.slow
    def test_random_align(self):
        """ Wrapper to set up the test for aligning a large number of randomly
            selected fields (aka datasets) from a input ascii file (CSV)
            using the new algorithms developed for producing the HLA products.
        """
        inputListFiles = ['ACSList50R_v2.csv', 'WFC3List50R_v2.csv']
        #inputListFiles = ['ACSList2.csv']

        # Desired number of random entries for testing from each input CSV
        inputNumEntries = 50
        #inputNumEntries = 3

        # Seed for random number generator
        inputSeedValue = 1

        # Obtain the full path to the file containing the dataset field names
        self.input_loc  = 'master_lists'
        self.curdir     = os.getcwd()

        # Loop over the CSV files so the final candidate table is organized according to the instrument
        randomTable = Table()
        for inputList in inputListFiles:
            input_file_path = self.get_data(inputList)

            # Randomly select a subset of field names (each field represented by a row) from
            # the master CSV file and return as an Astropy table
            #
            # NOTE: This is not really necessary as the input files have already been 
            # randomized.  However, the function does return an Astropy table.
            randomCandidateTable = catutils.randomSelectFromCSV(input_file_path[0],
                inputNumEntries, inputSeedValue)

            randomTable = vstack([randomTable, randomCandidateTable])

        # Invoke the methods which will handle acquiring/downloading the data from
        # MAST and perform the alignment
        percentSuccess = 0.0
        try:
            percentSuccess = self.random_align(randomTable)
        except Exception:
            pass

        assert(percentSuccess >= 0.70)

    def random_align(self, randomTable):
        """ Process randomly selected fields (aka datasets) stored in an Astropy table.

            Each field is used as input to determine if it can be aligned to an
            astrometric standard.  The success or fail status for each test is retained
            as the overall success or fail statistic is the necessary output from
            this test.
        """

        numSuccess = 0
        numAllDatasets = 0

        # Read the table and extract a list of each dataset name in IPPSSOOT format
        # which is either an association ID or an individual filename
        dataset_list = get_dataset_list(randomTable)

        numAllDatasets = len(dataset_list)

        # Process the dataset names in the list
        #
        # If the dataset name represents an association ID, the multiplicity
        # of images within the association need to be processed.  Otherwise,
        # the dataset is a single image.
        #
        # If the "alignment" of a field/dataset fails for any reason, trap
        # the exception and keep going.
        allDatasetTable = Table()
        datasetKey = -1
        for dataset in dataset_list:

            datasetKey += 1
            outputName = dataset + '.ecsv'

            print("TEST_RANDOM. Dataset: ", dataset, ' DatasetKey: ', datasetKey)
            try:
                datasetTable = alignimages.perform_align([dataset])

                # Filtered datasets
                if datasetTable['doProcess'].sum() == 0:
                    print("TEST_RANDOM. Filtered Dataset: ", dataset, "\n")
                    numAllDatasets -= 1;
                # Datasets to process
                elif datasetTable['doProcess'].sum() > 0:
                    # Determine images in dataset to be processed and the number of images
                    # This is in case an image was filtered out (e.g., expotime = 0)
                    index = np.where(datasetTable['doProcess']==1)[0]
                    sumOfStatus = datasetTable['status'][index].sum()

                    # Update the table with the datasetKey which is really just a counter
                    datasetTable['datasetKey'][:] = datasetKey
                    datasetTable['completed'][:] = True
                    datasetTable.write(outputName, format='ascii.ecsv')
                    datasetTable.pprint(max_width=-1)
   
                    # Successful datasets
                    if (sumOfStatus == 0):
                        print("TEST_RANDOM. Successful Dataset: ", dataset, "\n")
                        numSuccess += 1
                    # Unsuccessful datasets
                    else:
                        print("TEST_RANDOM. Unsuccessful Dataset: ", dataset, "\n")

                # Append the latest dataset table to the summary table 
                allDatasetTable = vstack([allDatasetTable, datasetTable])

                # Perform some clean up
                if os.path.exists('ref_cat.ecsv'): os.remove('ref_cat.ecsv')
                if os.path.exists('refcatalog.cat'): os.remove('refcatalog.cat')
                for f in sorted(glob.glob('*fl?.fits')):
                    os.remove(f)

            # Catch anything that happens as this dataset will be considered a failure, but
            # the processing of datasets should continue.  Generate sufficient output exception
            # information so problems can be addressed.
            except Exception:
           
                exc_type, exc_value, exc_tb = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
                print("TEST_RANDOM. Exception Dataset: ", dataset, "\n")
                continue

        # Write out the table
        allDatasetTable.write('resultsBigTest.ecsv', format='ascii.ecsv')
        #allDatasetTable.pprint(max_width=-1)

        # Determine the percent success over all datasets processed
        percentSuccess = numSuccess/numAllDatasets
        print('TEST_RANDOM. Number of successful tests: ', numSuccess, ' Total number of tests: ', numAllDatasets, ' Percent success: ', percentSuccess*100.0)
 
        return percentSuccess

def get_dataset_list(tableName):
    """ Standalone function to read the Astropy table and get the dataset names

    Parameters
    ==========
    tableName : str
        Filename of the input master CSV file containing individual
        images or association names, as well as observational
        information regarding the images

    Returns
    =======
    datasetNames: list
        List of individual image or association base (IPPSSOOT) names
    """

    #dataFromTable = Table.read(filename, format='ascii')
    datasetIDs = tableName['observationID']
    asnIDs     = tableName['asnID']

    datasetNames = []

    # Determine if the data is part of an association or is an individual image
    for imgid,asnid in zip(datasetIDs,asnIDs):

        # If the asnID is the string NONE, this is an individual image,
        # and it is necessary to get the individual image dataset name.
        # Otherwise, this is an association dataset, so just add the asnID.
        if (asnid.upper() == "NONE"):
            datasetNames.append(imgid)
        else:
            datasetNames.append(asnid)

    return datasetNames
