import sys
import traceback
import os
import pytest
import numpy
from astropy.table import Table

from stwcs import updatewcs

from base_test import BaseHLATest
from hlapipeline import align_to_gaia
import hlapipeline.utils.catalog_utils as catutils

@pytest.mark.bigdata
class TestAlignMosaic(BaseHLATest):
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

    def test_align_ngc188(self):
        """ Verify whether NGC188 exposures can be aligned to an astrometric standard.

        Characeteristics of this test:
          * Input exposures include both ACS and WFC3 images of the same general field-of-view
            of NGC188 suitable for creating a combined mosaic using both instruments.
        """
        self.input_loc = 'mosaic_ngc188'
        input_filenames = ['iaal01hxq_flc.fits', 'iaala3btq_flc.fits',
                            'iaal01hyq_flc.fits', 'iaala3bsq_flc.fits',
                            'j8boa1m8q_flc.fits', 'j8boa1m4q_flc.fits',
                            'j8boa1maq_flc.fits', 'j8boa1m6q_flc.fits']
        self.output_shift_file = 'test_mosaic_ngc188_shifts.txt'
        shift_file = self.run_align(input_filenames)

        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    def test_align_47tuc(self):
        """ Verify whether 47Tuc exposures can be aligned to an astrometric standard.

        Characeteristics of this test:
          * Input exposures include both ACS and WFC3 images of the same general field-of-view
            of 47Tuc suitable for creating a combined mosaic using both instruments.
        """
        self.input_loc = 'mosaic_47tuc'
        input_filenames = ['ib6v06c4q_flc.fits','ib6v06c7q_flc.fits',
                                'ib6v25aqq_flc.fits','ib6v25atq_flc.fits',
                                'jddh02gjq_flc.fits','jddh02glq_flc.fits',
                                'jddh02goq_flc.fits']
        self.output_shift_file = 'test_mosaic_47tuc_shifts.txt'
        shift_file = self.run_align(input_filenames)

        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    @pytest.mark.parametrize("input_filenames", [['j8ura1j1q_flt.fits','j8ura1j2q_flt.fits',
                                                  'j8ura1j4q_flt.fits','j8ura1j6q_flt.fits',
                                                  'j8ura1j7q_flt.fits','j8ura1j8q_flt.fits',
                                                  'j8ura1j9q_flt.fits','j8ura1jaq_flt.fits',
                                                  'j8ura1jbq_flt.fits','j8ura1jcq_flt.fits',
                                                  'j8ura1jdq_flt.fits','j8ura1jeq_flt.fits',
                                                  'j8ura1jfq_flt.fits','j8ura1jgq_flt.fits',
                                                  'j8ura1jhq_flt.fits','j8ura1jiq_flt.fits',
                                                  'j8ura1jjq_flt.fits','j8ura1jkq_flt.fits'],
                                                 ['j92c01b4q_flc.fits','j92c01b5q_flc.fits',
                                                  'j92c01b7q_flc.fits','j92c01b9q_flc.fits'],
                                                 ['jbqf02gzq_flc.fits', 'jbqf02h5q_flc.fits',
                                                  'jbqf02h7q_flc.fits', 'jbqf02hdq_flc.fits',
                                                  'jbqf02hjq_flc.fits', 'jbqf02hoq_flc.fits',
                                                  'jbqf02hqq_flc.fits', 'jbqf02hxq_flc.fits',
                                                  'jbqf02i3q_flc.fits', 'jbqf02i8q_flc.fits',
                                                  'jbqf02iaq_flc.fits'],
                                                 ['ib2u12kaq_flt.fits', 'ib2u12keq_flt.fits',
                                                  'ib2u12kiq_flt.fits', 'ib2u12klq_flt.fits'],
                                                 ['ibjt01a1q_flc.fits', 'ibjt01a8q_flc.fits',
                                                  'ibjt01aiq_flt.fits', 'ibjt01amq_flt.fits',
                                                  'ibjt01aqq_flt.fits', 'ibjt01auq_flt.fits',
                                                  'ibjt01yqq_flc.fits', 'ibjt01z0q_flc.fits',
                                                  'ibjt01zwq_flc.fits', 'ibjt01a4q_flc.fits',
                                                  'ibjt01acq_flc.fits', 'ibjt01akq_flt.fits',
                                                  'ibjt01aoq_flt.fits', 'ibjt01asq_flt.fits',
                                                  'ibjt01avq_flt.fits', 'ibjt01yuq_flc.fits',
                                                  'ibjt01ztq_flc.fits'],
                                                 ['ibnh02coq_flc.fits','ibnh02cmq_flc.fits',
                                                  'ibnh02c7q_flc.fits','ibnh02c5q_flc.fits',
                                                  'ibnh02cpq_flc.fits','ibnh02c9q_flc.fits',
                                                  'ibnh02bfq_flc.fits','ibnh02beq_flc.fits']])
    def test_align_single_visits(self,input_filenames):
        """ Verify whether single-visit exposures can be aligned to an astrometric standard.

        Characteristics of these tests:
          * Input exposures include exposures from a number of single visit datasets to explore what impact differing
            observing modes (differing instruments, detectors, filters, subarray size, etc.) have on astrometry.

        The following datasets are used in these tests:

            * ACS dataset 10048_a1: 2x F344N, 1x F435W, 1x F475W, 2x F502N, 2x F550M, 1x F555W, 1x F606W, 1x F625W,
              2x F658N, 1x F775W, 1x F814W, 1x F850LP, and 2x F892N ACS/HRC images
            * ACS dataset 10265_01: 4x F606W full-frame ACS/WFC images
            * ACS dataset 12580_02: 5x F475W & 6x F814W ACS/WFC images
            * WFC3 dataset 11663_12: 4x F160W full-frame WFC3/IR images
            * WFC3 dataset 12219_01: 8x F160W full-frame WFC3/IR images, 9x F336W full-frame WFC3/UVIS images
            * WFC3 dataset 12379_02: 4X F606W, 4x F502N full-frame WFC3/UVIS images

        """
        self.input_loc = 'base_tests'

        self.curdir = os.getcwd()
        for infile in input_filenames:
            self.get_input_file(infile, docopy=True)
            updatewcs.updatewcs(infile)

        output_shift_file = 'test_single_visits_shifts.txt'
        align_to_gaia.align(input_filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    def test_astroquery(self):
        """Verify that new astroquery interface will work"""
        self.curdir = os.getcwd()
        self.input_loc = ''

        shift_file = self.run_align('ib6v06060')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    @pytest.mark.xfail
    @pytest.mark.slow
    def test_align_randomFields(self):
        """ Wrapper to set up the test for aligning a large number of randomly
            selected fields (aka datasets) from a input ascii file (CSV).

            The wrapper provides the parameter settings for the underlying test,
            as well as implements the criterion for the overall success or failure
            of the test.
        """
        inputListFile = 'ACSList50.csv'

        # Desired number of random entries for testing
        inputNumEntries = 50

        # Seed for random number generator
        inputSeedValue = 1

        # Obtain the full path to the file containing the dataset field names
        self.input_loc  = 'master_lists'
        self.curdir     = os.getcwd()
        input_file_path = self.get_data(inputListFile)

        # Randomly select a subset of field names (each field represented by a row) from
        # the master CSV file and return as an Astropy table
        randomCandidateTable = catutils.randomSelectFromCSV(input_file_path[0],
            inputNumEntries, inputSeedValue)

        # Invoke the methods which will handle acquiring/downloading the data from
        # MAST and perform the alignment
        percentSuccess = 0.0
        try:
            percentSuccess = self.align_randomFields (randomCandidateTable)
        except Exception:
            pass

        assert(percentSuccess >= 0.70)

    @pytest.mark.xfail
    def test_align_fewRandomFields(self):
        """ Wrapper to set up the test for aligning a *FEW* randomly
            selected fields (aka datasets) from a input ascii file (CSV).

            The wrapper provides the parameter settings for the underlying test,
            as well as implements the criterion for the overall success or failure
            of the test.

            This routine is strictly for realistic testing on a small number of
            datasets.
        """
        inputListFile = 'ACSList5.csv'

        # Desired number of random entries for testing
        inputNumEntries = 5

        # Seed for random number generator
        inputSeedValue = 1

        # Obtain the full path to the file containing the dataset field names
        self.input_loc  = 'master_lists'
        self.curdir     = os.getcwd()
        input_file_path = self.get_data(inputListFile)

        # Randomly select a subset of field names (each field represented by a row) from
        # the master CSV file and return as an Astropy table
        randomCandidateTable = catutils.randomSelectFromCSV(input_file_path[0],
            inputNumEntries, inputSeedValue)

        # Invoke the methods which will handle acquiring/downloading the data from
        # MAST and perform the alignment
        percentSuccess = 0.0
        try:
            percentSuccess = self.align_randomFields (randomCandidateTable)
        except Exception:
            pass

        assert(percentSuccess >= 0.70)

    def align_randomFields(self, randomTable):
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
        for dataset in dataset_list:

           try:
               shift_file = self.run_align([dataset])
               x_shift = numpy.alltrue(numpy.isnan(shift_file['col2']))
               rms_x = max(shift_file['col6'])
               rms_y = max(shift_file['col7'])

               if not x_shift and ((rms_x <= 0.25) and (rms_y <= 0.25)):
                   numSuccess += 1

           # Catch anything that happens as this dataset will be considered a failure, but
           # the processing of datasets should continue.  Generate sufficient output exception
           # information so problems can be addressed.
           except Exception:
               exc_type, exc_value, exc_tb = sys.exc_info()
               traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
               continue

        # Determine the percent success over all datasets processed
        percentSuccess = numSuccess/numAllDatasets
        print('Number of successful tests: ', numSuccess, ' Total number of tests: ', numAllDatasets, ' Percent success: ', percentSuccess*100.0)

        return percentSuccess
#
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
