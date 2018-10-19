import os
import pytest
from astropy.table import Table

from stwcs import updatewcs

from base_test import BaseHLATest
from hlapipeline import align_to_gaia

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

    def test_astroquery(self):
        """Verify that new astroquery interface will work"""
        self.curdir = os.getcwd()
        self.input_loc = ''

        filenames = self.get_input_file('ib6v06060')
        for infile in filenames:
            updatewcs.updatewcs(infile)

        output_shift_file = 'test_astroquery_shifts.txt'
        align_to_gaia.align(filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)
