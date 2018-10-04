import pytest
from astropy.table import Table

from stwcs import updatewcs

from base_test import BaseHLATest
from hlapipeline import align_to_gaia

@pytest.mark.bigdata
class TestAlignMosaic(BaseHLATest):
    input_loc = 'mosaic_ngc188'
    ref_loc = ['truth']
    input_filenames = ['iaal01hxq_flc.fits', 'iaala3btq_flc.fits',
                       'iaal01hyq_flc.fits', 'iaala3bsq_flc.fits',
                       'j8boa1m8q_flc.fits', 'j8boa1m4q_flc.fits',
                       'j8boa1maq_flc.fits', 'j8boa1m6q_flc.fits']

    def test_align_ngc188(self):

        for file in self.input_filenames:
            self.get_input_file(file, docopy=True)
            updatewcs.updatewcs(file)

        output_shift_file = 'test_mosaic_ngc188_shifts.txt'
        align_to_gaia.align(self.input_filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)

    def test_align_47tuc(self):
        self.input_filenames = ['ib6v06c4q_flc.fits','ib6v06c7q_flc.fits',
                                'ib6v25aqq_flc.fits','ib6v25atq_flc.fits',
                                'jddh02gjq_flc.fits','jddh02glq_flc.fits',
                                'jddh02goq_flc.fits']
        self.input_loc = 'mosaic_47tuc'

        for file in self.input_filenames:
            self.get_input_file(file, docopy=True)
            updatewcs.updatewcs(file)

        output_shift_file = 'test_mosaic_47tuc_shifts.txt'
        align_to_gaia.align(self.input_filenames, shift_name=output_shift_file)

        shift_file = Table.read(output_shift_file, format='ascii')
        rms_x = max(shift_file['col6'])
        rms_y = max(shift_file['col7'])

        assert (rms_x <= 0.25 and rms_y <= 0.25)
