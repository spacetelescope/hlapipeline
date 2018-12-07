""" Utility to filter datasets which cannot be aligned 

The function analyze_data opens an input list containing FLT or FLC FITS file 
names in order to access the primary header data.  Based upon the values of specified 
FITS keywords, the function determines whether or not each file within this dataset 
can be reconciled against an astrometric catalog and, for multiple images, used to create 
a mosaic.

"""
from astropy.io import fits
from astropy.table import Table
import math

__all__ = ['analyze_data']

def analyze_data(inputFileList, **kwargs):
    """
    Determine if images within the dataset can be aligned
   
    Parameters
    ==========
    inputFileList: list
        List containing FLT or FLC filenames for all input images which comprise an associated 
        dataset where 'associated dataset' may be a single image, multiple images, an association, 
        or a number of associations

    Returns
    =======
    outputTable: object
        Astropy Table object containing data pertaining to the associated dataset, including 
        the doProcess bool

    Notes 
    =====
    The keyword/value pairs below define the "cannot process categories".
    OBSTYPE : not IMAGING
    MTFLAG : T
    SCAN-TYP : C or D (or !N)
    FILTER : G*, *POL*
    FILTER1 : G*, *POL*
    FILTER2 : G*, *POL*
    APERTURE : *GRISM*, G*-REF, RAMP, *POL*
    TARGNAME : DARK, TUNGSTEN, BIAS, FLAT, EARTH-CALIB, DEUTERIUM
    EXPOTIME : 0 

    Keywords only for WFC3 data: SCAN_TYP and FILTER
    Keywords only for ACS data: FILTER1 and FILTER2
    """
    OBSKEY = 'OBSTYPE'
    MTKEY  = 'MTFLAG'
    SCNKEY = 'SCAN_TYP'
    FILKEY = 'FILTER'
    FILKEY1 = 'FILTER1'
    FILKEY2 = 'FILTER2'
    APKEY  = 'APERTURE'
    TARKEY = 'TARGNAME'
    EXPKEY = 'EXPTIME'
    DATEKEY = 'DATE-OBS'

    filtNameList = [FILKEY1, FILKEY2]

    catalog = None
    foundSources = 0 
    rms_x = None
    rms_y = None
    matchSources = 0 
    isSuccess = False
    dateObs = None
    namesArray = ('imageName', 'instrument', 'detector', 'filter', 'aperture', 'obstype', 
            'subarray', 'dateObs', 'doProcess', 'catalog', '# found sources', '# match sources', 
            'rms_x', 'rms_y', 'isSuccess')
    dataType = ('S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'b', 'S20', 'b', 'S20', 'i4', 'i4', 'f8', 'f8', 'b')

    # Create an astropy table
    outputTable = Table(names=namesArray,dtype=dataType)

    # Loop over the list of images to determine viability for alignment processing
    # Capture the data characteristics before any evaluation so the information is
    # available for the output table regardless of which keyword is used to 
    # to determine the data is not viable for alignment.

    for inputFile in inputFileList:

        header_hdu  = 0
        header_data = fits.getheader(inputFile, header_hdu)

        obstype  = (header_data[OBSKEY]).upper()
        mtflag   = (header_data[MTKEY]).upper()

        # Keywords to use potentially for analysis
        instrume = (header_data['INSTRUME']).upper()
        detector = (header_data['DETECTOR']).upper()
        subarray = header_data['SUBARRAY']
        dataObs  = header_data['DATE-OBS']
        
        scan_typ = ''
        if instrume == 'WFC3':
            scan_typ = (header_data[SCNKEY]).upper()

        sfilter = ''
        if instrume == 'WFC3':
            sfilter  = (header_data[FILKEY]).upper()
        # Concatenate the two ACS filter names together with an underscore
        # If the filter name is blank, skip it
        if instrume == 'ACS':
            for filtname in filtNameList:

                # The filter keyword value could be zero or more blank spaces 
                # Strip the name down to non-zero leading or trailing blanks
                if len(header_data[filtname].upper().strip()) > 0:

                    # If the current filter variable already has some content,
                    # need to append an underscore before adding more text
                    if len(sfilter) > 0:
                        sfilter += '_'
                    sfilter += header_data[filtname].upper().strip()

        aperture = (header_data[APKEY]).upper()
        targname = (header_data[TARKEY]).upper()
        expotime = header_data[EXPKEY]

        # Determine if the image has one of these conditions.  The routine
        # will exit processing upon the first satisfied condition.

        noProcKey   = None
        noProcValue = None
        doProcess = True
        # Imaging vs spectroscopic or coronagraphic
        if obstype != 'IMAGING':
            noProcKey   = OBSKEY
            noProcValue = obstype 

        # Moving target
        elif mtflag == 'T':
            noProcKey   = MTKEY
            noProcValue = mtflag 

        # Bostrophidon without or with dwell (WFC3 only)
        elif any ([scan_typ == 'C', scan_typ == 'D']):
            noProcKey   = SCNKEY
            noProcValue = scan_typ

        # Filter which begins with !F (e.g., 'POL*' or 'G*')
        # The sfilter variable may be the concatenation of two filters (f160_clear) - 
        # look at the first character of each filter to see if the character is !F
        elif sfilter[0] != 'F' and sfilter[0] != '' and sfilter[0] != 'C' : 
            noProcKey   = FILKEY
            noProcValue = sfilter

        elif '_' in sfilter:
            pos = sfilter.index('_')
            pos += 1

            if sfilter[pos] != 'F' and sfilter[pos] != '' and sfilter[pos] != 'C' : 
                noProcKey   = FILKEY
                noProcValue = sfilter

        # Ramp, polarizer, or grism 
        elif any (x in aperture for x in ['RAMP', 'POL', 'GRISM', '-REF']):
            noProcKey   = APKEY
            noProcValue = aperture 

        # Calibration target
        elif any (x in targname for x in ['DARK', 'TUNG', 'BIAS', 'FLAT', 'DEUT', 'EARTH-CAL']):
            noProcKey   = TARKEY
            noProcValue = targname

        # Exposure time of zero
        elif math.isclose(expotime, 0.0, abs_tol=1e-5):
            noProcKey   = EXPKEY
            noProcValue = expotime 

        if (noProcKey is not None):
            doProcess = False

            # Issue message to log file - reports the first issue which makes the file ineligible
            # for alignment
            issue_msg(inputFile, noProcKey, noProcValue)


        # Populate a row of the table
        outputTable.add_row([inputFile, instrume, detector, sfilter, aperture, obstype, 
                             subarray, dateObs, doProcess, catalog, foundSources, matchSources, 
                             rms_x, rms_y, isSuccess])

    return(outputTable)


def issue_msg(filename, key, value):
    """ Generate a message for the output log indicating the file/association will not
        be processed as the characteristics of the data are known to be inconsistent
        with alignment.
    """

    print('\nDataset ' + filename + ' has (keyword = value) of (' + key + ' = ' + str(value) + ').')
    print('Dataset cannot be aligned.\n')

