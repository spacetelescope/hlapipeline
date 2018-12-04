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

#def analyze_data(inputFileList, optionalCols = None):
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

    The keyword/value pairs below define the "cannot process categories".
    OBSTYPE: not IMAGING
    MTFLAG: T
    SCAN-TYP: C or D (or !N)
    FILTER: G*, *POL*
    APERTURE: *GRISM*, G*-REF, RAMP, *POL*
    TARGNAME: DARK, TUNGSTEN, BIAS, FLAT, EARTH-CALIB, DEUTERIUM
    EXPOTIME: 0 
    """
    OBSKEY = 'OBSTYPE'
    MTKEY  = 'MTFLAG'
    SCNKEY = 'SCAN_TYP'
    FILKEY = 'FILTER'
    APKEY  = 'APERTURE'
    TARKEY = 'TARGNAME'
    EXPKEY = 'EXPTIME'
    DATEKEY = 'DATE-OBS'

    # Create an astropy table
    #optionalCols = kwargs.get('optionalCols', None)

    catalog = None
    foundSources = 0 
    rms_x = None
    rms_y = None
    matchSources = 0 
    isSuccess = False
    #isSuccess = 0
    dateObs = None
    namesArray = ('imageName', 'instrument', 'detector', 'filter', 'aperture', 'obstype', 
            'subarray', 'dateObs', 'doProcess', 'catalog', '# found sources', '# match sources', 
            'rms_x', 'rms_y', 'isSuccess')
    #dataType = ('S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'i1', 'S20', 'i1', 'S20', 'i4', 'i4', 'f8', 'f8', 'i1')
    dataType = ('S20', 'S20', 'S20', 'S20', 'S20', 'S20', 'b', 'S20', 'b', 'S20', 'i4', 'i4', 'f8', 'f8', 'b')

    outputTable = Table(names=namesArray,dtype=dataType)

    # Loop over the list of images to determine viability for alignment processing
    for inputFile in inputFileList:

        header_hdu  = 0
        header_data = fits.getheader(inputFile, header_hdu)

        obstype  = (header_data[OBSKEY]).upper()
        mtflag   = (header_data[MTKEY]).upper()
        scan_typ = (header_data[SCNKEY]).upper()

        sfilter  = (header_data[FILKEY]).upper()
        aperture = (header_data[APKEY]).upper()
        targname = (header_data[TARKEY]).upper()
        expotime = header_data[EXPKEY]

        # Keywords to use potentially for analysis
        instrume = (header_data['INSTRUME']).upper()
        detector = (header_data['DETECTOR']).upper()
        subarray = header_data['SUBARRAY']
        dataObs  = header_data['DATE-OBS']

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

        # Bostrophidon without or with dwell
        elif any ([scan_typ == 'C', scan_typ == 'D']):
            noProcKey   = SCNKEY
            noProcValue = scan_typ

        # Filter which begins with !F (e.g., 'POL*' or 'G*')
        elif sfilter[0] != 'F': 
            noProcKey   = FILKEY
            noProcValue = sfilter

        # Ramp, polarizer, or grism 
        elif any (['RAMP' in aperture, 'POL' in aperture, 'GRISM' in aperture, 
                  '-REF' in aperture]):
            noProcKey   = APKEY
            noProcValue = aperture 

        # Calibration target
        elif any (['DARK' in targname, 'TUNG' in targname, 'BIAS' in targname, 
                  'FLAT' in targname, 'DEUT' in targname, 'EARTH-CAL' in targname]):
            noProcKey   = TARKEY
            noProcValue = targname

        # Exposure time of zero
        elif math.isclose(expotime, 0.0, abs_tol=1e-5):
            noProcKey   = EXPKEY
            noProcValue = expotime 

        if (noProcKey is not None):
            #doProcess = 0
            doProcess = False

        # Issue message to log file 
        #issue_msg(inputFileList, noProcKey, noProcValue)

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

    print('Dataset ' + filename + ' has keyword = value of ' + key + ' = ' + str(value) + '.\n')
    print('Dataset cannot be aligned.\n')

