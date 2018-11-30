.. highlight:: python
   :linenothreshold: 5

=======
Design
=======

The design of the HLA Pipeline code stems from a need to automate previously
manually-intensive processing for aligning HST data to astrometric coordinate frames
and then generating various products.  This document provides an overview of the
logic implemented within this code to achieve the best possible results in an
automated environment.


Image alignment
================
Image alignment processing focuses exclusively on doing everything possible to
correct the WCS specification for each HST observation so that it is as close to
the GAIA coordinate frame as possible.  There will be data which can never be
aligned to GAIA based on what was observed in the data itself (a posteriori
correction).  In those cases, the previously applied a priori corrections to the
guide star coordinates will remain unchanged.

The processing implemented to perform an a posteriori correction to the WCS
includes:

#. Interpret input data and optional parameters
#. Apply filter to input observations to insure that they meet minimum criteria
   for being able to be aligned using a posteriori method.  This filter will be
   implemented as :py:func:`filter.filter`.  These criteria will insure the image
   is NOT one of the following:

        * GRISM
        * moving target
        * spatial scan
        * any other types of observations which will not have measurable sources in the
          field-of-view like internal calibration observations (eg., DARKS, INTFLATs,...)

    #. If input fails this check, it simply returns without updating the WCS of any
       input data.

#.  Build WCS for full set of input observations using :py:mod:`stwcs` and :py:mod:`astropy`
#.  Retrieve list of astrometric sources from database using
    :py:func:`astrometry_utils.create_astrometric_catalog`

    * if too few from GAIADR2, try GSC24x or other database

    #. If still not enough reference sources available, return without updating WCS

#.  Extract catalog of observable sources from each input image using :py:mod:`astropy` :py:mod:`photutils`

    #. Detection code will evaluate DQ array and user-specified values to zero out any SCI array pixel
        which is flagged as 'BAD' in the DQ array and use zeroed-out SCI array for
        source detection
    #. Use of segmentation-based detection algorithm will allow for use of criteria
       to classify and subsequently ignore cosmic-rays in output detection catalog.

    * For WFC3/IR or ACS/SBC observations (data which do not include significant cosmic-rays),
      start with segmentation WITHOUT any filter with a high detection threshold (1000)
    * For all cosmic-ray affected observations, start with segmentation using Gaussian2DKernel
      and a high detection threshold

#.  Cross-match source catalog with astrometric reference source catalog using :py:mod:`tweakwcs`

    #. if no cross-matches, jump to iteration step below.

#.  Perform fit between source catalog and reference catalog using :py:class:`tweakwcs.TPMatch`
    with :py:func:`tweakwcs.tweak_wcs`

    #. The type of fit needs to be determined by the number of cross-matches

      * 3(?) or more matches, use full 'general' fit for offset, X/Y scale, and X/Y rotation
      * 1 or 2 matches, only perform 'offset' fit.

    #. If valid fit, archive current primary (a priori or OPUS) WCS using :py:mod:`stwcs.altwcs`,
       then apply fit to create updated WCS using :py:func:`drizzlepac.updatehdr`
    #. Apply same solution to both FLT and FLC images, if both types of images are present

.. note ::
  A **valid fit** should exhibit these qualities:

    * Consistent (within a few sigma) offsets/rotations for all observations in
      a single HST-pipeline-defined association
    * RMS is less than 0.5 pixels
    * rotation should be (much, much??) less than 0.01 degrees
    * scale should be within 1e-3 of 1. 

#. If no valid fit, iterate with the following changes:

    * return to source catalog generation step and derive a new source catalog
      using a lower threshold (100?), and try again.
    * Should no valid solution still be found:

        * return to source catalog generation step again
        * derive another source catalog which did NOT use the Gaussian2DKernel
          filter for most observations *or*
        * in the case of WFC3/IR, use the Gaussian2DKernel with segmentation, and try again.

    * If no valid solution can still be determined, return without updating the WCS

#. Upload new WCS to astrometry database using same code already implemented in :py:mod:`updatewcs`
