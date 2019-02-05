"""
Copyright (C) University of Adelaide, Australia - All Rights Reserved
This source code is protected under proprietary reference only license.
Please refer to the accompanying LICENSE file for further information.
"""

import glob
import os
import sys
from numpy import *
from numbers import Number
# from scipy import *

__blur_cache = {}

class HS_Image:
    """A class to interact with the images taken from a camera.

    ``HS_Image`` is a class that is used to read the image frames and to process
    them prior to centroiding.

    **Instantiation**

    The class has two main options to read the image frames, specified by the
    instance variable ``location_type``.

    If ``location_type`` is set to 'file', then ``no_of_frames`` is set to 1, and
    the instance is supposed to interact with a single file. The instance variable
    ``location`` should then be the corresponding full-path file name.

    If ``location_type`` is set to 'folder', more than one image frame file can be
    read and averaged to produce an image. Basic preparation for this process is
    performed by the method :func:`~HS_Image.set_acquisition_parameters`.

    An example way to use the class in 'folder' mode is, with ``location`` properly set
    to a folder in which the image frame files are stored, to set ``fprefix`` (e.g.,
    'image') to be the file prefix of the files to be read, and set ``no_of_frames`` to
    a positive *integer* value (e.g., 10).

    Alternatively one can also set ``first_frame`` and ``last_frame`` as positive
    *integers*, ensuring that ``first_frame`` is smaller than ``last_frame``. The
    instance will then read and average that many files (from ``first_frame`` to
    ``last_frame``).

    **Instance Variables**

    - ``location_type``: the type of the ``location``, either 'file' or 'folder'

      :Type: *string*
      :Required by: :func:`~HS_Image.read_image`

    - ``location``: the path to a folder with image files to be read or the full path
      name of an image file to be read.

      :Type: *string*
      :Required by: :func:`~HS_Image.read_image`

    - ``first_frame``: a positive *integer* specifying the index of the first image
      frame file to be read

      This variable is applicable only when ``location_type`` is 'folder'.

      When ``fprefix`` is also given, ``first_frame`` is the index of the first file
      to be read with that file prefix, from a folder specified by ``location``.

      When ``fprefix`` is not given ``first_frame`` is the index of the first file to
      be read from a folder specified by ``location``.

      :Type: positive *integer*
      :Required by: (optional) :func:`~HS_Image.set_acquisition_parameters`,
                    :func:`~HS_Image.read_and_average_frames`
      :Updated by: (optional) :func:`~HS_Image.set_acquisition_parameters`,
                   :func:`~HS_Image.read_single_frame`

    - ``last_frame``: a positive *integer* specifying the index of the last image
      frame file to be read

      This variable is applicable only when ``location_type`` is 'folder'.

      When ``fprefix`` is also given, ``last_frame`` is the index of the last file
      to be read with that file prefix, from a folder specified by ``location``.

      When ``fprefix`` is not given ``last_frame`` is the index of the last file to
      be read from a folder specified by ``location``.

      :Type: positive *integer*
      :Required by: (optional) :func:`~HS_Image.set_acquisition_parameters`,
                    :func:`~HS_Image.read_and_average_frames`
      :Updated by: (optional) :func:`~HS_Image.set_acquisition_parameters`,
                   :func:`~HS_Image.read_single_frame`

    - ``no_of_frames``: the number of the image frame files to be read.

      :Type: positive *integer*
      :Required by: :func:`~HS_Image.set_acquisition_parameters` (optional),
                    :func:`~HS_Image.read_and_average_frames`
      :Updated by: :func:`~HS_Image.set_acquisition_parameters` (optional),
                   :func:`~HS_Image.read_single_frame`
    
    - ``fprefix``: a file name prefix of the image frame files to be read

      :Type: *string*
      :Required by: (optional) :func:`~HS_Image.set_acquisition_parameters`,
                    (optional) :func:`~HS_Image.read_and_average_frames`

    - ``original_image``: an image read from a file, read and averaged from multiple
      files

      :Type: two-dimensional *ndarray*
      :Required by: :func:`~HS_Image.process_image`,
                    :func:`~HS_Image.locate_bad_pixels`,
                    :func:`~HS_Image.fix_bad_pixels`,
                    :func:`~HS_Image.calculate_sorted_average`,
                    :func:`HS_Centroids.determine_resolution() <HWS.HS_Centroids.HS_Centroids.determine_resolution>`
      :Updated by: :func:`~HS_Image.clear_images`,
                   :func:`~HS_Image.read_single_frame`,
                   :func:`~HS_Image.read_and_average_frames`

    - ``modified_image``: an image obtained by processing ``original_image``

      :Type: *float* or two-dimensional *ndarray* of *floats*
      :Required by: :func:`HS_Centroids.find_centroids_from_image() <HWS.HS_Centroids.HS_Centroids.find_centroids_from_image>`,
                    :func:`HS_Centroids.find_centroids_using_template() <HWS.HS_Centroids.HS_Centroids.find_centroids_using_template>`
      :Updated by: :func:`~HS_Image.process_image`

    - ``background``: The back ground pixel value or an image to be subtracted from
      ``original_image``

      ``background`` is normally obtained from an image acquired without light.

      :Type: non-negative *float* or *ndarray* of non-negative *floats*
      :Required by: :func:`~HS_Image.process_image`

    - ``valid_filesize``: The expected file size of an image file.

      :Type: positive *integer*
      :Required by: :func:`~HS_Image.read_raw`

    - ``bad_pixels``: An array of the x and y coordinates of the pixels tagged to be
      "bad".

      See :func:`~HS_Image.locate_bad_pixels` for an example of how some pixels can
      be tagged to be "bad", and :func:`~HS_Image.fix_bad_pixels` for an example of
      how to deal with them.

      :Type: two-dimensional *ndarray*
      :Requried by: :func:`~HS_Image.fix_bad_pixels`
      :Updated by: :func:`~HS_Image.locate_bad_pixels`

    - ``to_fix_bad_pixels``: a *boolean* to indicate whether to fix bad pixels or not

      If ``True``, :func:`~HS_Image.process_image` will call
      :func:`~HS_Image.fix_bad_pixels` before subtracting the background from
      ``original_image``.

      :Required by: :func:`~HS_Image.process_image`

    - ``threshold_bad_pixels``: a threshold pixel value used to determine which
      pixels are bad

      :func:`~HS_Image.locate_bad_pixels` will tag any pixels whose values are higher
      than ``threshold_bad_pixels`` as bad pixels. The image for this operation must
      be an dark image in most cases.

      If a user specifies the input parameter ``pv_threshold`` when calling
      :func:`~HS_Image.locate_bad_pixels`, the method will use that value and update
      ``threshold_bad_pixels`` as well.

      :Type: no-negative *float*
      :Required by: (optional) :func:`~HS_Image.locate_bad_pixels`
      :Updated by: (optional) :func:`~HS_Image.locate_bad_pixels`

    - ``total_pv``: the sum of all pixel values in ``modified_image``

      :Type: *float*
      :Updated by: :func:`~HS_Image.process_image`

    .. note:: The variables below are currently not in use, and may be removed in
       future revisions unless found useful.

    - ``first_moment``
    - ``second_moment``
    - ``hi_filter_width``
    - ``lo_filter_width``
    - ``hi_filter_power``
    - ``lo_filter_power``
    - ``gpstime``
    - ``max_valid_pv``
    
    """

    def __init__(self, 
                 location_type = None, 
                 location = None,
                 fprefix = None,
                 first_frame = None,
                 last_frame = None,
                 no_of_frames = None,
                 background = 0):

        self.location_type = location_type
        self.location = location
        self.first_frame = first_frame
        self.last_frame = last_frame
        self.no_of_frames = no_of_frames
        self.fprefix = fprefix
        self.original_image = None
        self.modified_image = None
        self.background = None
        self.max_valid_pv = None
        self.valid_filesize = 2097152
        self.bad_pixels = None
        self.to_fix_bad_pixels = None
        self.to_mask_image = None
        self.mask_radius = None
        self.mask_center = None
        self.mask_blur   = None
        self.threshold_bad_pixels = None
        self.first_moment = None
        self.second_moment = None
        self.hi_filter_width = 10.
        self.hi_filter_power = 10.
        self.lo_filter_width = 20.
        self.lo_filter_power = 10.
        self.total_pv = None
        self.gpstime = None

        if self.location_type == 'file':
            self.no_of_frames = 1

        if self.location_type == 'folder':
            self.set_acquisition_parameters()


    def mask_image(self,h=1024,w=1024,center=None,radius=None,blur=None):
        if center is None: # use the middle of the image
            center = [int(w/2), int(h/2)]
        if radius is None: # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w-center[0], h-center[1])

        Y, X = ogrid[:h, :w]
        dist_from_center = sqrt((X - center[0])**2 + (Y-center[1])**2)

        mask = dist_from_center<=radius
        mask_out = mask.astype(float)

        mask_out[mask_out>0] = 1
        
        from scipy.ndimage.filters import gaussian_filter
        
        if blur is not None:
            mask_out = gaussian_filter(mask_out, blur)    
        
        return mask_out


    def set_acquisition_parameters(self):
        """Set the various instance variables in preparation for reading the
        image frames.

        This method assigns values to the instance variables based on a user's
        choice when instantiating the class, and raises an exception when there
        is a inconsistency.

        The method checks the following.

        If ``no_of_frames`` is not set by a user:

        - Both ``first_frame`` and ``last_frame`` must be given as positive *integers*.
        - ``first_frame`` must not be larger than ``last_frame``.
        - Then ``no_of_frames`` is set to ``last_frame`` - ``first_frame`` + 1

        If ``no_of_frames`` is set by a user:

        - ``no_of_frames`` must be a positive *integer*.
        - ``first_frame`` can either be ``None`` or a positive *integer*. If ``None``,
          it is set to 1.
        - ``last_frame`` can either be ``None`` or a positive *integer*. If ``None``, it
          is set to ``first_frame`` + ``no_of_frames`` - 1. If not ``None``, the value
          must be equal to that number.

        :Reruires: (optional) ``first_frame``,
                   (optional) ``last_frame``,
                   (optional) ``no_of_frames``
        :Updates: (optional) ``first_frame``,
                  (optional) ``last_frame``,
                  (optional) ``no_of_frames``

        """
        if self.no_of_frames is None:
            if (self.first_frame is None) and (self.last_frame is None):
                raise \
                  Exception('first_frame and last_frame must be set ' + \
                              'if no_of_frames is not set.')        

            if self.first_frame is None:
                raise \
                  Exception('first_frame must be set if no_of_frames is not set.')

            if self.is_pint(self.first_frame) is False:
                raise \
                  Exception('first_frame must be a positive integer.')

            if self.last_frame is None:
                raise \
                  Exception('last_frame must be set if no_of_frames is not set.')

            if self.is_pint(self.last_frame) is False:
                raise \
                  Exception('last_frame must be a positive integer.')
                
            if self.first_frame > self.last_frame:
                raise \
                    Exception('last_frame must be equal or bigger than first_frame.')
            self.no_of_frames = self.last_frame - self.first_frame + 1

        else:
            if self.is_pint(self.no_of_frames) is False:
                    raise \
                      Exception('no_of_frames must be a positive integer.')
            if self.first_frame is not None:
                if self.is_pint(self.first_frame) is False:
                    raise \
                      Exception('first_frame must be a positive integer.')
            else:
                self.first_frame = 1

            if self.last_frame is not None:
                if self.is_pint(self.first_frame) is False:
                    raise \
                      Exception('first_frame must be a positive integer.')

                if self.no_of_frames != self.last_frame - self.first_frame + 1:
                    raise \
                        Exception('first_frame, last_frame, and no_of_frames ' + \
                                  'are not consistent.')
            else:
                self.last_frame = self.first_frame + self.no_of_frames - 1

    def is_pint(self,n):
        """Check if a parameter is a positive *integer*.

        :param n: the parameter to be checked
        :Returns: ``True`` or ``False``
        :rtype: *boolean*

        """
        is_number = isinstance(n,Number)
        if is_number is True:
            return (n == int(n)) and (n > 0)
        else:
            return False

    def read_image(self):
        """Read an image file, or read multiple image files and average them.

        :Requires: ``location``, ``location_type``

        """
        if self.location == None:
            raise \
              Exception('Instance variable location is not set.')

        if self.location_type == 'file':
            self.read_single_frame()
        elif self.location_type == 'folder':
            self.read_and_average_frames()
        else:
            raise \
              Exception('Instance variable location_type must be file or folder.')

    def clear_images(self):
        """Set ``original_image`` and ``modified_image`` to ``None``.

        :Updates: ``original_image``, ``modified_image``

        """
        self.original_image = None
        self.modified_image = None

    def process_image(self):
        """Process ``original_image`` and set the resulting image to
        ``modified_image``.

        If the instance variable ``to_fix_bad_pixels`` is ``True``, the method calls
        :func:`~HS_Image.fix_bad_pixels` to replace the pixel value of the bad pixels
        with the average of the surrounding pixel values.

        This method then subtracts the background from ``original_image``, and sets any
        pixel values that became negative to zero.

        Finally it sums all the pixel values of ``modified_image`` and sets the value
        to ``total_pv``.

        :Requires: ``background``, ``original_image``
        :Updates: ``modified_image``, ``total_pv``

        """
        if self.original_image is None:
            raise Exception('There is no original_image to process.')

        if self.background is None:
            raise Exception('Instance variable background is not set.')

        tbg = type(self.background)

        if (tbg in [int, float, ndarray]) == False:
            raise \
              Exception('Instance variable background must be a number ' + \
                        'or a numpy array of the same shape as that ' + \
                        'of original_image.')

        if tbg == ndarray:
            if (len(self.background) != 1) and \
               (self.background.shape != self.original_image.shape):
                raise \
                  Exception('The shape of Instance variable background does not ' + \
                            'match the shape of original_image.')

        if self.to_fix_bad_pixels is True:
            fim = self.fix_bad_pixels()
        else:
            fim = self.original_image

        self.modified_image = fim - self.background
        self.modified_image[self.modified_image < 0] = 0

        if self.to_mask_image is True:
            h = self.modified_image.shape[0]
            w = self.modified_image.shape[1]
            mask = self.mask_image(h=h, w=w,
                                   center=array(self.mask_center),
                                   radius=self.mask_radius,
                                   blur=self.mask_blur)
            self.modified_image = multiply(self.modified_image, mask)

        self.total_pv = sum(self.modified_image)

    def read_single_frame(self):
        """Readd a single image frame file whose full path name is gien by the
        instance variable ``location``.

        This method also ensures ``first_frame`` and ``last_frame`` is set to
        ``None``, and ``no_of_frames`` to 1.

        :Requires: ``location``
        :Uses: :func:`~HS_Image.read_raw`
        :Updates: ``first_frame``, ``last_frame``, ``no_of_frames``
                  ``original_image``

        """
        self.first_frame = None
        self.last_frame = None
        self.no_of_frames = 1

        if self.location is None:
            raise \
              Exception('the instance variable location must be the ` + \
                        `full path file name to be read.')

        self.original_image = self.read_raw(self.location)

    def read_and_average_frames(self):
        """Read multiple image frame files and average them to make an iamge.

        :Requires: ``location``, ``first_frame``, ``last_frame``, ``no_of_frames``,
                   (optional) ``fprefix``
        :Uses: :func:`~HS_Image.read_raw`
        :Updates: ``no_of_frames``, ``original_image``

        """
        if self.fprefix is None:
            files = [os.path.basename(f) \
                     for f in glob.glob(self.location + '*.raw')]
            fpnames = [self.location + f for f in files]
        else:
            files = [os.path.basename(f) \
                     for f in glob.glob(self.location + self.fprefix + '*.raw')]

        if len(files) == 0:
            raise Exception('There are no raw files to read.')

        fpnames = [self.location + f for f in files]

        # Since it is possible to bypass checking mechanism in
        # set_acquisition_parameters, I do provide some flexibility here that
        # the method will still work as long as first_frame and last_frame are
        # given. However, they must still be positive integers.
        if self.is_pint(self.first_frame) is False:
            raise \
              Exception('first_frame must be a positive integer.')

        if self.is_pint(self.last_frame) is False:
            raise \
              Exception('last_frame must be a positive integer.')

        self.no_of_frames = abs(self.first_frame - self.last_frame) + 1
        self.original_image = zeros(shape(self.read_raw(fpnames[0])));

        for fpn in fpnames[self.first_frame-1:self.last_frame]:
            self.original_image += self.read_raw(fpn)

        self.original_image = self.original_image/float(self.no_of_frames)

    def read_and_average_frames_old(self):
        """An old version of :func:`~HS_Image.read_and_average_frames`

        It is not recommended to use this method, but is here for debugging purposes
        in case something goes wrong.

        May be removed in a future revision.

        """
        if self.fprefix is None:
            files = [os.path.basename(f) \
                     for f in glob.glob(self.location + '*.raw')]
            fpnames = [self.location + f for f in files]
        else:
            files = [os.path.basename(f) \
                     for f in glob.glob(self.location + self.fprefix + '*.raw')]

        if len(files) == 0:
            raise Exception('There are no raw files to read.')

        fpnames = [self.location + f for f in files]

        # If first_frame is not defined, assume it is 1
        if self.first_frame is None:
            self.first_frame = 1
        elif type(self.first_frame) != int:
            raise \
              Exception('Instance variable first_frame must be an integer.')

        # If last_frame is not defined, assume it is the last file
        if self.last_frame is None:
            self.last_frame = len(fpnames)
            self.no_of_frames = self.last_frame
        elif type(self.last_frame) != int:
            raise \
              Exception('Instance variable last_frame must be an integer.')

        if (self.first_frame > len(fpnames)):
            raise \
              Exception('first_frame ' + str(self.first_frame) + ' is too big.' + \
                        ' There are only ' + str(len(fpnames)) + ' files avialable.')

        if (self.last_frame > len(fpnames)):
            raise \
              Exception('last_frame ' + str(self.last_frame) + ' is too big.' + \
                        ' There are only ' + str(len(fpnames)) + ' files avialable.')

        self.no_of_frames = abs(self.first_frame - self.last_frame) + 1
        self.original_image = zeros(shape(self.read_raw(fpnames[1])));

        for fpn in fpnames[self.first_frame-1:self.last_frame]:
            self.original_image += self.read_raw(fpn)

        self.original_image = self.original_image/self.no_of_frames

    def read_raw(self,fname):
        """Read a single image frame file.

        This method is device dependent, and is supposed to work with Dalsa
        Pantera 1M60 cameras. 

        :param fname: the full path name of the file to be read
        :type fname: *string*
        :Requires: ``valid_filesize``
        :Returns: an image read from the file
        :rtype: two-dimensional *ndarray*

        """
        fsize = int(os.stat(fname).st_size)

        if os.stat(fname).st_size != self.valid_filesize:
            raise Exception('The file size of ' + fname + ' is not valid:' \
                            + '\n The size should be ' \
                            + str(self.valid_filesize) \
                            + ' but was ' + str(fsize) + '.')
        # Below is taken from python tutorial:
        # https://docs.python.org/2/library/os.html
        try:
            fp = open(fname)
        except IOError as e:
            print "I/O error({0}): {1}".format(e.errno, e.strerror)
            raise
        else:
            with fp:
                ri = fromfile(fp,uint16).reshape(1024,1024).astype(float)
                return ri

    def locate_bad_pixels(self,pv_threshold):
        """Locate the pixels whose values exceed ``pv_threshold``

        :param pv_threshold: threshold pixel value
        :type pv_threshold: positive *float*
        :Requires: ``original_image``
        :Updates: ``threshold_bad_pixels``, ``bad_pixels``

        """
        self.bad_pixels = array(where(self.original_image > pv_threshold)).T
        # Convert to xy format from row column by swapping columns
        self.bad_pixels[:,[0,1]] = self.bad_pixels[:,[1,0]]
        self.threshold_bad_pixels = pv_threshold

    def fix_bad_pixels(self):
        """Replace the values of the bad pixels with the averaged values of their
        surrounding pixels.

        For each pixel from the instance variable ``bad_pixels``, this method
        calculates the average of the pixel values that surround it, and sets that
        average value as the value of the pixel.

        :Requires: ``original_image``, ``bad_pixels``
        :Returns: an image with bad pixels fixed
        :rtype: two-dimensional *ndarray* of the shape same as that of
                ``original_image``

        """
        fim = self.original_image.copy()
        if self.bad_pixels is None:
            raise Exception('Instance variable bad_pixels is not defined:' \
                            + 'There are no pixels to fix.')

        np = shape(self.bad_pixels)[0]
        bp = self.bad_pixels.copy()

        # Convert xy to row column format
        bp[:,[0,1]] = bp[:,[1,0]]
        rmax,cmax = shape(fim)
        for k in arange(np):
            lbr = max(0,bp[k,0]-1)
            ubr = min(rmax,bp[k,0]+1)
            lbc = max(0,bp[k,1]-1)
            ubc = min(cmax,bp[k,1]+1)
            fpc = 0
            fim[tuple(bp[k])] = 0
            nop = size(fim)-1
            pv_mean = fim[lbr:ubr+1,lbc:ubc+1].sum()/float(nop)
            fim[tuple(bp[k])] = pv_mean

        return fim

    def calculate_sorted_average(self,no_of_pixels=100,image=None):
        """Sort the pixel values, take a number of values from the highest, and
        calculate their average.

        :param no_of_pixels: the number of pixels to take after sorting
        :type no_of_pixels: positive *integer*
        :param image: the image whose pixel values are to be sorted
        :type image: *ndarray*
        :Requires: (optional) ``original_image``
        :Returns: an average of the sorted pixel values
        :rtype: *float*

        """
        if image is None:
            im = self.original_image.ravel()
        else:
            im = image.ravel()
        
        idx = im.argsort()[::-1][:no_of_pixels]
        return mean(im[idx])
