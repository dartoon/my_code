#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 11:57:05 2018

@author: dartoon
"""

def extract_stars(data, catalogs, size=(11, 11)):
    """
    Extract cutout images centered on stars defined in the input
    catalog(s).

    Stars where the cutout array bounds partially or completely lie
    outside of the input ``data`` image will not be extracted.

    Parameters
    ----------
    data : `~astropy.nddata.NDData` or list of `~astropy.nddata.NDData`
        A `~astropy.nddata.NDData` object or a list of
        `~astropy.nddata.NDData` objects containing the 2D image(s) from
        which to extract the stars.  If the input ``catalogs`` contain
        only the sky coordinates (i.e. not the pixel coordinates) of the
        stars then each of the `~astropy.nddata.NDData` objects must
        have a valid ``wcs`` attribute.

    catalogs : `~astropy.table.Table`, list of `~astropy.table.Table`
        A catalog or list of catalogs of sources to be extracted from
        the input ``data``.  To link stars in multiple images as a
        single source, you must use a single source catalog where the
        positions defined in sky coordinates.

        If a list of catalogs is input (or a single catalog with a
        single `~astropy.nddata.NDData` object), they are assumed to
        correspond to the list of `~astropy.nddata.NDData` objects input
        in ``data`` (i.e. a separate source catalog for each 2D image).
        For this case, the center of each source can be defined either
        in pixel coordinates (in ``x`` and ``y`` columns) or sky
        coordinates (in a ``skycoord`` column containing a
        `~astropy.coordinates.SkyCoord` object).  If both are specified,
        then the pixel coordinates will be used.

        If a single source catalog is input with multiple
        `~astropy.nddata.NDData` objects, then these sources will be
        extracted from every 2D image in the input ``data``.  In this
        case, the sky coordinates for each source must be specified as a
        `~astropy.coordinates.SkyCoord` object contained in a column
        called ``skycoord``.  Each `~astropy.nddata.NDData` object in
        the input ``data`` must also have a valid ``wcs`` attribute.
        Pixel coordinates (in ``x`` and ``y`` columns) will be ignored.

        Optionally, each catalog may also contain an ``id`` column
        representing the ID/name of stars.  If this column is not
        present then the extracted stars will be given an ``id`` number
        corresponding the the table row number (starting at 1).  Any
        other columns present in the input ``catalogs`` will be ignored.

    size : int or array_like (int), optional
        The extraction box size along each axis.  If ``size`` is a
        scalar then a square box of size ``size`` will be used.  If
        ``size`` has two elements, they should be in ``(ny, nx)`` order.
        The size must be greater than or equal to 3 pixel for both axes.

    Returns
    -------
    stars : `Stars` instance
        A `Stars` instance containing the extracted stars.
    """

    if isinstance(data, NDData):
        data = [data]

    if isinstance(catalogs, Table):
        catalogs = [catalogs]

    for img in data:
        if not isinstance(img, NDData):
            raise ValueError('data must be a single or list of NDData '
                             'objects.')

    for cat in catalogs:
        if not isinstance(cat, Table):
            raise ValueError('catalogs must be a single or list of Table '
                             'objects.')

    if len(catalogs) == 1 and len(data) > 1:
        if 'skycoord' not in catalogs[0].colnames:
            raise ValueError('When inputting a single catalog with multiple '
                             'NDData objects, the catalog must have a '
                             '"skycoord" column.')

        if any([img.wcs is None for img in data]):
            raise ValueError('When inputting a single catalog with multiple '
                             'NDData objects, each NDData object must have '
                             'a wcs attribute.')
    else:
        for cat in catalogs:
            if 'x' not in cat.colnames or 'y' not in cat.colnames:
                if 'skycoord' not in cat.colnames:
                    raise ValueError('When inputting multiple catalogs, '
                                     'each one must have a "x" and "y" '
                                     'column or a "skycoord" column.')
                else:
                    if any([img.wcs is None for img in data]):
                        raise ValueError('When inputting catalog(s) with '
                                         'only skycoord positions, each '
                                         'NDData object must have a wcs '
                                         'attribute.')

        if len(data) != len(catalogs):
            raise ValueError('When inputting multiple catalogs, the number '
                             'of catalogs must match the number of input '
                             'images.')

    size = np.atleast_1d(size)
    if len(size) == 1:
        size = np.repeat(size, 2)

    min_size = 3
    if size[0] < min_size or size[1] < min_size:
        raise ValueError('size must be >= {} for x and y'.format(min_size))

    if len(catalogs) == 1:    # may included linked stars
        use_xy = True
        if len(data) > 1:
            use_xy = False    # linked stars require skycoord positions

        stars = []
        # stars is a list of lists, one list of stars in each image
        for img in data:
            stars.append(_extract_stars(img, catalogs[0], size=size,
                                        use_xy=use_xy))

        # transpose the list of lists, to associate linked stars
        stars = list(map(list, zip(*stars)))

        # remove 'None' stars (i.e. no or partial overlap in one or more
        # images) and handle the case of only one "linked" star
        stars_out = []
        n_input = len(catalogs[0]) * len(data)
        n_extracted = 0
        for star in stars:
            good_stars = [i for i in star if i is not None]
            n_extracted += len(good_stars)
            if len(good_stars) == 0:
                continue    # no overlap in any image
            elif len(good_stars) == 1:
                good_stars = good_stars[0]  # only one star, cannot be linked
            else:
                good_stars = LinkedStar(good_stars)

            stars_out.append(good_stars)
    else:    # no linked stars
        stars_out = []
        for img, cat in zip(data, catalogs):
            stars_out.extend(_extract_stars(img, cat, size=size, use_xy=True))

        n_input = len(stars_out)
        stars_out = [star for star in stars_out if star is not None]
        n_extracted = len(stars_out)

    n_excluded = n_input - n_extracted
    if n_excluded > 0:
        warnings.warn('{} star(s) were not extracted because their cutout '
                      'region extended beyond the input image.'
                      .format(n_excluded), AstropyUserWarning)

    return Stars(stars_out)