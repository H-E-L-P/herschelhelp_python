# -*- coding: utf-8 -*-

import logging
import os
from collections import OrderedDict

from .database import get_field

logging.basicConfig(level=logging.getLevelName(os.getenv('LOG_LEVEL', 'INFO')))
LOGGER = logging.getLogger(__name__)


def compute_coverages(footprint):
    """Compute the coverage of a footprint over HELP fields.

    Given a Multi-Order Coverage map, this function compares the footprint to
    HELP coverage giving the coverage percentage in each field.

    Parameters
    ----------
    footprint : pymoc.MOC
        Compared footprint.

    Returns
    -------
    OrderedDict
        The keys of the dictionary are the field names and the values the
        percentage of the field covered.

    """

    result = OrderedDict()

    for field in get_field():
        LOGGER.debug("Computing %s coverage...", field.name)
        intersection_area = (
            100. * footprint.intersection(field.footprint).area /
            field.footprint.area)
        if intersection_area < 0.1:
            result[field.name] = round(intersection_area, ndigits=2)
        else:
            result[field.name] = round(intersection_area, ndigits=1)

    return result
