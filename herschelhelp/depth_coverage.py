# -*- coding: utf-8 -*-

import logging
import os

import pyvo as vo
from pymoc import MOC

logging.basicConfig(level=logging.getLevelName(
    os.getenv('LOG_LEVEL', 'WARNING')))
LOGGER = logging.getLogger(__name__)


def get_depth_coverage(coverage_limits):
    """Generate a Multi-Order Coverage maps corresponding to depth constraints
    This method connects to HELP Virtual Observatory server to get the list of
    HEALPix cells corresponding to depth constraints. The constraints are
    expressed as a dictionnary associating to a band name the maximum error
    value in this band.
    This method returns the Multi-Order Coverage maps of the corresponding
    area.
    Parameters
    ----------
    coverage_limits : dict
        Dictionary associating to a band name the required maximum error
        in this band flux, in μJy.
    Returns
    -------
    pymoc.MOC
        The Multi-Order Coverage map corresponding to the depth constraints in
        HELP data.
    """

    where_clause = ' and '.join([
        f"ferr_{band}_mean < {error} " for band, error in
        coverage_limits.items()
    ])

    query = ("select top 100000000 hp_idx_O_13 from  depth.main where " +
             where_clause)
    print(query)

    vo_result = vo.tablesearch(
        "http://vohedamtest.lam.fr/__system__/tap/run/tap", query).table

    moc = MOC()
    if len(vo_result):
        moc.add(13, vo_result['hp_idx_o_13'])

    return moc
