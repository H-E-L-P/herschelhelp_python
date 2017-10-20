# -*- coding: utf-8 -*-

import logging
import os

import pyvo as vo
from pymoc import MOC

logging.basicConfig(level=logging.getLevelName(
    os.getenv('LOG_LEVEL', 'WARNING')))
LOGGER = logging.getLogger(__name__)


def get_depth_coverage(coverage_dict):
    """Generate a Multi-Order Coverage maps corresponding to depth constraints

    This method connects to HELP Virtual Observatory server to get the list of
    HEALPix cells corresponding to depth constraints. The constraints are
    expressed as a dictionnary associated to a band name the minimum depth
    value in this band.

    This method returns the Multi-Order Coverage maps of the corresponding
    area.

    Parameters
    ----------
    coverage_dict : dict
        Dictionnary associating to a band name the required minimum depth
        value.

    Returns
    -------
    pymoc.MOC
        The Multi-Order Coverage map corresponding to the depth constraints in
        HELP data.

    """

    moc_list = []

    for band in coverage_dict:

        query = "select top 100000000 healpix_order, healpix_npix from " \
                "depth.main where band = '{}' and depth >= {}".format(
                    str(band), float(coverage_dict[band]))

        LOGGER.debug("Query sent to VO server: %s", query)
        vo_result = vo.tablesearch(
            "http://vohedamtest.lam.fr/__system__/tap/run/tap", query).table

        moc = MOC()
        if len(vo_result):
            for group in vo_result.group_by('healpix_order').groups:
                order = group['healpix_order'][0]
                cells = group['healpix_npix']
                moc.add(order, cells)

        moc_list.append(moc)

    result_moc = moc_list.pop()
    for moc in moc_list:
        result_moc = result_moc.intersection(moc)

    return result_moc
