# -*- coding: utf-8 -*-

import os

import click
from astropy.table import Table
from pymoc import MOC

from .depth_coverage import get_depth_coverage
from .footprints import compute_coverages


@click.group()
def cli():
    """Herschel Extra-galactic Legacy Programme tool.

    Execute "herschelhelp <COMMAND> --help" to display the help of a command.
    """
    pass


@cli.command(short_help="Prints % of HELP covered by a MOC")
@click.argument("moc_file", metavar="<moc_file>")
def covbymoc(moc_file):
    """Prints the percentage of HELP covered by a MOC.

    Given a FITS Multi-order coverage map (MOC) footprint, this command
    displays the percentage of each HELP field covered by it.

    <moc_file> must be a FITS file containing a MOC.

    """

    footprint = MOC()
    footprint.read(moc_file)

    coverages = compute_coverages(footprint)

    output = Table(
        data=[list(coverages.keys()), list(coverages.values())],
        names=["Field", "Coverage"]
    )

    print(output)


@cli.command(short_help="Gets the MOC corresponding to depth criteria")
@click.option('-k', default=None, type=float, help="K limiting magnitude.")
@click.option('--output', default="coverage.fits",
              help="Name of the output FITS file containing the MOC "\
              "(coverage.fits by default).")
def mocfromdepth(k, output):
    """Produce the MOC corresponding to depth constraints.

    Given some constraints on the depth of the data (see options below) this
    command produce a Multi-Order Coverage map (MOC) of the corresponding area.
    The programme connects to the HELP Virtual Observatory server to get the
    information.

    """

    if os.path.isfile(output):
        raise click.BadOptionUsage("The file {} already exists".format(output))

    coverage_dict = {}

    if k:
        coverage_dict['K'] = k

    if len(coverage_dict) == 0:
        raise click.BadParameter("You must provide at least one depth limit. "
                                 "See herschelhelp mocfromdepth --help.")
    else:
        result = get_depth_coverage(coverage_dict)
        result.write(output)
        print("The corresponding area was saved to {}. It covers {} square " \
              "degrees.".format(output, result.area_sq_deg))
