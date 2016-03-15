# -*- coding: utf-8 -*-

import click
from astropy.table import Table
from pymoc import MOC

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

