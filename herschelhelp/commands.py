# -*- coding: utf-8 -*-

import os

import click
from astropy.io.registry import IORegistryError
from astropy.table import Table
from pymoc import MOC

from .depth_coverage import get_depth_coverage
from .filters import (correct_galactic_extinction, export_to_cigale,
                      export_to_eazy, get_filter_meta_table)
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
@click.option('--psw', default=None, type=float, help="SPIRE250 5 sigma "\
              "limiting magnitude.")
@click.option('--output', default="coverage.fits",
              help="Name of the output FITS file containing the MOC "\
              "(coverage.fits by default).")
def mocfromdepth(k, psw, output):
    """Produce the MOC corresponding to depth constraints.

    Given some constraints on the depth of the data (see options below) this
    command produce a Multi-Order Coverage map (MOC) of the corresponding area.
    The programme connects to the HELP Virtual Observatory server to get the
    information.

    """

    if os.path.isfile(output):
        raise click.BadOptionUsage("The file {} already exists".format(output))

    coverage_dict = {}

    if k is not None:
        coverage_dict['K'] = k
    if psw is not None:
        coverage_dict['PSW'] = psw

    if len(coverage_dict) == 0:
        raise click.BadParameter("You must provide at least one depth limit. "
                                 "See herschelhelp mocfromdepth --help.")
    else:
        result = get_depth_coverage(coverage_dict)
        result.write(output)
        print("The corresponding area was saved to {}. It covers {} square " \
              "degrees.".format(output, result.area_sq_deg))


@cli.command(short_help="Generate CSV file with filter information.")
@click.option('-f', default="help_filters.csv", type=str,
              help="Name of the CSV file (help_filters.csv by default).")
def filter_meta(f):
    """Save filter information to a CSV file.

    This command creates a CSV file containing the information on the various
    filters used in HELP.

    """
    if os.path.isfile(f):
        raise click.BadOptionUsage("The file {} already exists.".format(f))

    filter_table = get_filter_meta_table()
    filter_table.sort('mean_wavelength')
    filter_table.write(f, format='ascii.csv')


@cli.command(short_help="Generate filter files for CIGALE.")
@click.option('-d', default="filters", type=str,
              help="Directory where to store the filter files."
                   " Must not exist.")
def filter_export_cigale(d):
    """Export HELP filters to be used with CIGALE.

    This command export the filters used within HELP to files that can be
    imported in CIGALE for SED fitting.

    """
    if os.path.isdir(d):
        raise click.BadOptionUsage(
            "The directory {} already exists.".format(d))
    if os.path.isfile(d):
        raise click.BadOptionUsage("A file named {} exists".format(d))

    os.mkdir(d)
    export_to_cigale(d)


@cli.command(short_help="Generate filter files for EAZY.")
@click.argument("filename", metavar="<filename>")
def filter_export_eazy(filename):
    """Export HELP filters to be used with EAZY.

    This command export the HELP filters used within HELP to files that can be
    used with EAZY for photometric redshift computation.  It's argument is the
    name of the catalogue file that will be processed by EASY.  It must be in
    the HELP format with the aperture flux and error columns as f_ap_<band>
    and ferr_ap_<band>.

    Three new files will be created based on the filename of this table:
    - <filename>.res containing the responses of the used filters;
    - <filename>.translate associating the table columns to their filter;
    - <filename>.info containing additional information on the filters.

    """
    if os.path.exists("{}.res".format(filename)):
        raise click.BadArgumentUsage(
            "The file {}.res already exists".format(filename))
    if os.path.exists("{}.translate".format(filename)):
        raise click.BadArgumentUsage(
            "The file {}.translate already exists".format(filename))
    if os.path.exists("{}.info".format(filename)):
        raise click.BadArgumentUsage(
            "The file {}.info already exists".format(filename))

    export_to_eazy(filename)


@cli.command(short_help="Correct a catalogue for galactic extinction.")
@click.argument("filename", metavar="<filename>")
def correct_for_extinction(filename):
    """Correct HELP catalogue for galactic extinction.

    This command takes a HELP formatted catalogue and correct its fluxes and
    magnitudes for galactic extinction using the E(B-V) information associated
    with each source.  The catalogue must have an ebv column and the fluxes and
    magnitudes columns are identified using the HELP column format (f_<filter>,
    ferr_filter, f_ap_<filter>, ...).

    Column associated with filters that are not in the HELP database will not
    be corrected and an error message will be logged with the name of the
    missing filters.

    The corrected catalogue will be save as <filename>_corrected.fits.

    """
    new_name = "{}_corrected.fits".format(os.path.splitext(filename)[0])

    if os.path.exists(new_name):
        raise click.UsageError("{} already exists.".format(new_name))

    # TODO: Find a way to deal with not readable ascii files
    try:
        catalogue = Table.read(filename)
    except IORegistryError:
        catalogue = Table.read(filename, format='ascii')

    catalogue = correct_galactic_extinction(catalogue, inplace=True)
    catalogue.write(new_name)
