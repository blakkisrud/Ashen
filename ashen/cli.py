"""

Command line interface for the ashen application.

"""

import click
import csv
from ashen.ashen_utils import (
    make_decay_chain_db,
    load_icrp_107,
    energy_in_decay_chain,
    get_daughters,
)

@click.group()
def cli():
    """Ashen: Extract radionuclide data from ICRP 107."""
    pass

@cli.command()
@click.argument('nuclide_name')
def get_half_life(nuclide_name):
    """Get the half-life of a radionuclide."""
    # Load decay chains and emission data
    emission_energy = load_icrp_107()
    db = make_decay_chain_db(emission_data=emission_energy)
    # Get the half-life of the specified nuclide
    nuclide = db.get_decay_info(nuclide_name)
    result = nuclide.halflife

    click.echo(f"The half-life of {nuclide_name} is {result}.")

@cli.command()
@click.argument('nuclide_name')
@click.option('--energy-type', default="n", help='Type of energy, (a)lpha, (b)eta, (n)on-penetrative.')
@click.option('--verbose', is_flag=True, help='Enable verbose output.')

def energy(nuclide_name, energy_type, verbose):
    # Load decay chains and emission data
    emission_energy = load_icrp_107()
    db = make_decay_chain_db(emission_data=emission_energy)
    result = energy_in_decay_chain(
        db=db,
        nuc_name=nuclide_name,
    )

    if energy_type == "a":
        if verbose:
            click.echo(f"Verbose mode enabled. Showing detailed energy information for {nuclide_name}.")
            click.echo(f"Total Alpha Energy: {result['Total Alpha Energy']} MeV")
        else:
            click.echo(f"{result['Total Alpha Energy']}")
    elif energy_type == "b":
        if verbose:
            click.echo(f"Verbose mode enabled. Showing detailed energy information for {nuclide_name}.")
            click.echo(f"Total Beta Energy: {result['Total Beta Energy']} MeV")
        else:
            click.echo(f"{result['Total Beta Energy']}")
    elif energy_type == "n":
        if verbose:
            click.echo(f"Verbose mode enabled. Showing detailed energy information for {nuclide_name}.")
            click.echo(f"Total Non-Penetrative Energy: {result['Total Non-Penetrative Energy']} MeV")
        else:
            click.echo(f"{result['Total Non-Penetrative Energy']}")
    else:
        click.echo("Invalid type. Use 'a' for alpha, 'b' for beta, or 'n' for non-penetrative.")

@cli.command(name="batch-energy")
@click.argument('input_csv', type=click.Path(exists=True))
@click.option('--energy-type', default="n", help="Type of energy: (a)lpha, (b)eta, (n)on-penetrative.")
@click.option('--verbose', is_flag=True, help='Enable verbose output.')
def batch_energy(input_csv, energy_type, verbose):
    """
    Calculate energy for multiple nuclides listed in a CSV file.
    CSV should have one nuclide name per line (no header required).
    """
    # Load the database once
    emission_energy = load_icrp_107()
    db = make_decay_chain_db(emission_data=emission_energy)

    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if not row:
                continue  # skip empty lines
            nuclide_name = row[0].strip()

            try:
                result = energy_in_decay_chain(db=db, nuc_name=nuclide_name)
            except KeyError:
                click.echo(f"Nuclide {nuclide_name} not found in database. Skipping.")
                continue

            if verbose:
                click.echo(f"Results for {nuclide_name}:")
            if energy_type == "a":
                if verbose:
                    click.echo(f"  Verbose mode enabled. Showing detailed energy information for {nuclide_name}.")
                    click.echo(f"  Total Alpha Energy: {result['Total Alpha Energy']}")
                else:
                    click.echo(f"{result['Total Alpha Energy']}")
            elif energy_type == "b":
                if verbose:
                    click.echo(f"  Verbose mode enabled. Showing detailed energy information for {nuclide_name}.")
                    click.echo(f"  Total Beta Energy: {result['Total Beta Energy']}")
                else:
                    click.echo(f"{result['Total Beta Energy']}")
            elif energy_type == "n":
                if verbose:
                    click.echo(f"  Verbose mode enabled. Showing detailed energy information for {nuclide_name}.")
                    click.echo(f"  Total Non-Penetrative Energy: {result['Total Non-Penetrative Energy']}")
                else:
                    click.echo(f"{result['Total Non-Penetrative Energy']}")
            else:
                click.echo("Invalid energy type specified.")

if __name__ == "__main__":
    cli()


