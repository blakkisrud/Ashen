"""

Main source file for shared functions

"""

"""

Support functions to retreive the half life from a file from ICRP 107

Script can also contain the loading of different radinuclide data and 
characteristics

"""

import pandas as pd
import numpy as np
import os
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
import sys
import re
from collections import defaultdict


@dataclass
class RadiationEmission:
    energy: float  # in MeV
    yield_fraction: float  # per decay
    radiation_type: str  # 'X', 'α', 'β', 'γ', etc.

# Define a dataclass for the radionuclide data


@dataclass
class Nuclide:
    name: str
    # Keeping as a string to preserve units (e.g., "8.275h", "21.1m", "5.591d")
    halflife: str
    daughters: List[Tuple[str, float]] = field(
        default_factory=list)  # (Daughter name, branching fraction)
    emissions: List[RadiationEmission] = field(default_factory=list)

    def calculate_alpha_energy(self) -> float:
        """Calculate the total alpha energy for the decay of this radionuclide."""
        total_alpha_energy = 0.0
        for emission in self.emissions:
            if emission.radiation_type == "A":
                total_alpha_energy += emission.energy * emission.yield_fraction
        return total_alpha_energy

    def calculate_beta_energy(self) -> float:
        """Calculate the total beta energy for the decay of this radionuclide."""
        total_beta_energy = 0.0
        for emission in self.emissions:
            if emission.radiation_type == "B-":
                total_beta_energy += emission.energy * emission.yield_fraction
        return total_beta_energy

    def calculate_recoil_energy(self) -> float:
        """Calculate the total recoil energy for the decay of this radionuclide."""
        total_recoil_energy = 0.0
        for emission in self.emissions:
            if emission.radiation_type == "AR":
                total_recoil_energy += emission.energy * emission.yield_fraction
        return total_recoil_energy

    def calculate_ic_electrons(self) -> float:
        """Calculate the total internal conversion electrons for the decay of this radionuclide."""
        total_ic_electrons = 0.0
        for emission in self.emissions:
            if emission.radiation_type == "IE":
                total_ic_electrons += emission.energy * emission.yield_fraction
        return total_ic_electrons

    def calculate_auger_electrons(self) -> float:
        """Calculate the total Auger electrons for the decay of this radionuclide."""
        total_auger_electrons = 0.0
        for emission in self.emissions:
            if emission.radiation_type == "AE":
                total_auger_electrons += emission.energy * emission.yield_fraction
        return total_auger_electrons

    def calculate_low_energy_photons(self, cut_off=1e-3) -> float:
        """Calculate the total low-energy photons for the decay of this radionuclide."""
        total_low_energy_photons = 0.0
        for emission in self.emissions:
            if emission.radiation_type == "X" or emission.radiation_type == "G":
                if emission.energy <= cut_off:
                    total_low_energy_photons += emission.energy * emission.yield_fraction
        return total_low_energy_photons

    def calculate_total_non_penetrative_energy(self, cut_off=1e-3) -> float:
        """Calculate the total non-penetrative energy for the decay of this radionuclide."""
        return (self.calculate_alpha_energy() + self.calculate_ic_electrons() +
                self.calculate_auger_electrons() + self.calculate_low_energy_photons(cut_off) + self.calculate_beta_energy())
    
    def calculate_decay_constant(self) -> float:
        """Calculate the decay constant for this radionuclide."""

        half_life_hours = convert_to_hours(float(self.halflife[:-1]), self.halflife[-1])
        return np.log(2) / half_life_hours


class DecayDatabase:
    def __init__(self):
        self.nuclides: Dict[str, Nuclide] = {}

    def add_nuclide(self, nuclide: Nuclide):
        """Adds an already created Nuclide object to the database."""
        self.nuclides[nuclide.name] = nuclide

    def get_decay_info(self, nuclide_name: str):
        """Retrieves decay information for a given nuclide."""
        return self.nuclides.get(nuclide_name, None)

    def display_decay_chain(self, nuclide_name: str, level=0):
        """Recursively prints the decay chain from the given nuclide."""
        nuclide = self.get_decay_info(nuclide_name)
        if not nuclide:
            print("Nuclide not found.")
            return

        indent = "  " * level
        print(f"{indent}{nuclide.name} (Half-life: {nuclide.halflife})")
        for daughter, fraction in nuclide.daughters:
            print(f"{indent}  ├── {daughter} (Branching: {fraction:.3e})")
            self.display_decay_chain(daughter, level + 1)


def convert_to_hours(num, unit):
    if unit == "s":
        return num/3600
    elif unit == "m":
        return num/60
    elif unit == "h":
        return num
    elif unit == "d":
        return num*24
    elif unit == "y":
        return num*24*365
    elif unit == "us":
        return num/3600/1e6
    elif unit == "ns":
        return num/3600/1e9
    elif unit == "ms":
        return num/3600/1e3
    else:
        return None

def decay_part_process(split_line):

    if len(split_line) == 5:
        mother_name = split_line[1]
        half_life = split_line[2]
        branching_fraction = float(split_line[3])
        daughter_name = split_line[4]

        nuclide = Nuclide(mother_name, half_life, [
                          (daughter_name, branching_fraction)])

    elif len(split_line) == 7:

        mother_name = split_line[1]
        half_life = split_line[2]
        branching_fraction_1 = float(split_line[3])
        daughter_name_1 = split_line[4]
        branching_fraction_2 = float(split_line[5])
        daughter_name_2 = split_line[6]

        nuclide = Nuclide(mother_name, half_life, [(daughter_name_1, branching_fraction_1),
                                                   (daughter_name_2, branching_fraction_2)])

    elif len(split_line) == 9:

        mother_name = split_line[1]
        half_life = split_line[2]
        branching_fraction_1 = float(split_line[3])
        daughter_name_1 = split_line[4]
        branching_fraction_2 = float(split_line[5])
        daughter_name_2 = split_line[6]
        branching_fraction_3 = float(split_line[7])
        daughter_name_3 = split_line[8]

        nuclide = Nuclide(mother_name, half_life, [(daughter_name_1, branching_fraction_1),
                                                   (daughter_name_2,
                                                    branching_fraction_2),
                                                   (daughter_name_3, branching_fraction_3)])

    elif len(split_line) == 11:

        mother_name = split_line[1]
        half_life = split_line[2]
        branching_fraction_1 = split_line[3]
        daughter_name_1 = split_line[4]
        branching_fraction_2 = split_line[5]
        daughter_name_2 = split_line[6]
        branching_fraction_3 = split_line[7]
        daughter_name_3 = split_line[8]
        branching_fraction_4 = split_line[9]
        daughter_name_4 = split_line[10]

        nuclide = Nuclide(mother_name, half_life, [(daughter_name_1, branching_fraction_1),
                                                   (daughter_name_2,
                                                    branching_fraction_2),
                                                   (daughter_name_3,
                                                    branching_fraction_3),
                                                   (daughter_name_4, branching_fraction_4)])

    return nuclide

def get_daughters(db: DecayDatabase, nuclide_name: str, parent_fraction=1.0, daughters=None):
    """ 
    Recursively collects all daughters from a given nuclide and their effective branching fractions.

    Args:
        db (DecayDatabase): The database containing nuclide decay data.
        nuclide_name (str): The starting radionuclide.
        parent_fraction (float): The cumulative fraction from the parent (default: 1.0).
        daughters (dict): Dictionary to store results (default: None, initializes internally).

    Returns:
        dict: A dictionary where keys are daughter nuclide names and values are their effective branching ratios.
    """
    if daughters is None:
        daughters = {}

    nuclide = db.get_decay_info(nuclide_name)

    if not nuclide:
        if nuclide_name[-1] != "$":
            print(
                f"Nuclide '{nuclide_name}' not found in database! (Possibly a terminal product)")
        return daughters

    for daughter, fraction in nuclide.daughters:
        # Compute the effective branching fraction
        effective_fraction = parent_fraction * float(fraction)

        if daughter in daughters:
            # Sum contributions from different decay paths
            daughters[daughter] += effective_fraction
        else:
            daughters[daughter] = effective_fraction

        # Recursively find daughters of this daughter nuclide
        get_daughters(db, daughter, effective_fraction, daughters)

    return daughters

def half_life_from_icrp107(path="resources/FULL_RAD_LIST.RAD"):

    halflife_data = {}

    with open(path, 'r') as file:
        for line in file:
            line = ' '.join(line.split())
            split_line = line.split(" ")

            if len(split_line) == 3:
                nuc_name = split_line[0]
                half_life = split_line[1]
                halflife_data[nuc_name] = half_life

    return halflife_data

def make_decay_chain_db(path="resources/DECAY_CHAINS.TXT", emission_data=None, add_single_nuclides=True):

    # Retreive the decay chains - potentially finniky

    with open(path, 'r') as file:
        lines = file.readlines()

    # Find the indexes of the lines that contain "Daughter Products"

    indexes = [i for i, line in enumerate(
        lines) if "Daughter Products" in line]

    decay_chains = []

    for i, index_val in enumerate(indexes):

        if i == len(indexes)-1:
            break

        start_index = index_val-1
        end_index = indexes[i+1]-1

        chunk = lines[start_index:end_index]

        decay_chains.append(chunk)

    db = DecayDatabase()

    for i in range(len(decay_chains)):

        single_chain = decay_chains[i]
        decay_part = single_chain[3:-1]

        for line in decay_part:

            line = ' '.join(line.split())
            split_line = line.split(" ")
            nuc = decay_part_process(split_line=split_line)
            # Here add the emission data
            if emission_data is None:
                pass
            else:
                nuc_emissions = emission_data.get(nuc.name, [])
                nuc.emissions = [RadiationEmission(
                    energy, yield_fraction, radiation_type) for energy, yield_fraction, radiation_type in nuc_emissions]
            db.add_nuclide(nuc)

    if add_single_nuclides:
        half_life_data = half_life_from_icrp107()
        nuc_names_in_db = set(db.nuclides.keys())
        nuc_names_in_edata = set(emission_data.keys())

        non_overlapping = nuc_names_in_edata - nuc_names_in_db

        for nuc_name in non_overlapping:
            name = nuc_name
            halflife = half_life_data[nuc_name]
            daughters = []
            emissions = emission_data[nuc_name]
            nuclide = Nuclide(name, halflife, daughters)
            nuclide.emissions = [RadiationEmission(
                energy, yield_fraction, radiation_type) for energy, yield_fraction, radiation_type in emissions
            ]
            db.add_nuclide(nuclide)




    return db

def nuc_energy_lookup(nuc_name, beta_energy_data, alpha_energy_data):

    # See if the nuclide is in the beta data

    beta_nuc = beta_energy_data[beta_energy_data["Nuclide"] == nuc_name]
    alpha_nuc = alpha_energy_data[alpha_energy_data["Nuclide"] == nuc_name]

    decay_energy_radiated = {"Name": nuc_name,
                             "Alpha Energy (MeV)": 0,
                             "Beta Energy (MeV)": 0,
                             "Total Energy (MeV)": 0}

    total_energy = 0  # In MeV

    if beta_nuc.empty and alpha_nuc.empty:
        print("Not in list of beta or alpha decay data", nuc_name)
        return decay_energy_radiated

    if not beta_nuc.empty:

        beta_yield = float(beta_nuc["Yield"].values[0])

        beta_energy = float(beta_nuc["Mean E (MeV)"].values[0])
        decay_energy_radiated["Beta Energy (MeV)"] = beta_energy*beta_yield
        total_energy += beta_energy

    if not alpha_nuc.empty:

        alpha_energy = float(alpha_nuc["Energy (MeV/nt)"].values[0])
        decay_energy_radiated["Alpha Energy (MeV)"] = alpha_energy
        total_energy += alpha_energy

    decay_energy_radiated["Total Energy (MeV)"] = total_energy

    return decay_energy_radiated

def energy_in_decay_chain(db, nuc_name, parent_fraction=1.0):

    total_alpha_energy = 0
    total_beta_energy = 0
    total_non_penetrative_energy = 0

    if db.get_decay_info(nuc_name) is None:
        print("Nuclide not found in database", nuc_name)
        return {"Total Alpha Energy": 0,
                "Total Beta Energy": 0,
                "Total Non-Penetrative Energy": 0}

    total_alpha_energy += db.get_decay_info(
        nuc_name).calculate_alpha_energy() * parent_fraction
    total_beta_energy += db.get_decay_info(
        nuc_name).calculate_beta_energy() * parent_fraction
    total_non_penetrative_energy += db.get_decay_info(
        nuc_name).calculate_total_non_penetrative_energy() * parent_fraction

    decays = (get_daughters(db, nuc_name, parent_fraction=parent_fraction))

    for daughter_name, fraction in decays.items():
        if db.get_decay_info(daughter_name) is None:
            continue

        total_alpha_energy += db.get_decay_info(
            daughter_name).calculate_alpha_energy() * fraction
        total_beta_energy += db.get_decay_info(
            daughter_name).calculate_beta_energy() * fraction
        total_non_penetrative_energy += db.get_decay_info(
            daughter_name).calculate_total_non_penetrative_energy() * fraction

    energy_dict = {"Total Alpha Energy": total_alpha_energy,
                   "Total Beta Energy": total_beta_energy,
                   "Total Non-Penetrative Energy": total_non_penetrative_energy}

    return energy_dict

def load_icrp_107(path="resources/FULL_RAD_LIST.RAD"):

    energy_data = defaultdict(list)
    current_nuclide = None

    # List of possible radiation types to match
    valid_radiation_types = ['X', 'A', 'G', 'IE', 'AE', 'AR']

    # Create a regex pattern that matches the valid radiation types
    # radiation_type_pattern = '|'.join(valid_radiation_types)  # e.g., "X|A|G|IE|AE|AR"
    radiation_type_pattern = r'\b(?:' + \
        '|'.join(valid_radiation_types) + r')\b'

    with open(path, 'r') as file:
        for line in file:
            # Match a nuclide entry (e.g., "Ac-223        2.10m")
            # nuclide_match = re.match(r"^([A-Za-z\-0-9]+)\s+([\d\.]+[mhd])", line)
            nuclide_match = None
            line_no_repeating_white_space = ' '.join(line.split())
            split_line = line_no_repeating_white_space.split(" ")
            if len(split_line) == 3:
                nuclide_match = True
                current_nuclide = split_line[0]  # Extract nuclide name
                line_no_repeating_white_space = ' '.join(line.split())
                split_line_nuc = line_no_repeating_white_space.split(" ")
                continue  # Move to next line

            # Match radiation data (e.g., "2 8.38497E+00 6.06280E-05  X")
            # radiation_match = re.match(r"^\s*2\s+([\d\.E+-]+)\s+([\d\.E+-]+)\s+([A-Za-z]+)", line)
            radiation_match = not nuclide_match
            if radiation_match and current_nuclide:

                if nuclide_match == None:

                    line = ' '.join(line.split())
                    split_line = line.split(" ")
                    if len(split_line) != 3:
                        yield_fraction = float(split_line[1])
                        energy = float(split_line[2])
                        radiation_type = split_line[3]

                        # Store energy information in the temporary dictionary
                        energy_data[current_nuclide].append(
                            (energy, yield_fraction, radiation_type))

    return energy_data

if __name__ == "__main__":

    # Load energy record data testing

    emission_energy = load_icrp_107()
    #half_life_data = half_life_from_icrp107()
    
    db = make_decay_chain_db(emission_data=emission_energy)

    nucs_in_db = db.nuclides.keys()

    # Hours
    half_life_lower = 2.0
    half_life_upper = 12*24

    for nuc_name in nucs_in_db:
        half_life = db.get_decay_info(nuc_name).halflife

        match = re.match(r"([0-9.]+(?:[eE][+-]?[0-9]+)?)([a-zA-Z]+)", half_life)

        half_life_hours = convert_to_hours(float(match.group(1)), match.group(2))
        if half_life_hours > half_life_lower and half_life_hours < half_life_upper:
            print(nuc_name, half_life_hours)

    print("Have a nice day")
