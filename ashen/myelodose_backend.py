"""

This is the back-end code for myelodose, that is going to spring out 
from ashen in order to handle red marrow dosimetry calculations.

This will inherit some functionality from the red marrow module,
but quite a lot will need to be rewritten to fit the new structure.

"""

import yaml
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sys
import seaborn as sns
import json
import pickle as pkl
from dataclasses import dataclass, field
import math
from datetime import datetime as dt

from jinja2 import Environment, PackageLoader, select_autoescape

from ashen.ashen_utils import (
    get_daughters,
    make_decay_chain_db,
    load_icrp_107,
    convert_to_hours,
    )

from ashen.beta_spectrum_analysis import (
    replace_simple_beta_with_full_spectrum
)

# --- Load predefined data -------------------------------------------------------

# Read json file for checkbox data

with open("predefined_data.json", "r") as f:
    json_data = json.load(f)

# --- Simulated data -------------------------------------------------------------

#SEVEN_FIELD_VALUES = ["Ac-225", "Pb-212", "At-211"]
#PREDEFINED_VALUES = ["Lu-177", "Y-90", "I-131"]
#PREDEFINED_VALUES += list(SEVEN_FIELD_VALUES)

# --- Real data ---------------------------------------------------------------

PREDEFINED_VALUES = list(json_data.get("CHECKBOX_TITLES", {}).keys())
SEVEN_FIELD_VALUES = json_data.get("SEVEN_FIELD_VALUES", []) # These are the alpha-emitters with 7 fields

# Field name lists for forms
#ELECTRON_SITES = [f"Site_electrons_{i}" for i in range(1, 14)]

ELECTRON_SITES = [
        "craniofacial_bones",
        "mandible",
        "scapulae",
        "clavicles",
        "sternum",
        "ribs",
        "cervical_vertebrae",
        "thoracic_vertebrae",
        "lumbar_vertebrae",
        "sacrum",
        "os_coxae",
        "proximal_humeri",
        "proximal_femora",
]

ALPHA_SITES = [
        "cervical_vertebrae",
        "femur_head",
        "femur_neck",
        "iliac_crest",
        "lumbar_vertebrae",
        "ribs",
        "parietal_bone"
]

SITES_BOTH_ALPHA_ELECTRON = list(set(ELECTRON_SITES) & set(ALPHA_SITES))

ELECTRON_SURROGATES= {
    "femur_head": "proximal_femora",
    "femur_neck": "proximal_femora",
    "iliac_crest": "os_coxae",
    "parietal_bone": "craniofacial_bones"
}



# -- Flags ---------------------------------------------------------------

REMAKE_DB = False
SILENCE_WARNING = True
USE_UNITY_PHI_FOR_ELECTRON_SURROGATE = True

with open("combined_saf.yaml", "r") as f:
    SKELETAL_SITE_DATA_ELECTRONS = yaml.safe_load(f)

# --- Skeletal site data - alphas ------------

with open("combined_alpha_saf.yaml", "r") as f:
    SKELETAL_SITE_DATA_ALPHAS = yaml.safe_load(f)

if REMAKE_DB:
    emission_energy = load_icrp_107()
    decay_chain_db = make_decay_chain_db(emission_data=emission_energy)

    # Pickle the db for later use
    with open("decay_chain_db.pkl", "wb") as f:
        pkl.dump(decay_chain_db, f)

else:
    # Load the pickled db
    with open("decay_chain_db.pkl", "rb") as f:
        decay_chain_db = pkl.load(f)

CUT_OFF_DAUGHTER_HOURS = 2

# Global dict with ICRP-values for skeletal sites

ICRP_SKELETAL_SITE_VALUES_ELECTRONS = {
    "craniofacial_bones": 38,
    "mandible": 38,
    "scapulae": 38,
    "clavicles": 33,
    "sternum": 70,
    "ribs": 70,
    "cervical_vertebrae": 70,
    "thoracic_vertebrae": 70,
    "lumbar_vertebrae": 70,
    "sacrum": 70,
    "os_coxae": 48,
    "proximal_humeri": 25,
    "proximal_femora": 25,
}

source_tissue_dict = {
    "RM": "red_marrow",
    "TBS": "tbs"
} # TODO: Should be rewritten to match other parts of the code

# --- Make a dataclass to hold result of calculation
#     data ------------

@dataclass
class CalculationResult:
    """
    Dataclass to hold the results of a single calculation
    where "single" mean for a specific nuclide in a specific site.
    """
    alpha_RBE_value: float = 0.0
    parent_nuclide: str = ""
    radionuclide: str = ""
    source_tissue: str = ""
    name: str = ""
    CF: str = ""
    total_corrected_cumulative_activity_MBqhrs: float = 0.0
    only_electron_calculation: bool = False
    sites_compatible: bool = True
    include_in_final_results: bool = True
    absorbed_dose_Gy_electrons: float = 0.0
    absorbed_dose_Gy_alpha: float = 0.0
    absorbed_dose_Gy_alpha_rbe_adjusted: float = 0.0
    fraction_energy_absorbed_electrons: float = 0.0
    fraction_energy_absorbed_alpha: float = 0.0
    surrogate_electron_site: str = ""
    surrogate_electron_site_used: bool = False
    electron_unity_saf_used: bool = False

@dataclass
class CombinedCalculationResults:
    results: list = field(default_factory=list)
    ordered_daughters: list = field(default_factory=list)

    def save_to_json(self, filename: str):
        with open(filename, 'w') as f:
            json.dump([result.__dict__ for result in self.results], f, indent=4)

    def prepare_rows_for_plotting(self) -> list:
        rows = []

        for result in self.results:
            if not result.include_in_final_results:
                continue
            rows.append({
                "radionuclide": result.radionuclide,
                "dose_alpha": result.absorbed_dose_Gy_alpha,
                "dose_alpha_rbe_adjusted": result.absorbed_dose_Gy_alpha_rbe_adjusted,
                "dose_electron": result.absorbed_dose_Gy_electrons,
                "source_tissue": result.source_tissue,
                "CF": result.CF,
                "alpha_RBE_value": result.alpha_RBE_value,
                "site": result.name,
                "surrogate_electron_site_used": result.surrogate_electron_site_used,
                "electron_unity_saf_used": result.electron_unity_saf_used,
                "parent_nuclide": result.parent_nuclide
            })

        return rows


def retrieve_reference_mass_target(path_to_data = "resources\site_volumes.xlsx"):

    mass_df = pd.read_excel(path_to_data)

    reference_masses = {}

    for index, row in mass_df.iterrows():
        site = row['Site']
        mass_g = row['Reference_RM']
        reference_masses[site] = mass_g

    return reference_masses

def mass_data_sites(path_to_data = "resources\site_volumes.xlsx"):

    mass_df = pd.read_excel(path_to_data)

    # Make a dict of dicts where site is key and sub-keys 
    # are volume and mass

    mass_data = {}

    for index, row in mass_df.iterrows():
        site = row['Site']
        spongiosa_volume_ml = row['Spongiosa volume']
        tbv_fraction = row['Trab bone fraction']
        marrow_mass_g = row['Marrow Mass']
        icrp_cf = row['ICRP_CF']

        mass_data[site] = {
            "spongiosa_volume_ml": spongiosa_volume_ml,
            "tbv_fraction": tbv_fraction,
            "total_marrow_mass_g": marrow_mass_g,
            "ICRP_CF": icrp_cf

        }

    return mass_data

def get_saf_data_electrons(site: str, 
                           source_tissue: str,
                           CF: int,
                           skeletal_data = SKELETAL_SITE_DATA_ELECTRONS):

    """

    Return the relevant SAF data for a given:

    - site
    - CF
    - Source tissue
    
    :param site: Name of the skeletal site
    :param source_tissue: Name of the source tissue
    :param CF: Cellularity factor
    :param skeletal_data: The full skeletal data dictionary
    
    :return: SAF data in terms of specific absorbed fractions
    """

    # Check the CF is either ircrp or a valid number
    # Check also that CF is between 10 and 100 in 10 increments


    if CF != "icrp":
        try:
            cf_value = int(CF)
            if cf_value < 10 or cf_value > 100 or cf_value % 10 != 0:
                raise ValueError("CF must be between 10 and 100 in increments of 10, or 'icrp'")
        except ValueError:
            raise ValueError("CF must be an integer between 10 and 100 in increments of 10, or 'icrp'")
    
    if source_tissue == "tbs": # Ugly hack to handle TBS case
        saf_data = skeletal_data[site]["electrons"][source_tissue][0]
    else:
        saf_data = skeletal_data[site]["electrons"][source_tissue][CF]

    return saf_data

def get_af_data_alphas(site: str,
                          source_tissue: str,
                          CF: str,
                          skeletal_data = SKELETAL_SITE_DATA_ALPHAS):
    
    """
    
    Return the relevant SAF data for a given:
    
    - site
    - CF
    - Source tissue
    
    :param site: Name of the skeletal site
    :param source_tissue: Name of the source tissue
    :param CF: Cellularity factor
    :param skeletal_data: The full skeletal data dictionary
    
    :return: SAF data in terms of specific absorbed fractions
    """

    # Check if CF is in the correct format
    # 

    if CF == "icrp":
        CF = str(ICRP_SKELETAL_SITE_VALUES_ELECTRONS.get(site, None)) # Using the electron CFs


    if CF not in ["10", "20", "30", "40", "50", "60", "70", "80", "90", "100"]:
    
        raise ValueError("CF must be between 10 and 100 in increments of 10")
    
    af_data = skeletal_data[site]["alphas"][source_tissue][CF]
    
    return af_data

def interpolate_phi(energy: float, 
                    saf_data: dict,
                    interpolation_technique: str) -> float:

    if interpolation_technique not in ["closest", "linear", "spline", "loglog"]:
        raise ValueError(f"Interpolation technique {interpolation_technique} not recognized.")

    # Make numpy arrays from the saf_data dict

    E = saf_data.keys()
    Phi = saf_data.values()

    # Turn into numpy arrays

    E = np.array(list(E))
    Phi = np.array(list(Phi))

    if (energy < E.min() or energy > E.max()) and not SILENCE_WARNING:
        print("Warning: Energy is outside the range of the SAF data provided.")
    
    if interpolation_technique == "closest":
        closest_energy = E[np.abs(E - energy).argmin()]
        phi_value = saf_data[closest_energy]
        return phi_value

    elif interpolation_technique == "spline":
        from scipy.interpolate import interp1d
        f = interp1d(E, Phi, kind='cubic', bounds_error=False, fill_value="extrapolate")
        return float(f(energy))
    
    elif interpolation_technique == "loglog":
        return float(np.exp(np.interp(np.log(energy), np.log(E), np.log(Phi))))
    
    elif interpolation_technique == "linear":
        return float(np.interp(energy, E, Phi))
    
    else:
        raise ValueError("Interpolation technique not recognized.")

    return None

def _saf_from_emission_data(emission_data, 
                           site: str,
                           source_tissue: str,
                           CF: str,
                           interpolation_technique: str):
    
    """

    This is potentially a legacy function
    that will be removed later.
    
    :param emission_data: Description
    :param site: Description
    :type site: str
    :param source_tissue: Description
    :type source_tissue: str
    :param CF: Description
    :type CF: str
    :param interpolation_technique: Description
    :type interpolation_technique: str
    """
    
    electron_emission_types = ["B-", "IC", "AE"]

    saf_values = []

    saf_data = get_saf_data_electrons(site=site,
                                     source_tissue=source_tissue,
                                     CF=CF)

    for emission in emission_data:
        emission_type = emission.radiation_type
        if emission_type not in electron_emission_types:
            continue
        energy = emission.energy
        phi = interpolate_phi(energy=energy,
                              saf_data=saf_data,
                              interpolation_technique=interpolation_technique)
        saf_values.append((energy, phi))

    return saf_values

def correct_cumulative_activity(sites, source_tissue, input_units):

    errors = []

    mass_data = mass_data_sites()

    # TODO - Handle ICRP CF values here
    # Handle the units of the cumulative activity concentration here as well

    if input_units == "MBqhrs_per_ml":

        for site in sites:
            print(f"Processing site: {site}")
            cumulative_activity_conc = float(site.get("value", 0))

            cf_value = site.get("CF", None)

            if cf_value == "icrp":
                cf_value = ICRP_SKELETAL_SITE_VALUES_ELECTRONS.get(site.get("name", ""), None)/100.0
            else:
                cf_value = float(site.get("CF", None))/100.0 if site.get("CF", None) is not None else None

            print(f"Cumulative activity concentration: {cumulative_activity_conc} MBq·hrs/ml")
            print(f"Cellularity factor: {site.get('CF', 'N/A')}")

            spongiosa_volume = mass_data.get(site.get("name", ""), {}).get("spongiosa_volume_ml", None)

            if spongiosa_volume is None:
                print(f"Error: Spongiosa volume data not found for site {site.get('name', '')}. Cannot correct cumulative activity.")
                errors.append(f"Spongiosa volume data not found for site {site.get('name', '')}.")
                continue

            # Make a correction if source tissue is RM

            if source_tissue == "RM":
                print("Making correction for red marrow source tissue.")
                tbv_fraction = mass_data.get(site.get("name", ""), {}).get("tbv_fraction", None)
                marrow_fraction = (1 - tbv_fraction)*cf_value
                marrow_mass = marrow_fraction*1.03 # Hacky - must place this value someplace else
                corrected_cumulative_activity = cumulative_activity_conc/marrow_mass
                print(f"Corrected cumulative activity concentration: {corrected_cumulative_activity} MBq·hrs/ml")
                total_marrow_mass = mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", None)
                total_cumulative_activity = corrected_cumulative_activity * total_marrow_mass

                print(f"Total cumulative activity in site: {total_cumulative_activity} MBq·hrs")

            elif source_tissue == "TBS":
                    print("Using total bone source tissue - no correction applied.")

                    # Total trabecular bone source tissue
                    total_cumulative_activity = cumulative_activity_conc*mass_data.get(site.get("name", ""), {}).get("spongiosa_volume_ml", None) 

                    #total_marrow_mass = mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0)

            print(f"Total cumulative activity in site: {total_cumulative_activity} MBq·hrs")

            site['total_corrected_cumulative_activity_MBqhrs'] = total_cumulative_activity

        return sites
    
    if input_units == "MBqhrs":

        for site in sites:
            cumulative_activity = float(site.get("value", 0))
            site['total_corrected_cumulative_activity_MBqhrs'] = cumulative_activity

        return sites


def _calculate_absorbed_dose_electron(corr_sites, 
                           source_tissue: str,
                           nuclide,
                           calculation_input,
                           mass_data):
    
    """
    This is a potential legacy function that may be removed later.
    """

    for site in corr_sites:

        energy_emitted_in_site = 0
        energy_absorbed_in_site = 0

        if source_tissue == "RM":

            total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))
            total_red_marrow_mass = total_marrow_mass * (float(site.get("CF", 0))/100.0)

        elif source_tissue == "TBS":

            total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))
            icrp_cf = mass_data.get(site.get("name", ""), {}).get("ICRP_CF", None)
            total_red_marrow_mass = total_marrow_mass * icrp_cf

        electron_saf_data = get_saf_data_electrons(site=site.get("name", ""),
                                                    source_tissue=source_tissue_dict[source_tissue],
                                                    CF=site.get("CF", None))

        for em in nuclide.emissions:
            if em.radiation_type not in ["B-", "IE", "AE"]:
                continue
            energy = em.energy
            Phi = interpolate_phi(energy=energy,
                                  saf_data=electron_saf_data,
                                  interpolation_technique=calculation_input.get("interpolation_technique", "linear"))
            #print(f"Energy: {energy} MeV - Phi: {Phi} for site {site.get('name', '')}")

            energy_absorbed_in_site += energy * Phi * em.yield_fraction*total_red_marrow_mass
            energy_emitted_in_site += energy * em.yield_fraction

        total_energy_absorbed_in_MeV = energy_absorbed_in_site * site.get('total_corrected_cumulative_activity_MBqhrs', 0) * 3600*1e6 
        total_energy_absorbed_in_J = total_energy_absorbed_in_MeV * 1.60218e-13 # TODO: Magic number - place elsewhere

        absorbed_dose_Gy = total_energy_absorbed_in_J / (total_red_marrow_mass * 1e-3)  # mass in kg

        site["absorbed_dose_Gy_electrons"] = absorbed_dose_Gy
        site["fraction_energy_absorbed_electrons"] = energy_absorbed_in_site/energy_emitted_in_site if energy_emitted_in_site > 0 else 0

    return corr_sites

def calculate_absorbed_dose_alpha(corr_sites,
                                  source_tissue: str,
                                  nuclide,
                                  calculation_input,
                                  mass_data):
    
    for site in corr_sites:

        energy_emitted_in_site = 0
        energy_absorbed_in_site = 0

        alpha_af_data = get_af_data_alphas(site=site.get("name", ""),
                                             source_tissue=source_tissue_dict[source_tissue],
                                             CF=site.get("CF", None))

        print(f"Calculating alpha dose for site: {site.get('name', '')}")

        total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))

        if total_marrow_mass == 0:
            print(f"Warning: Total marrow mass for site {site.get('name', '')} is zero. Skipping dose calculation.")
            continue

        total_red_marrow_mass = total_marrow_mass * (float(site.get("CF", 0))/100.0)

        for em in nuclide.emissions:
            if em.radiation_type not in ["A"]:
                continue

            energy = em.energy

            AF = interpolate_phi(energy=energy,
                                  saf_data=alpha_af_data,
                                  interpolation_technique=calculation_input.get("interpolation_technique", "linear"))

            print(f"Energy: {energy} MeV - AF: {AF} for site {site.get('name', '')}")

            energy_absorbed_in_site += energy * AF * em.yield_fraction
            energy_emitted_in_site += energy * em.yield_fraction

        total_energy_absorbed_in_MeV = energy_absorbed_in_site * site.get('total_corrected_cumulative_activity_MBqhrs', 0) * 3600*1e6 
        total_energy_absorbed_in_J = total_energy_absorbed_in_MeV * 1.60218e-13 # TODO: Magic number - place elsewhere

        absorbed_dose_Gy = total_energy_absorbed_in_J / (total_red_marrow_mass * 1e-3)  # mass in kg
        site["absorbed_dose_Gy_alpha"] = absorbed_dose_Gy
        site["fraction_energy_absorbed_alpha"] = energy_absorbed_in_site/energy_emitted_in_site if energy_emitted_in_site > 0 else 0

    return corr_sites

def electron_dose_to_site(
    site: dict,
    nuclide,
    source_tissue: str,
    mass_data,
    branching_ratio: float = 1.0,
    use_unity_saf: bool = False
):
    
    # This should return absorbed dose from electrons 
    # and the fraction of energy absorbed

    # Placeholder for future implementation

    energy_emitted_in_site = 0
    energy_absorbed_in_site = 0

    energy_emitted_by_type = {}
    energy_absorbed_by_type = {}

    cf = site.get("CF", None)

    if cf == "icrp":
        cf_value = ICRP_SKELETAL_SITE_VALUES_ELECTRONS.get(site.get("name", ""), None)/100.0
    else:
        cf_value = float(site.get("CF", None))/100.0 if site.get("CF", None) is not None else None

    if source_tissue == "RM":

        total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))
        total_red_marrow_mass = total_marrow_mass * cf_value

    elif source_tissue == "TBS":

        total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))
        icrp_cf = mass_data.get(site.get("name", ""), {}).get("ICRP_CF", None)
        total_red_marrow_mass = total_marrow_mass * icrp_cf

    electron_saf_data = get_saf_data_electrons(site=site.get("name", ""),
                                                source_tissue=source_tissue_dict[source_tissue],
                                                CF=cf)
    
    for em in nuclide.emissions:
        if em.radiation_type not in ["B-", "IE", "AE"]:
            continue
        energy = em.energy

        if use_unity_saf: # Override for surrogate sites
            Phi = 1.0/total_red_marrow_mass
        else:
            Phi = interpolate_phi(energy=energy, # TODO: Should get this from calculation_input
                                  saf_data=electron_saf_data,
                                  interpolation_technique="linear")
        

        energy_absorbed_in_site += energy * Phi * em.yield_fraction*total_red_marrow_mass
        energy_emitted_in_site += energy * em.yield_fraction

        energy_emitted_by_type[em.radiation_type] = energy_emitted_by_type.get(em.radiation_type, 0) + energy * em.yield_fraction
        energy_absorbed_by_type[em.radiation_type] = energy_absorbed_by_type.get(em.radiation_type, 0) + energy * Phi * em.yield_fraction*total_red_marrow_mass

    total_energy_absorbed_in_MeV = energy_absorbed_in_site * site.get('total_corrected_cumulative_activity_MBqhrs', 0) * 3600*1e6*branching_ratio
    total_energy_absorbed_in_J = total_energy_absorbed_in_MeV * 1.60218e-13 # TODO: Magic number - place elsewhere
    absorbed_dose_Gy = total_energy_absorbed_in_J / (total_red_marrow_mass * 1e-3)  # mass in kg

    result_dict = {}

    result_dict['absorbed_dose_Gy'] = absorbed_dose_Gy
    result_dict['fraction_energy_absorbed'] = energy_absorbed_in_site/energy_emitted_in_site if energy_emitted_in_site > 0 else 0
    result_dict['energy_emitted_by_type'] = energy_emitted_by_type
    result_dict['energy_absorbed_by_type'] = energy_absorbed_by_type

    return result_dict

def alpha_dose_to_site(
    site: dict,
    nuclide,
    source_tissue: str,
    mass_data,
    branching_ratio: float = 1.0,
    rbe_alpha_value: float = 1.0
):
    
    energy_emitted_in_site = 0
    energy_absorbed_in_site = 0

    cf = site.get("CF", None)

    if cf == "icrp":
        cf_value = ICRP_SKELETAL_SITE_VALUES_ELECTRONS.get(site.get("name", ""), None)/100.0
    else:
        cf_value = float(site.get("CF", None))/100.0 if site.get("CF", None) is not None else None

    alpha_af_data = get_af_data_alphas(site=site.get("name", ""),
                                         source_tissue=source_tissue_dict[source_tissue],
                                            CF=cf)
    
    total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))

    if total_marrow_mass == 0:
        print(f"Warning: Total marrow mass for site {site.get('name', '')} is zero. Skipping dose calculation.")
        return None
    
    total_red_marrow_mass = total_marrow_mass * cf_value

    for em in nuclide.emissions:
        if em.radiation_type not in ["A"]:
            continue

        energy = em.energy

        AF = interpolate_phi(energy=energy,
                              saf_data=alpha_af_data,
                              interpolation_technique="linear")

        print(f"Energy: {energy} MeV - AF: {AF} for site {site.get('name', '')}")

        energy_absorbed_in_site += energy * AF * em.yield_fraction
        energy_emitted_in_site += energy * em.yield_fraction

    total_energy_absorbed_in_MeV = energy_absorbed_in_site * site.get('total_corrected_cumulative_activity_MBqhrs', 0) * 3600*1e6*branching_ratio
    total_energy_absorbed_in_J = total_energy_absorbed_in_MeV * 1.60218e-13 # TODO: Magic number - place elsewhere

    absorbed_dose_Gy = total_energy_absorbed_in_J / (total_red_marrow_mass * 1e-3)  # mass in kg
    absorbed_dose_Gy_rbe_adjusted = absorbed_dose_Gy * rbe_alpha_value

    result_dict = {}

    result_dict['absorbed_dose_Gy'] = absorbed_dose_Gy
    result_dict['fraction_energy_absorbed'] = energy_absorbed_in_site/energy_emitted_in_site if energy_emitted_in_site > 0 else 0
    result_dict["energy_emitted_in_site"] = energy_emitted_in_site
    result_dict["energy_absorbed_in_site"] = energy_absorbed_in_site
    result_dict["absorbed_dose_Gy_rbe_adjusted"] = absorbed_dose_Gy_rbe_adjusted
    
    return result_dict
    
def calculate_absorbed_dose_to_site(
    site: dict,
    nuclide,
    calculation_input,
    mass_data,
    branching_ratio: float = 1.0,
    rbe_alpha_value: float = 1.0
):
    
    # Should return a calculation result dataclass instance
    
    calculation_results = CalculationResult()

    if calculation_input.get("only_electron_calculation", False):

        result_tmp_electrons = electron_dose_to_site(
            site=site,
            nuclide=nuclide,
            source_tissue=calculation_input.get("source_tissue", ""),
            mass_data=mass_data,
            branching_ratio=branching_ratio
        )

    else:

        # TODO - here we need to handle when there are incompatible sites
        # double check that this logic works as intended

        # Switch statement where one out three cases are handled:
        # Case 1: Site is compatible with both alpha and electron dose calculation - calculate both as normal
        # Case 2: Site is not compatible, and a surrogate is used
        # Case 3: Site is not compatible, a surrogate is used, and the surrogate is calculated with unity SAF (for electrons)

        if site.get("name", "") in ELECTRON_SURROGATES.keys():
            surrogate_site_name = ELECTRON_SURROGATES[site.get("name", "")]
            print(f"Using surrogate electron site {surrogate_site_name} for alpha site {site.get('name', '')}")
            surrogate_site = site.copy()
            surrogate_site["name"] = surrogate_site_name
            calculation_results.surrogate_electron_site = surrogate_site_name
            calculation_results.surrogate_electron_site_used = True

            if USE_UNITY_PHI_FOR_ELECTRON_SURROGATE:
                print("Using unity AF for surrogate electron site calculation.")
                calculation_results.electron_unity_saf_used = True

            result_tmp_electrons = electron_dose_to_site(
                site=surrogate_site,
                nuclide=nuclide,
                source_tissue=calculation_input.get("source_tissue", ""),
                mass_data=mass_data,
                branching_ratio=branching_ratio,
                use_unity_saf=USE_UNITY_PHI_FOR_ELECTRON_SURROGATE
            )

        else:

            result_tmp_electrons = electron_dose_to_site(
                site=site,
                nuclide=nuclide,
                source_tissue=calculation_input.get("source_tissue", ""),
                mass_data=mass_data,
                branching_ratio=branching_ratio
            )

        result_tmp_alpha = alpha_dose_to_site(
            site=site,
            nuclide=nuclide,
            source_tissue=calculation_input.get("source_tissue", ""),
            mass_data=mass_data,
            branching_ratio=branching_ratio,
            rbe_alpha_value=rbe_alpha_value
        )

    # Populate the calculation results dataclass

    calculation_results.parent_nuclide = calculation_input.get("radionuclide", "")
    calculation_results.radionuclide = nuclide.name
    calculation_results.source_tissue = calculation_input.get("source_tissue", "")
    calculation_results.name = site.get("name", "")
    calculation_results.CF = site.get("CF", "")

    calculation_results.total_corrected_cumulative_activity_MBqhrs = site.get("total_corrected_cumulative_activity_MBqhrs", 0)

    calculation_results.only_electron_calculation = calculation_input.get("only_electron_calculation", False)
    calculation_results.sites_compatible = calculation_input.get("sites_compatible", True)

    if result_tmp_electrons is not None:
        calculation_results.absorbed_dose_Gy_electrons = result_tmp_electrons.get("absorbed_dose_Gy", 0)
        calculation_results.fraction_energy_absorbed_electrons = result_tmp_electrons.get("fraction_energy_absorbed", 0)

    if not calculation_input.get("only_electron_calculation", False):
        if result_tmp_alpha is not None:
            calculation_results.absorbed_dose_Gy_alpha = result_tmp_alpha.get("absorbed_dose_Gy", 0)
            calculation_results.absorbed_dose_Gy_alpha_rbe_adjusted = result_tmp_alpha.get("absorbed_dose_Gy_rbe_adjusted", 0)
            calculation_results.fraction_energy_absorbed_alpha = result_tmp_alpha.get("fraction_energy_absorbed", 0)

    return calculation_results

def calculate_absorbed_dose_to_chain(
    site: dict,
    calculation_input,
    mass_data,
    rbe_alpha_value: float = 1.0
):
    
    # Calculate the absorbed dose from the entire decay chain

    if len(calculation_input.get("daughters", [])) == 0:
        return calculate_absorbed_dose_to_site(
            site=site,
            nuclide=decay_chain_db.get_decay_info(calculation_input.get("radionuclide", "")),
            calculation_input=calculation_input,
            mass_data=mass_data,
            rbe_alpha_value=rbe_alpha_value
        )
    
    else:

        combined_result = CombinedCalculationResults()

        daughter_dict = calculation_input.get("daughters", {})
        parent_nuclide = decay_chain_db.get_decay_info(calculation_input.get("radionuclide", ""))

        branching_ratios = get_daughters(decay_chain_db, parent_nuclide.name)

        for daughter, calc_daughter in zip(daughter_dict.keys(), daughter_dict.values()):

            if calc_daughter == 1:

                daughter_nuclide = decay_chain_db.get_decay_info(daughter)

                daughter_result = calculate_absorbed_dose_to_site(
                    site=site,
                    nuclide=daughter_nuclide,
                    calculation_input=calculation_input,
                    mass_data=mass_data,
                    branching_ratio=branching_ratios.get(daughter, None),
                    rbe_alpha_value=rbe_alpha_value
                )

                combined_result.results.append(daughter_result)

            print(f"Doing daughter: {daughter} - Included: {calc_daughter}")

        # Add parent nuclide calculation as well

        parent_result = calculate_absorbed_dose_to_site(
            site=site,
            nuclide=parent_nuclide,
            calculation_input=calculation_input,
            mass_data=mass_data,
            branching_ratio=1.0,
            rbe_alpha_value=rbe_alpha_value
        )

        combined_result.results.append(parent_result)


        # Adjust with RBE-value for alpha dose
        # TODO - why is the comment above here?

        return combined_result

def calculate_absorbed_dose_from_input_data(
    calculation_input: dict
):

    sites = calculation_input.get("fields", [])
    source_tissue = calculation_input.get("source_tissue", None)

    corr_sites = correct_cumulative_activity(sites, source_tissue, calculation_input.get("input_unit", None))

    nuclide = decay_chain_db.get_decay_info(calculation_input.get("radionuclide", ""))

    if nuclide.is_beta_emitter: # TODO IMPORTANT! This should not be done here, but deeper down the chain

        nuclide = replace_simple_beta_with_full_spectrum(nuclide)

    mass_data = mass_data_sites()

    try :
        rbe_alpha_value = float(calculation_input.get("alpha_RBE_value", 1.0))
    except TypeError:
        rbe_alpha_value = 1.0

    combined_results = CombinedCalculationResults()

    for site in corr_sites:

        calc_result = calculate_absorbed_dose_to_chain(
            site=site,
            calculation_input=calculation_input,
            mass_data=mass_data,
            rbe_alpha_value=rbe_alpha_value
        )

        if isinstance(calc_result, CombinedCalculationResults):
            combined_results.results.extend(calc_result.results)
        else:
            combined_results.results.append(calc_result)

    combined_results.ordered_daughters = [calculation_input.get("radionuclide", "")] + list(calculation_input.get("daughters", {}).keys())

    return combined_results

def post_process_back_end_inputs(input):

    #pros_fields = [f for f in input["fields"] if f["MBqhrs_per_ml"]]
    pros_fields = [f for f in input["fields"] if f["value"]]

    input["fields"] = pros_fields
    
    for k, v in input.items():
        print(f"{k}: {v}")

    if len(pros_fields) == 0:
        raise ValueError("No valid fields with cumulative activity concentration provided.")
    
    if not input["radionuclide"]:
        raise ValueError("No radionuclide specified in input.")

    # Check that the CF values are valid

    for field in input["fields"]:
        cf_value = field.get("CF", None)

        if field.get("name", "") in ELECTRON_SURROGATES.keys():
            print(f"Warning: Site {field.get('name', '')} is an alpha-only site. Using surrogate electron site {ELECTRON_SURROGATES[field.get('name', '')]}.")
            field["CF"] = field.get("CF", None)  # Keep the same CF for surrogate site
            continue
        else:
            if cf_value == str(ICRP_SKELETAL_SITE_VALUES_ELECTRONS[field.get("name", "")]):
                field["CF"] = "icrp"
                continue

        if cf_value is None:
            raise ValueError(f"Cellularity factor (CF) not specified for site {field.get('name', '')}.")
        try:
            cf_int = int(cf_value)
            if (cf_int < 10 or cf_int > 100 or cf_int % 10 != 0):
                raise ValueError(f"Cellularity factor (CF) for site {field.get('name', '')} must be between 10 and 100 in increments of 10.")
        except ValueError:
            raise ValueError(f"Cellularity factor (CF) for site {field.get('name', '')} must be an integer between 10 and 100 in increments of 10.")

    nuclide = input["radionuclide"]
    daughters = input["daughters"]
    full_chain = [nuclide] + list(daughters.keys())

    if all(d not in SEVEN_FIELD_VALUES for d in full_chain):
        print("\nNote: Selected radionuclide is not an alpha emitter. No alpha-calculation needed.")
        input["only_electron_calculation"] = True
        return input
    
    else:
        input["only_electron_calculation"] = False
        # Now check if sites are compatible with alpha calculation

        if "iliac_crest" in [f['name'] for f in input['fields'] if f['value']]:
            print("Warning: 'iliac_crest' site is extra problematic")
            input['incompatible_sites'] = ['iliac_crest']
            input['sites_compatible'] = False
            return input

        sites = [f['name'] for f in input['fields'] if f['value']]

        incompatible_sites = [s for s in sites if s not in SITES_BOTH_ALPHA_ELECTRON]

        if len(incompatible_sites) > 0:
            print("\nWarning: The following selected sites are incompatible with alpha and electron calculations ")
            for s in incompatible_sites:
                print(f" - {s}")
            input['incompatible_sites'] = incompatible_sites
            input['sites_compatible'] = False
            return input

        else:
            input['sites_compatible'] = True
            return input

    return input

def make_plot_figure(calc_results,
                     plot_rbe_adjusted_alpha: bool = False):

    df = pd.DataFrame(calc_results.prepare_rows_for_plotting())

    ordered_daughters = calc_results.ordered_daughters
    df = df.copy()
    df['radionuclide'] = pd.Categorical(df['radionuclide'],
                                        categories=ordered_daughters,
                                        ordered=True)

    df = df.sort_values(['site', 'radionuclide'])

    # Check if there is only a single nuclide

    if df['radionuclide'].nunique() == 1 and df["dose_alpha"].sum() == 0:
        print("Only a single radionuclide present in results.")
        generate_single_plot = True
    else:
        generate_single_plot = False

    surrogate_axes = []

    if generate_single_plot:
        fig, ax = plt.subplots(figsize=(8, 6))

        x = np.arange(len(df))

        ax.bar(
            x,
            df['dose_electron'],
            label='Electron AD'
        )

        ax.set_xticks(x)
        ax.set_xticklabels(df['site'], rotation=45)
        ax.set_title("Absorbed Dose per Site " + "(" + df['radionuclide'].iloc[0] + ")")
        ax.set_ylabel('Absorbed Dose (Gy)')
        ax.legend(frameon=False)

        fig.tight_layout()

        return fig

    sites = df['site'].unique()

    num_sites = len(sites)

    sqrt_num_sites = math.sqrt(num_sites)

    if sqrt_num_sites.is_integer():
        ncols = int(sqrt_num_sites)
    else:
        ncols = int(sqrt_num_sites) + 1

    fig, axes = plt.subplots(
        nrows=math.ceil(num_sites / ncols),
        ncols=ncols,
        figsize=(5 * ncols, 4 * math.ceil(num_sites / ncols)),
        sharey=True
    )

    alpha_dose_to_plot = "dose_alpha_rbe_adjusted" if plot_rbe_adjusted_alpha else "dose_alpha"

    if len(sites) == 1:
        axes = [axes]

    for ax, site in zip(fig.get_axes(), sites):
        site_df = df[df['site'] == site]

        x = np.arange(len(site_df))

        ax.bar(
            x,
            site_df[alpha_dose_to_plot],
            label='Alpha AD'
        )

        if site_df['surrogate_electron_site_used'].any():

            if site_df['electron_unity_saf_used'].any():

                ax.bar(
                    x,
                    site_df['dose_electron'],
                    bottom=site_df[alpha_dose_to_plot],
                    label='Electron AD (unity AF)',
                    hatch='//',
                    color='tab:gray'
                )

            else:
                ax.bar(
                    x,
                    site_df['dose_electron'],
                    bottom=site_df[alpha_dose_to_plot],
                    label='Electron AD (surrogate)',
                    hatch='//',
                    edgecolor='black',
                    fill=False
                )

            surrogate_axes.append(ax)

        else:

            ax.bar(
                x,
                site_df['dose_electron'],
                bottom=site_df[alpha_dose_to_plot],
                label='Electron AD',
            )

        ax.set_xticks(x)
        ax.set_xticklabels(site_df['radionuclide'], rotation=45)
        ax.set_title(site)
        ax.set_ylabel('Absorbed Dose (Gy)')

    fig.get_axes()[0].legend(frameon=False)
    fig.tight_layout()

    if len(surrogate_axes) > 0:
        for ax in surrogate_axes:
            handles, labels = ax.get_legend_handles_labels()
            new_handles = []
            new_labels = []
            for handle, label in zip(handles, labels):
                if label not in new_labels:
                    new_handles.append(handle)
                    new_labels.append(label)
            ax.legend(new_handles, new_labels, frameon=False)

    return fig

def check_for_warnings(calc_result):

    print("\nChecking for warnings with back-end-function...")
    print(decay_chain_db.get_decay_info(calc_result.ordered_daughters[0]))
    # TODO: Add the daughters when back
    warnings = []

    rows = calc_result.prepare_rows_for_plotting()

    for row in rows:
        if row.get("surrogate_electron_site_used", False) and row.get("electron_unity_saf_used", False):
            total_ad = row["dose_alpha_rbe_adjusted"] + row["dose_electron"]
            electron_fraction = row["dose_electron"] / total_ad if total_ad > 0 else 0
            warnings.append(f"Surrogate electron site and unity SAF used for {row['site']} - check if this is appropriate. Electrons contribute {electron_fraction:.1%} of total absorbed dose.")

    parent_nuclide = rows[0]['parent_nuclide'] if len(rows) > 0 else "N/A"

    parent_nuclide_info = decay_chain_db.get_decay_info(parent_nuclide)


    included_daughters = set()
    included_daughters.add(parent_nuclide)

    half_lives_of_daughters = decay_chain_db.get_half_lives_of_daughters(parent_nuclide)
    half_lives_of_daughters[parent_nuclide] = parent_nuclide_info.halflife

    for i in range(0, len(rows)):
        daughter_nuclide = rows[i]["radionuclide"]
        included_daughters.add(daughter_nuclide)

    for daughter in included_daughters:
        print(f"Daughter nuclide: {daughter}")

    daughter_half_lives_hours = {}

    for daughter in included_daughters:
        half_life = half_lives_of_daughters.get(daughter, "Unknown")
        half_life_hours = convert_to_hours(float(half_life[:-1]), half_life[-1])
        daughter_half_lives_hours[daughter] = half_life_hours

    print("\nDaughter nuclides and their half-lives in hours:")

    parent_half_life_hours = daughter_half_lives_hours[parent_nuclide]

    if parent_half_life_hours is not None:
        for daughter in included_daughters:
            daughter_half_life_hours = daughter_half_lives_hours.get(daughter)
            if daughter_half_life_hours is not None and daughter_half_life_hours > parent_half_life_hours:
                print(f"Warning: Daughter nuclide {daughter} has a longer half-life ({daughter_half_life_hours:.20f} hours) than the parent nuclide ({parent_half_life_hours:.20f} hours). This may affect the accuracy of the dose calculations.")

    # Check if any daughters have half life above cutoff of some set value

    for daughter in included_daughters:
        daughter_half_life_hours = daughter_half_lives_hours.get(daughter)
        if daughter_half_life_hours is not None and daughter_half_life_hours > CUT_OFF_DAUGHTER_HOURS:
            if daughter != parent_nuclide:
                warnings.append(f"Warning: Daughter nuclide {daughter} has a half-life of {daughter_half_life_hours:.2f} hours, which exceeds the cutoff of {CUT_OFF_DAUGHTER_HOURS} hours. Carefully consider this daughter in the calculation.")

    #other_radiation_type_warnings = []

    for daughter in included_daughters:
        daughter_info = decay_chain_db.get_decay_info(daughter)
        if daughter_info is not None:
            #warnings.append(check_non_electron_or_alpha_energies(daughter_info))
            # Concatenate lists instead
            warnings = warnings + check_non_electron_or_alpha_energies(daughter_info)

    return warnings

def check_non_electron_or_alpha_energies(nuclide, warning_threshold = 0.001):

    non_electron_alpha_energy_warnings = []

    emission_energies = {}
    emission_energies_included_parties = {}

    for em in nuclide.emissions:
        if em.radiation_type not in ["A", "B-", "IE", "AR", "AE"] and em.energy > 0:
            if em.radiation_type not in emission_energies.keys():
                emission_energies[em.radiation_type] = 0
            energy_emitted = em.energy*em.yield_fraction
            emission_energies[em.radiation_type] += energy_emitted

        if em.radiation_type in ["A", "B-", "IE", "AE"] and em.energy > 0:
            if em.radiation_type not in emission_energies_included_parties.keys():
                emission_energies_included_parties[em.radiation_type] = 0
            energy_emitted = em.energy*em.yield_fraction
            emission_energies_included_parties[em.radiation_type] += energy_emitted

    print(f"Nuclide: {nuclide.name}, Emission energies: {emission_energies}")

    # Fraction of energy from non-electron and non-alpha emissions compared to total energy from all emissions

    total_energy = sum(emission_energies_included_parties.values()) + sum(emission_energies.values())

    for key in emission_energies.keys():
        energy_fraction = emission_energies[key] / total_energy
        if energy_fraction > warning_threshold:
            non_electron_alpha_energy_warnings.append(f"Warning: {key} emissions contribute {energy_fraction:.2%} of the total energy emitted for nuclide {nuclide.name}, which exceeds the threshold of {warning_threshold:.2%}. Consider this in the absorbed dose calculations.")

    return non_electron_alpha_energy_warnings

def build_calculation_report(
        calc_results):
    
    report_body = ""

    datetime = dt.now().strftime("%Y-%m-%d %H:%M:%S")
    
    return None
    



#calc_result = calculate_absorbed_dose_to_chain(
#    site=corr_sites[0],
#    calculation_input=calculation_input,
#    mass_data=mass_data
#)
#
#print("Calculation result for first site:")
#print(calc_result)
#calc_result.save_to_json("calculation_result_test.json")

# Have the elements to perform the electron calculation now

#corr_sites = correct_cumulative_activity(sites, source_tissue)
#
#nuclide = decay_chain_db.get_decay_info(calculation_input.get("radionuclide", ""))
#
#if nuclide.is_beta_emitter:
#
#    nuclide = replace_simple_beta_with_full_spectrum(nuclide)
#
#print("Corrected sites data:")
#print(corr_sites)
#
#mass_data = mass_data_sites()

#calculate_absorbed_dose_alpha(corr_sites,
#                              source_tissue,
#                              nuclide,
#                              calculation_input,
#                              mass_data)

# Dump the results to a JSON file

#with open("calculation_results.json", "w") as f:
#    json.dump(calculation_input, f, indent=4)


#final_sites = calculate_absorbed_dose_electron(corr_sites,
#                                      source_tissue,
#                                        nuclide,
#                                        calculation_input,
#                                        mass_data)
#
#print("Final absorbed dose results:")
#for site in final_sites:
#    print(f"Site: {site.get('name', '')} - Absorbed dose (Gy): {site.get('absorbed_dose_Gy', 0)} - Fraction energy absorbed: {site.get('fraction_energy_absorbed', 0)}")

#ircrp_masses = retrieve_reference_mass_target()

#--- Skeletal site data - electrons ------------




    


# Test the saf-function

#saf_values = get_saf_data_electrons(site="lumbar_vertebrae",
#                       source_tissue="red_marrow",
#                          CF="20")
#
#foo = interpolate_phi(energy=0.00005,
#                saf_data=saf_values,
#                interpolation_technique="spline")
#
#print(f"Interpolated phi value: {foo}")
#
#print(saf_values)

#print(test_energy_inpolation())
#plot_skeletal_site_data_electrons(sites = sites)
