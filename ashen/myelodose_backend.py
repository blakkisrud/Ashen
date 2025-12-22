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

from ashen.ashen_utils import (
    get_daughters,
    make_decay_chain_db,
    load_icrp_107)

from ashen.beta_spectrum_analysis import (
    replace_simple_beta_with_full_spectrum
)

REMAKE_DB = False
SILENCE_WARNING = True

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

source_tissue_dict = {
    "RM": "red_marrow",
    "TBS": "tbs"
} # TODO: Should be rewritten to match other parts of the code

# --- Make a dataclass to hold result of calculation
#     data ------------

@dataclass
class CalculationResult:
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
    fraction_energy_absorbed_electrons: float = 0.0
    fraction_energy_absorbed_alpha: float = 0.0

@dataclass
class CombinedCalculationResults:
    results: list = field(default_factory=list)

    def save_to_json(self, filename: str):
        with open(filename, 'w') as f:
            json.dump([result.__dict__ for result in self.results], f, indent=4)



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

        mass_data[site] = {
            "spongiosa_volume_ml": spongiosa_volume_ml,
            "tbv_fraction": tbv_fraction,
            "total_marrow_mass_g": marrow_mass_g
        }

    return mass_data

def get_saf_data_electrons(site: str, 
                           source_tissue: str,
                           CF: str,
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

def saf_from_emission_data(emission_data, 
                           site: str,
                           source_tissue: str,
                           CF: str,
                           interpolation_technique: str):
    
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


# Function to parse input from GUI and run calculations

#with open("calculation_input.json", "r") as f:
#    calculation_input = json.load(f)
#
#sites = calculation_input.get("fields", [])
#
#source_tissue = calculation_input.get("source_tissue", None)

def correct_cumulative_activity(sites, source_tissue):

    errors = []

    mass_data = mass_data_sites()

    for site in sites:
        print(f"Processing site: {site}")
        cumulative_activity_conc = float(site.get("MBqhrs_per_ml", 0))
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


def calculate_absorbed_dose_electron(corr_sites, 
                           source_tissue: str,
                           nuclide,
                           calculation_input,
                           mass_data):

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
    branching_ratio: float = 1.0
):
    
    # This should return absorbed dose from electrons 
    # and the fraction of energy absorbed

    # Placeholder for future implementation

    energy_emitted_in_site = 0
    energy_absorbed_in_site = 0

    energy_emitted_by_type = {}
    energy_absorbed_by_type = {}

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
    branching_ratio: float = 1.0
):
    
    energy_emitted_in_site = 0
    energy_absorbed_in_site = 0

    alpha_af_data = get_af_data_alphas(site=site.get("name", ""),
                                         source_tissue=source_tissue_dict[source_tissue],
                                            CF=site.get("CF", None))
    
    total_marrow_mass = float(mass_data.get(site.get("name", ""), {}).get("total_marrow_mass_g", 0))

    if total_marrow_mass == 0:
        print(f"Warning: Total marrow mass for site {site.get('name', '')} is zero. Skipping dose calculation.")
        return None
    
    total_red_marrow_mass = total_marrow_mass * (float(site.get("CF", 0))/100.0)

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

    result_dict = {}

    result_dict['absorbed_dose_Gy'] = absorbed_dose_Gy
    result_dict['fraction_energy_absorbed'] = energy_absorbed_in_site/energy_emitted_in_site if energy_emitted_in_site > 0 else 0
    result_dict["energy_emitted_in_site"] = energy_emitted_in_site
    result_dict["energy_absorbed_in_site"] = energy_absorbed_in_site
    
    return result_dict
    
def calculate_absorbed_dose_to_site(
    site: dict,
    nuclide,
    calculation_input,
    mass_data,
    branching_ratio: float = 1.0
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
            branching_ratio=branching_ratio
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
            calculation_results.fraction_energy_absorbed_alpha = result_tmp_alpha.get("fraction_energy_absorbed", 0)

    return calculation_results

def calculate_absorbed_dose_to_chain(
    site: dict,
    calculation_input,
    mass_data
):
    
    # Calculate the absorbed dose from the entire decay chain

    if len(calculation_input.get("daughters", [])) == 0:
        return calculate_absorbed_dose_to_site(
            site=site,
            nuclide=decay_chain_db.get_decay_info(calculation_input.get("radionuclide", "")),
            calculation_input=calculation_input,
            mass_data=mass_data
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
                    branching_ratio=branching_ratios.get(daughter, 1.0)
                )

                combined_result.results.append(daughter_result)

            print(f"Doing daughter: {daughter} - Included: {calc_daughter}")

        # Add parent nuclide calculation as well

        parent_result = calculate_absorbed_dose_to_site(
            site=site,
            nuclide=parent_nuclide,
            calculation_input=calculation_input,
            mass_data=mass_data,
            branching_ratio=1.0
        )

        combined_result.results.append(parent_result)

        return combined_result

def calculate_absorbed_dose_from_input_data(
    calculation_input: dict
):

    sites = calculation_input.get("fields", [])
    source_tissue = calculation_input.get("source_tissue", None)

    corr_sites = correct_cumulative_activity(sites, source_tissue)

    nuclide = decay_chain_db.get_decay_info(calculation_input.get("radionuclide", ""))

    if nuclide.is_beta_emitter:

        nuclide = replace_simple_beta_with_full_spectrum(nuclide)

    mass_data = mass_data_sites()

    combined_results = CombinedCalculationResults()

    for site in corr_sites:

        calc_result = calculate_absorbed_dose_to_chain(
            site=site,
            calculation_input=calculation_input,
            mass_data=mass_data
        )

        if isinstance(calc_result, CombinedCalculationResults):
            combined_results.results.extend(calc_result.results)
        else:
            combined_results.results.append(calc_result)

    return combined_results

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
