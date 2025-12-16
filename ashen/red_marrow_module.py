""""

Module to calculate red marrow absorbed doses using different skeletal models.

To implement:
Electron sources from Hough and O'Reilly
Alpha sources from Watchman et al.

Input   - Cumulative activity in a skeletal site - with or without cortical bone
        - Cellularity Factor (CF) - representing the fraction of red marrow in the skeletal site
        - Radionuclide decay data - from ICRP 107 or user defined
Output - Absorbed dose to red marrow in that skeletal site

For the physical macroscopic part of the calculation - we can use the UF-phantom, and 
put absorbed fractions for the alpha-particles in the same set of of sites

The overlap between the UF-phantom and the Watchman-publication is somewhat limited, but
these are safe to use:

Lumbar vertebrae - L1-L5
Cervical vertebrae - C1-C7
Pareietal bone - skull

Somewhat less safe is the illiac crest, this might be substituted with the "pelvis" in the UF-phantom



"""

import os
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import pandas as pd
import sys
from ashen.ashen_utils import load_icrp_107, make_decay_chain_db, get_daughters, Nuclide, RadiationEmission, EXTRA_ICRP_107_PATH, DecayDatabase, convert_to_hours
import seaborn as sns
import pickle
from ashen.beta_spectrum_analysis import compute_yields, extract_data, PATH_TO_BETA_FILES, replace_simple_beta_with_full_spectrum
from scipy.interpolate import interp1d
import uuid
from datetime import datetime
from tqdm import tqdm

DEBUG_MODE = True
USE_FULL_BETA_SPECTRA = True
RELOAD_DECAY_DB = False
USE_DEBUG_DB = False

benchmark_values_path = "ashen/resources/JH_benchmark.xlsx"
JH_dose_factors = pd.read_excel(benchmark_values_path, sheet_name='Ark1')

# Site list for electron-emitters

def return_list_of_sites_by_nuclide(nuclide: Nuclide) -> list[str]:
    """Return a list of skeletal sites for a given nuclide based on its emission types."""
    
    ELECTRON_EMITTER_LIST = [
        "cranifacial_bones",
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

    ALPHA_EMITTER_LIST = [
        "cervical_vertebrae",
        "femur_head",
        "femur_neck",
        "iliac_crest",
        "lumbar_vertebrae",
        "ribs",
        "parietal_bone"]

    beta_substitutes = {
        "ribs" : "ribs",
        "lumbar_vertebrae" : "lumbar_vertebrae",
        "cervical_vertebrae" : "cervical_vertebrae",
        "femur_head" : "proximal_femora",
        "femur_neck" : "proximal_femora",
        "iliac_crest" : "os_coxae",
    }

    sites = []
    if not nuclide.is_alpha:
        if nuclide.is_electron_emitter:
            sites.extend(ELECTRON_EMITTER_LIST)
            return sites
        else:
            return None
    if nuclide.is_alpha:
        sites.extend(ALPHA_EMITTER_LIST)
        if nuclide.is_electron_emitter:
            for site in ALPHA_EMITTER_LIST:
                if site in beta_substitutes:
                    sites.append(beta_substitutes[site])
    return list(set(sites))

TB_FRACTIONS = {
    "Lumbar": 0.102,
    "Ribs": 0.115,
}

if USE_DEBUG_DB:

    dummy_daughter = Nuclide(
        name="CC-124",
        halflife="1.0h",
        daughters=[],
        emissions=[
            RadiationEmission(energy=0.1, yield_fraction=1.0, radiation_type="B-"),
            RadiationEmission(energy=0.2, yield_fraction=1.0, radiation_type="B-"),
            RadiationEmission(energy=0.3, yield_fraction=1.0, radiation_type="B-"),
        ]
    )

    dummy_nuclide = Nuclide(
        name="CC-123",
        halflife="1.0h",
        # daughters should be a list of (daughter_name, branching_fraction)
        # store the name (string) here, not the Nuclide object itself
        daughters=[],
        #daughters=[(dummy_daughter.name, 0.5)],
        emissions=[
            RadiationEmission(energy=0.1, yield_fraction=1.0, radiation_type="B-"),
            RadiationEmission(energy=0.2, yield_fraction=1.0, radiation_type="B-"),
            RadiationEmission(energy=0.3, yield_fraction=1.0, radiation_type="B-"),
        ]
    )

    dummy_db = DecayDatabase()
    dummy_db.add_nuclide(dummy_nuclide)
    dummy_db.add_nuclide(dummy_daughter)

    print("Dummy database created with nuclides:", list(dummy_db.nuclides.keys()))
    dummy_db.get_all_nuclide_names()

def extra_data_icrp107_lookup(nuclide_name: str) -> str:
    """Look up the decay type for a given nuclide from the extra ICRP 107 data."""
    with open(EXTRA_ICRP_107_PATH, 'r') as f:
        lines = f.readlines()

    for line in lines:
        parts = line.split()
        if parts[0] == nuclide_name:
            decay_type = parts[2]
            return decay_type

    raise ValueError(f"Nuclide {nuclide_name} not found in extra ICRP 107 data")

@dataclass
class AbsorbedDoseResultSingleNuclide:
    nuclide_name: str
    source_tissue: str
    target_tissue: str
    site: str
    cellularity_factor: float
    alpha_ad_mGy_per_MBq_s: float
    electron_ad_mGy_per_MBq_s: float
    alpha_ad_Gy_total: float
    electron_ad_Gy_total: float
    full_beta_spectrum_used: bool
    interpolation_technique: str
    unique_id: str = ""

    def convert_to_dataframe(self) -> pd.DataFrame:

        data = {
            'Nuclide': self.nuclide_name,
            'Source Tissue': self.source_tissue,
            'Target Tissue': self.target_tissue,
            'Site': self.site,
            'Cellularity Factor': self.cellularity_factor,
            'Alpha AD (mGy/MBq-s)': self.alpha_ad_mGy_per_MBq_s,
            'Electron AD (mGy/MBq-s)': self.electron_ad_mGy_per_MBq_s,
            'Alpha AD Total (Gy)': self.alpha_ad_Gy_total,
            'Electron AD Total (Gy)': self.electron_ad_Gy_total,
            'Full Beta Spectrum Used': self.full_beta_spectrum_used,
            'Interpolation Technique': self.interpolation_technique,
            'Unique ID': self.unique_id
        }
        return pd.DataFrame([data])

@dataclass
class AbsorbedDoseResult:

    source_tissue: str
    target_tissue: str
    site: str
    cellularity_factor: float
    interpolation_technique: str
    full_beta_spectrum_used: bool
    parent_nuclide_name: str

    results: list[AbsorbedDoseResultSingleNuclide]
    
    unique_id: str = ""

    def add_result(self, result: AbsorbedDoseResultSingleNuclide):
        self.results.append(result)

    def get_total_ad_alpha(self) -> float:
        total_ad = 0.0
        for res in self.results:
            total_ad += res.alpha_ad_Gy_total
        return total_ad

    def get_total_ad_electron(self) -> float:
        total_ad = 0.0
        for res in self.results:
            total_ad += res.electron_ad_Gy_total
        return total_ad
    
    def get_total_ad(self) -> float:
        total_ad = 0.0
        for res in self.results:
            total_ad += res.electron_ad_Gy_total + res.alpha_ad_Gy_total

        return total_ad

    def get_dose_factor(self) -> float:
        dose_factor = 0.0
        for res in self.results:
            dose_factor += res.electron_ad_mGy_per_MBq_s + res.alpha_ad_mGy_per_MBq_s

        return dose_factor
    
    def convert_to_dataframe(self) -> pd.DataFrame:
        df_list = [res.convert_to_dataframe() for res in self.results]
        return pd.concat(df_list, ignore_index=True)

    def combine_and_make_df(self) -> pd.DataFrame:
        """Combine all absorbed dose results into a single DataFrame."""

        # Get the other data from the results:

        for res in self.results:
            assert res.source_tissue == self.source_tissue, "Source tissue mismatch in results"
            assert res.target_tissue == self.target_tissue, "Target tissue mismatch in results"
            assert res.site == self.site, "Site mismatch in results"
            assert res.cellularity_factor == self.cellularity_factor, "Cellularity factor mismatch in results"
            assert res.interpolation_technique == self.interpolation_technique, "Interpolation technique mismatch in results"
            assert res.full_beta_spectrum_used == self.full_beta_spectrum_used, "Full beta spectrum used mismatch in results"

        combined_df = pd.DataFrame()

        combined_df['Dose Factor'] = [self.get_dose_factor()]
        combined_df['Total AD (Alpha)'] = [self.get_total_ad_alpha()]
        combined_df['Total AD (Electron)'] = [self.get_total_ad_electron()]
        combined_df['Total AD'] = [self.get_total_ad()]

        combined_df['Source Tissue'] = [self.source_tissue]
        combined_df['Target Tissue'] = [self.target_tissue]
        combined_df['Site'] = [self.site]
        combined_df['Cellularity Factor'] = [self.cellularity_factor]
        combined_df['Interpolation Technique'] = [self.interpolation_technique]
        combined_df['Full Beta Spectrum Used'] = [self.full_beta_spectrum_used]

        combined_df['Parent Nuclide'] = [self.parent_nuclide_name]


        return combined_df

"""
The source_target_pair dataclass will hold all the information needed
to calculate the absorbed dose to red marrow from a given source region in
a specific target region.
"""
@dataclass
class Source_Target_Pair:
    source: str
    target: str
    electron_phi_table: pd.DataFrame
    alpha_phi_table: pd.DataFrame
    target_mass: float
    site_trabecular_volume: float

    site_name: str

"""
We need a class to contain the different pairs for a specific skeletal site
so we can easily access them later and debug during development.
"""
@dataclass
class SkeletalSite:
    site_name: str
    source_target_pairs: list[Source_Target_Pair]

    def get_source_target_pair(self, source: str, target: str) -> Source_Target_Pair:
        for pair in self.source_target_pairs:
            if pair.source == source and pair.target == target:
                return pair
        return None
    
    def save_to_site_to_disk(self, filename: str):
        """Save the skeletal site data to disk."""
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
    @staticmethod
    def load_from_disk(filename: str) -> 'SkeletalSite':
        """Load the skeletal site data from disk."""
        with open(filename, 'rb') as f:
            site = pickle.load(f)
        return site

# Handle the mass of the source and target regions

RHO_RM = 1.03  # g/cm3 - From Hough and O'Reilly
J_from_MeV = 1.60218e-13  # Conversion factor from MeV to Joules

def create_red_marrow_source_target_pair(path_to_site_file: str, site_name: str,
                                         source_tissue: str,
                                         target_tissue: str,
                                         forced_mass = False,
                                         forced_reference_mass = None) -> Source_Target_Pair:
    """Create a source-target pair for red marrow in a specific site."""

    if forced_mass and forced_reference_mass is None:
        print("Forced mass is set to True, but no forced_reference_mass is provided")
        sys.exit()
    
    if forced_mass:
        mass_data = pd.read_excel(path_to_site_file, sheet_name='reference_masses')
        rm_mass_reference = forced_reference_mass
        print(f"Using forced reference mass for RM: {rm_mass_reference} g")
        rm_mass_ircp = forced_reference_mass*mass_data['ICRP_reference_CF']

    else:

        mass_data = pd.read_excel(path_to_site_file, sheet_name='reference_masses')

        rm_mass_reference = (mass_data['Homogenous_bone_volume']*mass_data['Spongiosa_fraction']*mass_data['Marrow_fraction'])[0]*RHO_RM
        rm_mass_ircp = rm_mass_reference*mass_data['ICRP_reference_CF']

        spongiosa_volume = mass_data['Homogenous_bone_volume'][0]*mass_data['Spongiosa_fraction'][0]

    print(f"Reference mass for RM in {site_name}: {rm_mass_reference} g")

    if target_tissue == "RM":

        alpha_sheet_name_by_source = {
            "TBE": 'Alpha_TBE_to_RM',
            "TBS": 'Alpha_TBS_to_RM',
            "TBV": 'Alpha_TBV_to_RM',
            "RM": 'Alpha_RM_to_RM'
        }
    else:
        print("Not implemented yet")
        sys.exit()
        
    CF_list = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    if source_tissue in ["TBE", "TBS", "TBV", "RM"]:

        alpha_source_sheet = alpha_sheet_name_by_source[source_tissue]

        print("Reading alpha source sheet:", alpha_source_sheet)

        Phi_table_alpha = pd.read_excel(path_to_site_file, sheet_name=alpha_source_sheet)

        # Make the alpha-table for phi

        phi_table_alpha = pd.DataFrame(columns=['Energy (MeV)', 'CF', 'Phi'])

        for index, row in Phi_table_alpha.iterrows():
            energy = row['Energy (MeV)']

            for cf in CF_list:
                phi = row[cf]

                if phi > 1.0:
                    phi = 1.0

                # Make a row for the phi table
                new_row = {'Energy (MeV)': energy, 'CF': cf, 'Phi': phi}
                phi_table_alpha = pd.concat([phi_table_alpha, pd.DataFrame([new_row])], ignore_index=True)
    

    # Make a phi table for the specific absorbed fractions calculated
    # from the Phi values in the table

    if target_tissue == "RM":
        print("Something")
        beta_sheet_name_by_source = {
            "RM": 'Beta_RM_to_RM_spesific',
            "TBS": 'Beta_TBS_to_RM_spesific',
            "TBV": 'Beta_TBV_to_RM_spesific',
        }
    else:
        print("Not implemented yet")
        sys.exit()

    if source_tissue == "RM":

        beta_source_sheet = beta_sheet_name_by_source[source_tissue]
        Phi_table = pd.read_excel(path_to_site_file, sheet_name=beta_source_sheet)

        phi_table_beta = pd.DataFrame(columns=['Energy (MeV)', 'CF', 'Phi'])

        # Now for each row in the table:

        for index, row in Phi_table.iterrows():
            energy = row['Energy (MeV)']

            for cf in CF_list:
                Phi = (row[cf])
                phi = rm_mass_reference*Phi*(cf/100)

                if phi > 1.0:
                    phi = 1.0

                # Make a row for the phi table
                new_row = {'Energy (MeV)': energy, 'CF': cf, 'Phi': phi}
                phi_table_beta = pd.concat([phi_table_beta, pd.DataFrame([new_row])], ignore_index=True)

    elif source_tissue in ["TBS", "TBV"]:
        print("Generating a CF-independent phi table for beta-emitters from Hough for source:", source_tissue)
        beta_source_sheet = beta_sheet_name_by_source[source_tissue]
        Phi_table = pd.read_excel(path_to_site_file, sheet_name=beta_source_sheet)

        phi_table_beta = pd.DataFrame(columns=['Energy (MeV)', 'Phi'])

        for index, row in Phi_table.iterrows():
            energy = row['Energy (MeV)']
            phi = row['Phi']*rm_mass_ircp[0]

            if phi > 1.0:
                phi = 1.0

            new_row = {'Energy (MeV)': energy, 'Phi': phi}
            phi_table_beta = pd.concat([phi_table_beta, pd.DataFrame([new_row])], ignore_index=True)

    else:
        print("Source tissue not implemented yet for beta-emitters, program exiting.")
        print("Not implemented yet")
        sys.exit()

    # I now have enough to construct the whole lumbar source-target pair for red marrow
    # as the target

    lumbar_red_marrow_to_red_marrow = Source_Target_Pair(
        source=source_tissue,
        target=target_tissue,
        electron_phi_table=phi_table_beta,
        alpha_phi_table=phi_table_alpha,
        target_mass=rm_mass_reference,
        site_name="Lumbar",
        site_trabecular_volume=spongiosa_volume
    )

    return lumbar_red_marrow_to_red_marrow

def plot_lumbar_red_marrow_phi_values():

    """
    Debug plot of the different phi values for lumbar red marrow
    from different source regions.
    """

    lumbar_TBS_to_RM = create_red_marrow_source_target_pair(path_to_site_file='ashen/resources/lumbar_fractions.xlsx',
                                             site_name="Lumbar", source_tissue="TBS", target_tissue="RM")

    lumbar_RM_to_RM = create_red_marrow_source_target_pair(path_to_site_file='ashen/resources/lumbar_fractions.xlsx',
                                             site_name="Lumbar", source_tissue="RM", target_tissue="RM")

    lumbar_TBV_to_RM = create_red_marrow_source_target_pair(path_to_site_file='ashen/resources/lumbar_fractions.xlsx',
                                             site_name="Lumbar", source_tissue="TBV", target_tissue="RM")


    beta_phi_TBS = lumbar_TBS_to_RM.electron_phi_table
    alpha_phi_TBS = lumbar_TBS_to_RM.alpha_phi_table

    beta_phi_rm = lumbar_RM_to_RM.electron_phi_table
    alpha_phi_rm = lumbar_RM_to_RM.alpha_phi_table

    beta_phi_TBV = lumbar_TBV_to_RM.electron_phi_table
    alpha_phi_TBV = lumbar_TBV_to_RM.alpha_phi_table

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(beta_phi_TBS['Energy (MeV)'], beta_phi_TBS['Phi'], "--o", label="TBS source")
    ax.set_xscale('log')
    ax.set_title("Beta phi values for TBS to RM")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for cf in [20, 40, 60, 80, 100]:
        single_cf = alpha_phi_TBS[alpha_phi_TBS['CF'] == cf]
        ax.plot(single_cf['Energy (MeV)'], single_cf['Phi'], "--o", label=f"CF={cf}")
    ax.set_title("Alpha phi values for TBS to RM")
    plt.legend()

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for cf in [20, 40, 60, 80, 100]:
        single_cf = alpha_phi_TBV[alpha_phi_TBV['CF'] == cf]
        ax.plot(single_cf['Energy (MeV)'], single_cf['Phi'], "--o", label=f"CF={cf}")

    ax.set_title("Alpha phi values for TBV to RM")
    plt.legend()

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for cf in [20, 40, 60, 80, 100]:
        single_cf = alpha_phi_rm[alpha_phi_rm['CF'] == cf]
        ax.plot(single_cf['Energy (MeV)'], single_cf['Phi'], "--o", label=f"CF={cf}")

    ax.set_title("Alpha phi values for RM to RM")
    plt.legend()

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for cf in [20, 40, 60, 80, 100]:
        single_cf = beta_phi_rm[beta_phi_rm['CF'] == cf]
        ax.plot(single_cf['Energy (MeV)'], single_cf['Phi'], "--o", label=f"CF={cf}")
    ax.set_title("Beta phi values for RM to RM")
    ax.set_xscale('log')
    plt.legend()


    plt.show()

def calculate_ad_from_nuclide_in_source_target_pair(nuclide: Nuclide,
                                              pair: Source_Target_Pair,
                                              source_name: str,
                                              CF: float,
                                              branch_fraction: float,
                                              cumulative_activity: float,
                                              interpolation_technique: str,
                                              uid_calc: str,
                                              use_full_beta_spectrum: bool,
                                              force_only_local_absorption: bool) -> float:
    """Calculate the absorbed dose to red marrow from a given nuclide in a source-target pair."""

    electron_phi_table = pair.electron_phi_table
    alpha_phi_table = pair.alpha_phi_table

    total_energy_absorbed = 0.0
    total_energy_emitted = 0.0

    total_alpha_energy_absorbed = 0.0
    total_electron_energy_absorbed = 0.0

    total_alpha_energy_emitted = 0.0
    total_electron_energy_emitted = 0.0

    for em in nuclide.emissions:
        particle_type = em.radiation_type

        if particle_type in ["B-", "IE", "AE"] and source_name == "RM":

            #closest_energy = electron_phi_table['Energy (MeV)'].iloc[(electron_phi_table['Energy (MeV)'] - em.energy).abs().argsort()[:1]].values[0]
            #phi_row = electron_phi_table[electron_phi_table['Energy (MeV)'] == closest_energy]
            single_cf_row = electron_phi_table[electron_phi_table['CF'] == CF]

            phi_interpolated = interpolate_phi(em.energy, single_cf_row, 
                                               interpolation_technique=interpolation_technique)
            
            if force_only_local_absorption:
                phi_interpolated = 1.0
            
            print(f"{em.energy} MeV -> interpolated phi: {phi_interpolated}")

            #assert single_cf_row.shape[0] == 1, "Expected exactly one row for the given energy and CF"

            #energy_absorbed = em.energy * em.yield_fraction * single_cf_row['Phi'].values[0]
            energy_absorbed = em.energy * em.yield_fraction * phi_interpolated
            
            energy_emitted = em.energy * em.yield_fraction

            total_electron_energy_absorbed += energy_absorbed
            total_electron_energy_emitted += energy_emitted

            total_energy_absorbed += energy_absorbed
            total_energy_emitted += energy_emitted

        elif particle_type in ["B-", "IE", "AE"] and source_name in ["TBS", "TBV"]:
            # Find the closes phi value in the table
            #closest_energy = electron_phi_table['Energy (MeV)'].iloc[(electron_phi_table['Energy (MeV)'] - em.energy).abs().argsort()[:1]].values[0]
            #phi_row = electron_phi_table[electron_phi_table['Energy (MeV)'] == closest_energy]
            single_cf_row = electron_phi_table

            #assert phi_row.shape[0] == 1, "Expected exactly one row for the given energy"

            phi_interpolated = interpolate_phi(em.energy, single_cf_row, interpolation_technique=interpolation_technique)

            if force_only_local_absorption:
                phi_interpolated = 1.0

            #energy_absorbed = em.energy * em.yield_fraction * phi_row['Phi'].values[0]
            energy_absorbed = em.energy * em.yield_fraction * phi_interpolated
            energy_emitted = em.energy * em.yield_fraction

            total_electron_energy_absorbed += energy_absorbed
            total_electron_energy_emitted += energy_emitted

            total_energy_absorbed += energy_absorbed
            total_energy_emitted += energy_emitted

        elif particle_type == "A" and (source_name == "RM" or source_name == "TBS" or source_name == "TBV"):

            # Find the closes phi value in the table
            #closest_energy = alpha_phi_table['Energy (MeV)'].iloc[(alpha_phi_table['Energy (MeV)'] - em.energy).abs().argsort()[:1]].values[0]
            #phi_row = alpha_phi_table[alpha_phi_table['Energy (MeV)'] == closest_energy]

            single_cf_row = alpha_phi_table[alpha_phi_table['CF'] == CF]
            phi_interpolated = interpolate_phi(em.energy, single_cf_row, interpolation_technique=interpolation_technique)

            if force_only_local_absorption:
                phi_interpolated = 1.0

            #assert single_cf_row.shape[0] == 1, "Expected exactly one row for the given energy and CF"

            energy_absorbed = em.energy * em.yield_fraction * phi_interpolated
            energy_emitted = em.energy * em.yield_fraction

            total_alpha_energy_absorbed += energy_absorbed
            total_alpha_energy_emitted += energy_emitted

            total_energy_absorbed += energy_absorbed
            total_energy_emitted += energy_emitted

            #print(f"Particle type: {particle_type}, Energy: {em.energy} MeV, Yield fraction: {em.yield_fraction}, phi: {single_cf_row['Phi'].values[0]}, Energy absorbed: {energy_absorbed} MeV")

        else:
            #print(f"Particle type {particle_type} not handled yet or {source_name}.")
            continue

    target_mass = pair.target_mass
    mass = target_mass*(CF/100)*1e-3  # Convert to kg

    total_energy_J = total_energy_absorbed*J_from_MeV*branch_fraction  # Convert MeV to Joules
    absorbed_dose_Gy = ((1e6*total_energy_J)/mass) # Convert to Gy for a MBq-s
    absorbed_dose_mGy_per_MBq_s = absorbed_dose_Gy*1e3

    print(f"Reference mass of RM target: {target_mass} g")
    print(f"Total absorbed energy in RM: {total_energy_absorbed} MeV")
    print(f"Mass of RM target: {mass} kg")
    print(f"Absorbed dose to RM: {absorbed_dose_mGy_per_MBq_s} mGy/MBq-s")
    # Print in scientific notation
    print(f"Absorbed dose to RM: {absorbed_dose_mGy_per_MBq_s:.3e} mGy/MBq-s")
    print(f"Total emitted energy in RM: {total_energy_emitted} MeV")
    print(f"Fraction of energy absorbed in RM: {total_energy_absorbed/total_energy_emitted:.3f}")

    # Do the alpha and electron separately

    total_alpha_energy_J = total_alpha_energy_absorbed*J_from_MeV*branch_fraction  # Convert MeV to Joules
    total_electron_energy_J = total_electron_energy_absorbed*J_from_MeV*branch_fraction  # Convert MeV to Joules

    absorbed_dose_alpha_Gy = ((1e6*total_alpha_energy_J)/mass) # Convert to Gy for a MBq-s
    absorbed_dose_electron_Gy = ((1e6*total_electron_energy_J)/mass) # Convert to Gy for a MBq-s

    absorbed_dose_alpha_mGy_per_MBq_s = absorbed_dose_alpha_Gy*1e3
    absorbed_dose_electron_mGy_per_MBq_s = absorbed_dose_electron_Gy*1e3

    absorbed_dose_alpha_Gy_total = absorbed_dose_alpha_Gy*cumulative_activity
    absorbed_dose_electron_Gy_total = absorbed_dose_electron_Gy*cumulative_activity


    result_ad_nuclide = AbsorbedDoseResultSingleNuclide(
        nuclide_name=nuclide.name,
        source_tissue=source_name,
        target_tissue="RM",
        site=pair.site_name,
        cellularity_factor=CF,
        alpha_ad_mGy_per_MBq_s=absorbed_dose_alpha_mGy_per_MBq_s,
        electron_ad_mGy_per_MBq_s=absorbed_dose_electron_mGy_per_MBq_s,

        alpha_ad_Gy_total=absorbed_dose_alpha_Gy_total,
        electron_ad_Gy_total=absorbed_dose_electron_Gy_total,
        full_beta_spectrum_used=use_full_beta_spectrum,
        interpolation_technique=interpolation_technique,
        unique_id=uid_calc

    )

    return result_ad_nuclide

def create_site_data():

    site_files  = {
        "Lumbar": 'ashen/resources/lumbar_fractions.xlsx',
    }

    source_types = ["RM", "TBS", "TBV"]

    site_output = {}


    for site_name, path_to_site_file in site_files.items():

        print(f"Creating source-target pairs for site: {site_name}")
        
        source_target_pairs = []

        for source in source_types:

            pair = create_red_marrow_source_target_pair(path_to_site_file=path_to_site_file,
                                                            site_name=site_name, source_tissue=source, target_tissue="RM")

            source_target_pairs.append(pair)

        site = SkeletalSite(
            site_name=site_name,
            source_target_pairs=source_target_pairs
        )

        site_output[site_name] = site

    return site_output

def calculate_absorbed_dose_to_red_marrow(
        input_site: str,
        input_cf: float,
        input_spongiosa_volume: float,
        input_cumulative_activity: float,
        source_name: str,
        target_name: str,
        nuclide_name: str,
        interpolation_technique: str = "loglog",
        full_beta_spectrum_used: bool = True,
        force_only_local_absorption: bool = False
    ) -> AbsorbedDoseResult:

    # Make a calculation UID from the UID library later

    uid_calc = f"{datetime.utcnow():%Y%m%dT%H%M%SZ}-{uuid.uuid4().hex[:8]}"

    CF = input_cf  # Cellularity Factor

    marrow_mass = input_spongiosa_volume*RHO_RM*(1-TB_FRACTIONS[input_site])  # g

    site = site_data[input_site]

    pair = site.get_source_target_pair(source=source_name, target=target_name)

    if source_name == "RM":

        input_rm_mass = marrow_mass*(CF/100)
        reference_mass = pair.target_mass

        cumulative_activity_corrected = (input_cumulative_activity/input_rm_mass)*reference_mass

    elif source_name in ["TBS"]:

        cumulative_activity_volume_conc = input_cumulative_activity/input_spongiosa_volume  # MBq-s/cm3
        cumulative_activity_corrected = cumulative_activity_volume_conc*pair.site_trabecular_volume

    else:
        print("Source tissue not implemented yet for cumulative activity correction, exiting.")
        sys.exit()

    print(f"Corrected cumulative activity in source region ({source_name}): {cumulative_activity_corrected} MBq-s")

    if RELOAD_DECAY_DB:

        emission_data = load_icrp_107()
        decay_db = make_decay_chain_db(emission_data = emission_data)
        with open('decay_db.pkl', 'wb') as f:
            pickle.dump(decay_db, f)

    with open('decay_db.pkl', 'rb') as f:
        loaded_decay_db = pickle.load(f)

    if USE_DEBUG_DB:
        loaded_decay_db = dummy_db

    parent_nuclide = loaded_decay_db.get_decay_info(nuclide_name)
    daughters = get_daughters(loaded_decay_db, nuclide_name)

    parent_dict = {parent_nuclide.name: 1.0}

    # Add the parent nuclide to the daughters list for dose calculation

    daughters[parent_nuclide.name] = 1.0

    # Remove stable nuclides from the daughters list

    for dau in list(daughters.keys()):
        if "$" in dau:
            del daughters[dau]

    ad_result = AbsorbedDoseResult(
        source_tissue=source_name,
        target_tissue=target_name,
        site=input_site,
        cellularity_factor=CF,
        results=[],
        interpolation_technique=interpolation_technique,
        full_beta_spectrum_used=full_beta_spectrum_used,
        unique_id=uid_calc,
        parent_nuclide_name=parent_nuclide.name
    )

    for nuclide_name in daughters.keys():
        print(f"Daughter nuclide: {nuclide_name}, Branching fraction: {daughters[nuclide_name]}")
        nuclide = loaded_decay_db.get_decay_info(nuclide_name)
        try:
            decay_type = extra_data_icrp107_lookup(nuclide.name)
        except ValueError:
            # Missing entry in extra ICRP 107 data (e.g., dummy/test nuclide)
            # Treat as no extra information available and continue gracefully
            print(f"extra ICRP 107 data not found for {nuclide.name}; assuming no beta emissions info")
            decay_type = ""

        if "B-" in decay_type:
            print(f"{nuclide.name} is a beta-emitter according to extra ICRP 107 data.")

            if full_beta_spectrum_used:
                nuclide = replace_simple_beta_with_full_spectrum(nuclide)

            else:
                print("Using simple beta spectrum for dose calculation. - Careful, this is not accurate!")

        ad_result_nuclide = calculate_ad_from_nuclide_in_source_target_pair(
            nuclide=nuclide,
            pair=pair,
            source_name=source_name,
            CF=CF,
            branch_fraction=daughters[nuclide.name],
            cumulative_activity=cumulative_activity_corrected,
            interpolation_technique=interpolation_technique,
            uid_calc=uid_calc,
            use_full_beta_spectrum=full_beta_spectrum_used,
            force_only_local_absorption=force_only_local_absorption
        )

        ad_result.add_result(ad_result_nuclide)

    print("Absorbed dose results for Lumbar RM from source:", source_name)
    for res in ad_result.results:
        print(f"Nuclide: {res.nuclide_name}, Alpha AD: {res.alpha_ad_mGy_per_MBq_s:.3e} mGy/MBq-s, Electron AD: {res.electron_ad_mGy_per_MBq_s:.3e} mGy/MBq-s")

    print(f"Total Alpha AD: {ad_result.get_total_ad_alpha():.3e} Gy")
    print(f"Total Electron AD: {ad_result.get_total_ad_electron():.3e} Gy")
    print(f"Total AD: {ad_result.get_total_ad():.3e} Gy")
    print(f"Dose factor: {ad_result.get_dose_factor():.3e} mGy/MBq-s")

    return ad_result

def is_chain_safe(db, parent_nuclide, daughter_absolute_warning_threshold_hours=0.25):

    daughters = get_daughters(db, parent_nuclide.name)

    # Check if there are only stable daughters

    unstable_daughters = []

    for dau_name in daughters.keys():
        if "$" in dau_name:
            continue
        else:
            unstable_daughters.append(dau_name)

    if len(unstable_daughters) == 0:
        print("All daughters are stable, chain is safe.")
        return True

    half_life_parent = parent_nuclide.halflife
    half_life_parent_hours = convert_to_hours(float(half_life_parent[0:-1]), half_life_parent[-1])

    daughter_half_lives = []

    for dau_name in daughters.keys():
        dau_nuclide = db.get_decay_info(dau_name)

        if dau_nuclide is None:
            if "$" in dau_name:
                continue
            else:
                print(f"Error: Daughter nuclide {dau_name} not found in database.")
                sys.exit()
        else:
            print(f"Checking daughter nuclide: {dau_name}")
            daughter_half_lives.append(convert_to_hours(float(dau_nuclide.halflife[0:-1]), dau_nuclide.halflife[-1]))


    print("Parent half-life (hours):", half_life_parent_hours)
    print("Daughter half-lives (hours):", daughter_half_lives)

    if half_life_parent_hours < min(daughter_half_lives, default=float('inf')):
        print("Chain is not safe: Parent half-life is shorter than daughter half-lives.")
        print(f"Parent: {parent_nuclide.name}, Half-life: {half_life_parent_hours} hours")
        longest_daughter = max(daughter_half_lives)
        print(f"Longest daughter half-life: {longest_daughter} hours")
        return False

    print("Chain is potentially safe.")
    longest_daughter = max(daughter_half_lives)
    print(f"Parent: {parent_nuclide.name}, Half-life: {half_life_parent_hours} hours")
    print(f"Longest daughter half-life: {longest_daughter} hours")

    if longest_daughter > daughter_absolute_warning_threshold_hours:
        print("Warning: Long-lived daughter detected in decay chain.")
        print(f"Longest daughter half-life: {longest_daughter} hours exceeds threshold of {daughter_absolute_warning_threshold_hours} hours. Consider using daughter migration")

    return True

def emission_profile(db, nuclide_name, np_photon_threshold=1e-3, do_plot = False,
                     plot_spectrum = False,
                     input_use_full_beta_spectra = True):

    nuclide = db.get_decay_info(nuclide_name)

    if nuclide is None:
        print(f"Error: Nuclide {nuclide_name} not found in database.")
        return None
    
    total_energy_be_type = {
        "Alpha": 0.0,
        "Beta-": 0.0,
        "Beta+": 0.0,
        "Gamma": 0.0,
        "X-ray": 0.0,
        "Auger": 0.0,
        "Internal Conversion": 0.0,
        "Alpha recoil": 0.0,
        "Spontaneous Fission": 0.0,
        "Prompt gamma": 0.0,
        "Delayed gamma": 0.0,
        "Neutron": 0.0
    }

    short_code_lookup = {
        "A": "Alpha",
        "B-": "Beta-",
        "B+": "Beta+",
        "G": "Gamma",
        "X": "X-ray",
        "AE": "Auger",
        "IE": "Internal Conversion",
        "AR": "Alpha recoil",
        "SF": "Spontaneous Fission",
        "PG": "Prompt gamma",
        "DG": "Delayed gamma",
        "N": "Neutron"
    }

    radiation_types_not_penetrating = [
        "Alpha", "Beta-", "Beta+", "Auger", "Internal Conversion",
    ]

    for em in nuclide.emissions:
        particle_type = em.radiation_type
        
        if particle_type not in short_code_lookup:
            print(f"Warning: Unknown radiation type {particle_type} for nuclide {nuclide_name}, skipping.")
            continue

        full_type = short_code_lookup[particle_type]

        energy_contribution = em.energy * em.yield_fraction

        total_energy_be_type[full_type] += energy_contribution

    # Tally up assumed non-penetrating radiation energy
    non_penetrating_energy = 0.0

    for rad_type in radiation_types_not_penetrating:
        non_penetrating_energy += total_energy_be_type[rad_type]

    print(f"Emission profile for nuclide: {nuclide_name}")
    print("Radiation Type\tTotal Energy Contribution (MeV)")

    for rad_type, total_energy in total_energy_be_type.items():
        if total_energy > 0: 
            print(f"{rad_type}\t{total_energy:.4f}")

    print(f"Fraction of energy per type:")
    total_energy_all = sum(total_energy_be_type.values())

    for rad_type, total_energy in total_energy_be_type.items():
        fraction = total_energy / total_energy_all if total_energy_all > 0 else 0
        if fraction > 0:
            print(f"{rad_type}\t{fraction:.4f}")

    print(f"Total non-penetrating energy: {non_penetrating_energy:.4f} MeV")
    print(f"Fraction of non-penetrating energy: {non_penetrating_energy / total_energy_all:.4f}")

    if do_plot:
        labels = []
        sizes = []
        for rad_type, total_energy in total_energy_be_type.items():
            if total_energy > 0:
                labels.append(rad_type)
                sizes.append(total_energy)

        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                shadow=True, startangle=90)
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

        plt.title(f"Emission profile for {nuclide_name}")
        plt.show()

    if plot_spectrum:

        # Lets also plot a spectrum of the emissions

        # First swap out with the full beta spectrum if needed
        if input_use_full_beta_spectra:
            nuclide = replace_simple_beta_with_full_spectrum(nuclide)
        else:
            print("Using simple beta spectrum for plotting. - Careful, this is not accurate!")

        # Combine beta, gamma and (auger + IC) spectra into a single 1x3 figure
        # Accept both legacy "B-" and verbose "Beta-" labels used in different parts
        beta_energies = [em.energy for em in nuclide.emissions if em.radiation_type in ("B-", "Beta-")]
        beta_yields = [em.yield_fraction for em in nuclide.emissions if em.radiation_type in ("B-", "Beta-")]

        gamma_energies = [em.energy for em in nuclide.emissions if em.radiation_type == "G"]
        gamma_yields = [em.yield_fraction for em in nuclide.emissions if em.radiation_type == "G"]

        auger_energies = [em.energy for em in nuclide.emissions if em.radiation_type == "AE"]
        auger_yields = [em.yield_fraction for em in nuclide.emissions if em.radiation_type == "AE"]

        ic_energies = [em.energy for em in nuclide.emissions if em.radiation_type == "IE"]
        ic_yields = [em.yield_fraction for em in nuclide.emissions if em.radiation_type == "IE"]

        fig, axs = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

        # Beta subplot (left)
        ax_beta = axs[0]
        if len(beta_energies) == 0:
            ax_beta.text(0.5, 0.5, 'No beta emissions', ha='center', va='center')
            ax_beta.set_xticks([])
            ax_beta.set_yticks([])
        else:
            ax_beta.stem(beta_energies, beta_yields, basefmt=" ", linefmt='C0-', markerfmt='C0o')
            ax_beta.set_xlabel("Energy (MeV)")
            ax_beta.set_ylabel("Yield Fraction")
        ax_beta.set_title(f"Beta Spectrum\n{nuclide_name}")

        # Gamma subplot (middle)
        ax_gamma = axs[1]
        if len(gamma_energies) == 0:
            ax_gamma.text(0.5, 0.5, 'No gamma emissions', ha='center', va='center')
            ax_gamma.set_xticks([])
            ax_gamma.set_yticks([])
        else:
            ax_gamma.stem(gamma_energies, gamma_yields, basefmt=" ", linefmt='C1-', markerfmt='C1s')
            ax_gamma.set_xlabel("Energy (MeV)")
            ax_gamma.set_ylabel("Yield Fraction")
        ax_gamma.set_title(f"Gamma Spectrum\n{nuclide_name}")

        # Auger + IC subplot (right)
        ax_aug = axs[2]
        if len(auger_energies) == 0 and len(ic_energies) == 0:
            ax_aug.text(0.5, 0.5, 'No Auger/IC emissions', ha='center', va='center')
            ax_aug.set_xticks([])
            ax_aug.set_yticks([])
        else:
            if len(auger_energies) > 0:
                ax_aug.stem(auger_energies, auger_yields, basefmt=" ", linefmt='C2-', markerfmt='C2^', label='Auger')
            if len(ic_energies) > 0:
                ax_aug.stem(ic_energies, ic_yields, basefmt=" ", linefmt='r-', markerfmt='ro', label='IC')
            ax_aug.set_xlabel("Energy (MeV)")
            ax_aug.set_ylabel("Yield Fraction")
            ax_aug.legend()
        ax_aug.set_title(f"Auger / IC\n{nuclide_name}")

        plt.suptitle(f"Emission Spectra for {nuclide_name}")
        plt.show()


    return None

def interpolate_phi(energy: float, 
                    single_cf_table: pd.DataFrame,
                    interpolation_technique: str) -> float:

    E = single_cf_table['Energy (MeV)'].values
    Phi = single_cf_table['Phi'].values

    if single_cf_table.shape[0] == 0:
        print(f"No phi data found")
        return None
    
    if interpolation_technique == "closest":
        closest_energy = single_cf_table['Energy (MeV)'].iloc[(single_cf_table['Energy (MeV)'] - energy).abs().argsort()[:1]].values[0]
        phi_value = single_cf_table[single_cf_table['Energy (MeV)'] == closest_energy]['Phi'].values[0]
        return phi_value
    
    elif interpolation_technique == "linear":
        idx = np.abs(E - energy).argmin()
        return float(Phi[idx])
    
    elif interpolation_technique == "spline":
        f = interp1d(E, Phi, kind='cubic', bounds_error=False, fill_value="extrapolate")
        return float(f(energy))
    
    elif interpolation_technique == "loglog":

        return float(np.exp(np.interp(np.log(energy), np.log(E), np.log(Phi))))

def batch_calculation_from_input_file(input_data_path: str,
                                      output_data_path: str,
                                      decay_chain_path: str):


    aggregated_df = pd.DataFrame()
    aggregated_chain_df = pd.DataFrame()

    input_df = pd.read_excel(input_data_path)

    print(input_df)

    # Rewrite with tqdm now with progress bar

    for row in tqdm(input_df.itertuples(index=False), desc="Processing input rows"):
        input_row = row
        print("Processing input row:", input_row)

        input_site = input_row.Site
        input_cf = input_row.CF
        input_spongiosa_volume = input_row.Spongiosa_Volume_cm3
        input_cumulative_activity = input_row.Cumulative_Activity_MBq_s
        source_name = input_row.Source
        target_name = input_row.Target
        nuclide_name = input_row.Nuclide
        input_interpolation_technique = input_row.Interpolation_Technique
        input_use_full_beta_spectra = bool(input_row.Use_Full_Beta)



        lumbar_site = site_data[input_site]

        result = calculate_absorbed_dose_to_red_marrow(
            input_site=input_site,
            input_cf=input_cf,
            input_spongiosa_volume=input_spongiosa_volume,
            input_cumulative_activity=input_cumulative_activity,
            source_name=source_name,
            target_name=target_name,
            nuclide_name=nuclide_name,
            interpolation_technique=input_interpolation_technique,
            full_beta_spectrum_used=input_use_full_beta_spectra
        )

        # Convert the result to a dataframe

        result_df = result.convert_to_dataframe()

        chain_df = result.combine_and_make_df()

        # Always append the chain results. Previously this only appended when
        # `aggregated_df` already had rows which meant for a single-input-row run
        # the chain dataframe was never added and remained empty. Concatenate
        # regardless; if `chain_df` is empty this is a no-op.
        aggregated_chain_df = pd.concat([aggregated_chain_df, chain_df], ignore_index=True)

        aggregated_df = pd.concat([aggregated_df, result_df], ignore_index=True)


        # Append to output Excel file if it exists, otherwise create a new file

        #if os.path.exists(output_data_path):
        #    with pd.ExcelWriter(output_data_path, mode='a', if_sheet_exists='replace') as writer:
        #        result_df.to_excel(writer, sheet_name='Results', index=False)
        #else:
        #    result_df.to_excel(output_data_path, sheet_name='Results', index=False)

        #print(result_df)

        #dose_factor = result.get_dose_factor()

    print("Writing aggregated results to output file:", output_data_path)
    print(aggregated_df)

    aggregated_df.to_excel(output_data_path, sheet_name='Results', index=False)
    print("Writing aggregated decay chain results to output file:", decay_chain_path)
    print(aggregated_chain_df)
    aggregated_chain_df.to_excel(decay_chain_path, sheet_name='Decay_Chain_Results', index=False)


if __name__ == "__main__":

    site_data = create_site_data()

    # Test of the chain safe function 

    if RELOAD_DECAY_DB:

        emission_data = load_icrp_107()
        decay_db = make_decay_chain_db(emission_data = emission_data)
        with open('decay_db.pkl', 'wb') as f:
            pickle.dump(decay_db, f)

    with open('decay_db.pkl', 'rb') as f:
        loaded_decay_db = pickle.load(f)

    #is_chain_safe(loaded_decay_db, loaded_decay_db.get_decay_info("Tb-149"))
    #emission_profile(loaded_decay_db, "Tb-149", plot_spectrum=True)

    input_site = "Lumbar"
    input_cf = 70
    input_spongiosa_volume = 210.261312
    input_cumulative_activity = 1
    source_name = "RM"
    target_name = "RM"
    nuclide_name = "Ac-225"
    input_interpolation_technique = "loglog"
    input_use_full_beta_spectra = True
    force_only_local_absorption = False

    result = calculate_absorbed_dose_to_red_marrow(
        input_site=input_site,
        input_cf=input_cf,
        input_spongiosa_volume=input_spongiosa_volume,
        input_cumulative_activity=input_cumulative_activity,
        source_name=source_name,
        target_name=target_name,
        nuclide_name=nuclide_name,
        interpolation_technique=input_interpolation_technique,
        full_beta_spectrum_used=input_use_full_beta_spectra,
        force_only_local_absorption=force_only_local_absorption # For debugging purposes
    )

    dose_factor = result.get_dose_factor()

    print(dose_factor)

    # Test for a file input

    #input_data_path = "ashen/resources/abstract_data.xlsx"
    #output_data_path = "ashen/red_marrow_output/tat13_abstract_data_absorbed_dose_results.xlsx"
    #decay_chain_path = "ashen/red_marrow_output/tat13_abstract_data_combined_decay_chain_results.xlsx"
    #batch_calculation_from_input_file(input_data_path, output_data_path, decay_chain_path)


    # Find the right benchmark value

#    jh_cf = JH_dose_factors[JH_dose_factors['CF'] == input_cf]
#    jh_nuclide = jh_cf[jh_cf['Nuclide'] == nuclide_name]
#    jh_site = jh_nuclide[jh_nuclide['Site'] == input_site]
#    jh_target = jh_site[jh_site['Target'] == target_name]
#    jh_source = jh_target[jh_target['Source'] == source_name]
#
#    if jh_source.empty:
#        print("No benchmark data found for the given parameters.")
#
#    else:
#        benchmark_dose_factor = jh_source['DoseFactor'].values[0]
#
#        # Check the result against the benchmark
#
#        print(f"Calculated dose factor: {dose_factor:.6e} mGy/MBq-s")
#        print(f"Benchmark dose factor: {benchmark_dose_factor:.6e} mGy/MBq-s")
#        print(f"Percent difference: {100.0*(dose_factor - benchmark_dose_factor)/benchmark_dose_factor:.3f} %")
#
#
#
#
#    print("Final absorbed dose result:", result)
#
    print("Have a nice day")



    #input_site = "Lumbar" # Dropdown
    #input_cf = 70 # Dropdown
    #input_interpolation_technique = "loglog"  # Dropdown

    #source_name = "RM" # Dropdown
    #target_name = "RM" # Dropdown

    #nuclide_name = 'Ac-225' # Text input or Dropdown

    #input_spongiosa_volume = 210.261312# cm3 # Numeric input
    #input_cumulative_activity = 1 # MBq-s # Numeric input

    #alpha_or_beta = "Beta"  # Dropdown

    #if alpha_or_beta == "Beta":

    #    phi_table = site_data[input_site].get_source_target_pair(source=source_name, target=target_name).electron_phi_table
    #    energies = np.arange(0.01, 3.01, 0.01)  # MeV

    #elif alpha_or_beta == "Alpha":
    #    phi_table = site_data[input_site].get_source_target_pair(source=source_name, target=target_name).alpha_phi_table
    #    
    #    energies = np.arange(3.01, 10.0, 0.1)  # MeV

    #if "CF" in phi_table.columns:
    #    Phi_table_parsed = phi_table[phi_table['CF'] == input_cf]
    #else:
    #    Phi_table_parsed = phi_table


    #interpolation_methods = ["closest", "linear", "spline",
    #                         "loglog"]


    #result_vector = np.zeros((len(energies), len(interpolation_methods)))

    #for i, energy in enumerate(energies):
    #    for j, method in enumerate(interpolation_methods):
    #        phi_value = interpolate_phi(energy=energy, single_cf_table=Phi_table_parsed,
    #                                    interpolation_technique=method)
    #        result_vector[i, j] = phi_value
    #        print(f"Energy: {energy} MeV, Method: {method}, Phi: {phi_value}")

    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1)
    #for j, method in enumerate(interpolation_methods):
    #    ax.plot(energies, result_vector[:, j], "--o", label=method)

    #ax.set_xscale('log')
    #ax.set_title(f"Interpolation methods for {alpha_or_beta} phi values")
    #plt.legend()
    #plt.show()
