"""
Script to save predefined values and simulated data to JSON files for easy access.

This would be the radionuclides, sites, and other predefined parameters used in the GUI.

"""

import json
import sys
import csv
import yaml
from os.path import exists
from collections import defaultdict
from pathlib import Path

def deep_merge(dst, src):
    for key, value in src.items():
        if (
            key in dst
            and isinstance(dst[key], dict)
            and isinstance(value, dict)
        ):
            deep_merge(dst[key], value)
        else:
            dst[key] = value

from ashen.ashen_utils import (
    make_decay_chain_db,
    load_icrp_107,
    energy_in_decay_chain,

)

from ashen.red_marrow_module  import (
    return_list_of_sites_by_nuclide,
)

from typing import Dict
import numpy as np
import pandas as pd
from pathlib import Path

def make_alpha_saf_yaml():

    csv_input_files = {
            "cervical_vertebrae": "resources/SAFs/ALPHAS/cervical_vertebrae",
            "femur_head": "resources/SAFs/ALPHAS/femur_head",
            "femur_neck": "resources/SAFs/ALPHAS/femur_neck",
            "iliac_crest": "resources/SAFs/ALPHAS/iliac_crest",
            "lumbar_vertebrae": "resources/SAFs/ALPHAS/lumbar_vertebrae",
            "ribs": "resources/SAFs/ALPHAS/ribs",
            "parietal_bone": "resources/SAFs/ALPHAS/parietal_bone",
    }

    yaml_output_files = {
            "cervical_vertebrae": "resources/SAFs/ALPHAS/cervical_vertebrae",
            "femur_head": "resources/SAFs/ALPHAS/femur_head",
            "femur_neck": "resources/SAFs/ALPHAS/femur_neck",
            "iliac_crest": "resources/SAFs/ALPHAS/iliac_crest",
            "lumbar_vertebrae": "resources/SAFs/ALPHAS/lumbar_vertebrae",
            "ribs": "resources/SAFs/ALPHAS/ribs",
            "parietal_bone": "resources/SAFs/ALPHAS/parietal_bone",
    }

    sites = csv_input_files.keys()

    source_tissues = ["rm", "tbs"]

    source_tissue_dict = {
        "rm": "red_marrow",
        "tbs": "tbs",
    }

    full_yaml_output_file_list = []

    for source_tissue in source_tissues:

        for site in sites:

            print(f"Processing site: {site}")

            csv_file = csv_input_files[site] + f"_{source_tissue}.csv"
            output_yaml = yaml_output_files[site] + f"_{source_tissue}.yaml"

            # Check if csv-file exists

            if not exists(csv_file):
                print(f"CSV file {csv_file} does not exist. Skipping.")
                continue

            site_name = site
            emission_type = "alphas"
            cfs = ["100", "90", "80", "70", "60", "50", "40", "30", "20", "10"]

            # ---------- READ CSV ----------
            df = pd.read_csv(csv_file, sep=";")

            # ---------- INIT NESTED DICT ----------
            saf_dict = {site_name: {emission_type: {}}}

            # ---------- Tissue ----------

            saf_dict[site_name][emission_type][source_tissue_dict[source_tissue]] = {}
            for cf in cfs:
                saf_dict[site_name][emission_type][source_tissue_dict[source_tissue]][cf] = {}
                for _, row in df.iterrows():
                    energy = float(row["energy"])
                    value = float(row[cf])
                    saf_dict[site_name][emission_type][source_tissue_dict[source_tissue]][cf][energy] = value
            # --------- WRITE YAML ----------
            full_yaml_output_file_list.append(output_yaml)
            
            with open(output_yaml, "w") as f:
                yaml.dump(saf_dict, f, sort_keys=False)

    # ---------- SETTINGS ----------

    list_of_yaml_files = full_yaml_output_file_list

    output_yaml = "combined_alpha_saf.yaml"

    # ---------- LOAD AND MERGE ----------

    skeletal_data = {}

    for fname in list_of_yaml_files:
        print(f"File: {fname}")
        with open(fname) as f:
            data = yaml.safe_load(f)
            deep_merge(skeletal_data, data)

    # ---------- WRITE COMBINED YAML ----------

    with open(output_yaml, "w") as f:
        yaml.safe_dump(skeletal_data, f, sort_keys=False)

def make_electron_saf_yaml():

    csv_input_files = {
            "craniofacial_bones": "resources/SAFs/ELECTRONS/craniofacial_bones.csv",
            "mandible": "resources/SAFs/ELECTRONS/mandibles.csv",
            "scapulae": "resources/SAFs/ELECTRONS/scapulae.csv",
            "clavicles": "resources/SAFs/ELECTRONS/clavicles.csv",
            "sternum": "resources/SAFs/ELECTRONS/sternum.csv",
            "ribs": "resources/SAFs/ELECTRONS/ribs.csv",
            "cervical_vertebrae": "resources/SAFs/ELECTRONS/cervical_vertebrae.csv",
            "thoracic_vertebrae": "resources/SAFs/ELECTRONS/thoracic_vertebrae.csv",
            "lumbar_vertebrae": "resources/SAFs/ELECTRONS/lumbar_vertebrae.csv",
            "sacrum": "resources/SAFs/ELECTRONS/sacrum.csv",
            "os_coxae": "resources/SAFs/ELECTRONS/os_coxae.csv",
            "proximal_humeri": "resources/SAFs/ELECTRONS/proximal_humeri.csv",
            "proximal_femora": "resources/SAFs/ELECTRONS/proximal_femora.csv",
    }

    yaml_output_files = {
            "craniofacial_bones": "resources/SAFs/ELECTRONS/craniofacial_bones.yaml",
            "mandible": "resources/SAFs/ELECTRONS/mandible.yaml",
            "scapulae": "resources/SAFs/ELECTRONS/scapulae.yaml",
            "clavicles": "resources/SAFs/ELECTRONS/clavicles.yaml",
            "sternum": "resources/SAFs/ELECTRONS/sternum.yaml",
            "ribs": "resources/SAFs/ELECTRONS/ribs.yaml",
            "cervical_vertebrae": "resources/SAFs/ELECTRONS/cervical_vertebrae.yaml",
            "thoracic_vertebrae": "resources/SAFs/ELECTRONS/thoracic_vertebrae.yaml",
            "lumbar_vertebrae": "resources/SAFs/ELECTRONS/lumbar_vertebrae.yaml",
            "sacrum": "resources/SAFs/ELECTRONS/sacrum.yaml",
            "os_coxae": "resources/SAFs/ELECTRONS/os_coxae.yaml",
            "proximal_humeri": "resources/SAFs/ELECTRONS/proximal_humeri.yaml",
            "proximal_femora": "resources/SAFs/ELECTRONS/proximal_femora.yaml",
    }

    sites = csv_input_files.keys()

    for site in sites:

        print(f"Processing site: {site}")

        csv_file = csv_input_files[site]
        output_yaml = yaml_output_files[site]
        site_name = site
        emission_type = "electrons"

        red_marrow_cfs = ["100", "90", "80", "70", "60", "50", "40", "30", "20", "10", "icrp"]
        other_tissue = ["tbs", "tbv"]

        # ---------- READ CSV ----------
        df = pd.read_csv(csv_file, sep=";")

        # ---------- INIT NESTED DICT ----------
        saf_dict = {site_name: {emission_type: {}}}

        # ---------- RED MARROW ----------
        saf_dict[site_name][emission_type]["red_marrow"] = {}
        for cf in red_marrow_cfs:
            saf_dict[site_name][emission_type]["red_marrow"][cf] = {}
            for _, row in df.iterrows():
                energy = float(row["energy"])
                value = float(row[cf])
                saf_dict[site_name][emission_type]["red_marrow"][cf][energy] = value

        # ---------- OTHER SOURCE TISSUES ----------
        for tissue in other_tissue:
            saf_dict[site_name][emission_type][tissue] = {0: {}}
            for _, row in df.iterrows():
                energy = float(row["energy"])
                value = float(row[tissue])
                saf_dict[site_name][emission_type][tissue][0][energy] = value

        print(saf_dict)

        # ---------- WRITE YAML ----------
        with open(output_yaml, "w") as f:
            yaml.dump(saf_dict, f, sort_keys=False)

    # ---------- SETTINGS ----------

    list_of_yaml_files = yaml_output_files.values()

    output_yaml = "combined_saf.yaml"

    # ---------- LOAD AND MERGE ----------
    combined_dict = {}

    for yaml_file in list_of_yaml_files:
        with open(yaml_file, "r") as f:
            site_dict = yaml.safe_load(f)
            # Merge into combined dict
            combined_dict.update(site_dict)  # site_name is top-level key in each file

    # ---------- WRITE COMBINED YAML ----------
    with open(output_yaml, "w") as f:
        yaml.dump(combined_dict, f, sort_keys=False)

    print(f"Combined YAML saved to {output_yaml}")

# Load decay chains and emission data
emission_energy: Dict[str, float] = load_icrp_107()
db = make_decay_chain_db(emission_data=emission_energy)

# Save the daughters to a JSON file in the dict "CHECKBOX_VALUES"

CHECKBOX_VALUES = {}
CHECKBOX_TITLES = {}

SEVEN_FIELD_VALUES = ["Ac-225", "Pb-212", "At-211"]
PREDEFINED_VALUES = ["Lu-177", "Y-90", "I-131"]
PREDEFINED_VALUES += list(SEVEN_FIELD_VALUES)

# Find a list of alpha and beta-emitters

nuc_list = list(db.nuclides.keys())

for nuc in nuc_list:
    nuclide = db.get_decay_info(nuc)

    if nuclide.is_alpha:
        print(f"Alpha emitter found: {nuc}")
        SEVEN_FIELD_VALUES.append(nuc)

    if nuclide.is_electron_emitter:
        print(f"Electron emitter found: {nuc}")
        if nuc not in PREDEFINED_VALUES:
            PREDEFINED_VALUES.append(nuc)

PREDEFINED_VALUES += list(SEVEN_FIELD_VALUES)

for nuclide in PREDEFINED_VALUES:

    daughters = db.get_half_lives_of_daughters(nuclide)
    list_of_daughters = list(daughters.keys())
    CHECKBOX_TITLES[nuclide] = list_of_daughters

    for da in daughters:
        CHECKBOX_VALUES[da] = daughters[da]


# Write both dicts into a single JSON object
with open("predefined_data.json", "w") as f:
    json.dump({
        "CHECKBOX_VALUES": CHECKBOX_VALUES,
        "CHECKBOX_TITLES": CHECKBOX_TITLES,
        "SEVEN_FIELD_VALUES": SEVEN_FIELD_VALUES,
        "PREDEFINED_VALUES": PREDEFINED_VALUES
    }, f, indent=4)

make_alpha_saf_yaml()

sys.exit()

alpha_emitter_sites = [
        "cervical_vertebrae",
        "femur_head",
        "femur_neck",
        "iliac_crest",
        "lumbar_vertebrae",
        "ribs",
        "parietal_bone"
]

electron_emitter_sites = [
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
