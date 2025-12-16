"""
One-off script to check, plot and aggregate red marrow dose results.

"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def generate_input_file_collection(path_to_input_files: str,
                                   list_of_parent_nuclides: list,
                                   interpolation_method_list: str,
                                   source_tissue_list: list,
                                   target_tissue_list: list,
                                   site_list: list,
                                   CF_list: list,
                                   beta_spectrum_list: list,
                                   spongiosa_list: list,
                                   Cumulative_Act_list: list) -> None:
    
    input_df = pd.DataFrame(columns=[
        "Parent",
        "Interpolation_Technique",
        "Source",
        "Target",
        "Site",
        "CF",
        "Use_Full_Beta",
        "Spongiosa_Volume_cm3",
        "Cumulative_Activity_MBq_s"
    ])

    # Make the combinations

    for parent in list_of_parent_nuclides:
        for interp in interpolation_method_list:
            for source in source_tissue_list:
                for target in target_tissue_list:
                    for site in site_list:
                        for CF in CF_list:
                            for use_full_beta in beta_spectrum_list:
                                for spongiosa_vol in spongiosa_list:
                                    for cum_act in Cumulative_Act_list:
                                        new_row = {
                                            "Nuclide": parent,
                                            "Interpolation_Technique": interp,
                                            "Source": source,
                                            "Target": target,
                                            "Site": site,
                                            "CF": CF,
                                            "Use_Full_Beta": use_full_beta,
                                            "Spongiosa_Volume_cm3": spongiosa_vol,
                                            "Cumulative_Activity_MBq_s": cum_act
                                        }
                                        input_df = input_df._append(new_row, ignore_index=True)

    # Save to Excel
    input_df.to_excel(path_to_input_files, index=False)

    return input_df


nuclide_list = ["Lu-177", 
                "Y-90", 
                "At-211", 
                "Ra-224", 
                "Ra-223", 
                "Ac-225", 
                "Th-227",
                "Tb-161"]

interpolation_techniques = ["linear", "loglog", "closest", "spline"]
source_tissues = ["RM", "TBS"]
use_full_beta_spectra_options = [True]
cf_values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
spongiosa_volumes = [200]
cumulative_activities = [1]

generate_input_file_collection(
    path_to_input_files="ashen/resources/abstract_data.xlsx",
    list_of_parent_nuclides=nuclide_list,
    interpolation_method_list=interpolation_techniques,
    source_tissue_list=source_tissues,
    target_tissue_list=["RM"],
    site_list=["Lumbar"],
    CF_list=cf_values,
    beta_spectrum_list=use_full_beta_spectra_options,
    spongiosa_list=spongiosa_volumes,
    Cumulative_Act_list=cumulative_activities
)


sys.exit()

input_data = "ashen/red_marrow_output/test_1_combined_decay_chain_results.xlsx"

data_df = pd.read_excel(input_data, sheet_name="Decay_Chain_Results")

# Check the contents of the DataFrame
print(data_df.head())

df_only_full_beta = data_df

print(f"Number of results with full beta spectrum used: {len(df_only_full_beta)}")

parent_nuclides = np.unique(df_only_full_beta['Parent Nuclide'])

for parent in parent_nuclides:
    df_parent = df_only_full_beta[df_only_full_beta['Parent Nuclide'] == parent]

    if len(df_parent) != 8:
        print(f"Warning: Expected 8 source-target pairs for parent nuclide {parent}, found {len(df_parent)}")
        continue

    # Calculate the mean and standard deviation of absorbed doses

    mean_ad = df_parent['Dose Factor'].mean()
    std_ad = df_parent['Dose Factor'].std()
    cov = 100*(std_ad / mean_ad) if mean_ad != 0 else 0
    print(f"Parent Nuclide: {parent} - CoV of Dose Factor: {cov:.4f}")


    