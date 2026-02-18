"""

Functions to calculate the yield from beta-emitting radionuclides.

This is needed as the group in Gothenburg have used full spectrum beta data
and we need to be able to compare our results to theirs.

Update when moving myelodose out of the ashen package
Third update, because I am an idiot


"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import pandas as pd
import sys
from ashen.ashen_utils import load_icrp_107, make_decay_chain_db, get_daughters, Nuclide, RadiationEmission
import pickle
from pathlib import Path

PATH_TO_BETA_FILES = Path(__file__).parent / 'resources' / 'ICRP-07.BET'

def extract_data(filename, keyword) -> pd.DataFrame:
    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith(keyword):
            # Split the line to get the number after the keyword
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"No number found after keyword {keyword}")
            n_entries = int(parts[1])
            # Extract the following n_entries lines
            data_lines = lines[i+1 : i+1+n_entries]
            # Convert to floats if desired
            data = [tuple(map(float, l.split())) for l in data_lines]
            df_spectrum = pd.DataFrame(data, columns=['E(MeV)', 'P(E)'])
            return df_spectrum

    raise ValueError(f"Keyword {keyword} not found in file")


def compute_yields(df):
    E = df['E(MeV)'].values
    P = df['P(E)'].values.astype(float)
    # trapezoidal integration for P dE and E*P dE
    A = np.trapz(P, E)                 # total area under P(E)
    Eweighted = np.trapz(E*P, E)      # ∫ E*P(E) dE
    meanE = Eweighted / A if A>0 else np.nan

    # build bin list: centre energy, fractional yield per bin
    energies = 0.5*(E[:-1] + E[1:])
    binareas = 0.5*(P[:-1] + P[1:]) * (E[1:] - E[:-1])
    fractions = binareas / A if A>0 else np.zeros_like(binareas)

    out = pd.DataFrame({
        'Energy_MeV': energies,
        'Yield': fractions,
        'BinArea': binareas
    })
    return {'A': A, 'meanE_MeV': meanE, 'meanE_J': meanE*1.602176634e-13 if not np.isnan(meanE) else np.nan,
            'bins': out}

def replace_simple_beta_with_full_spectrum(nuclide: Nuclide) -> Nuclide:
    # Remove existing beta emissions
    nuclide.emissions = [em for em in nuclide.emissions if em.radiation_type != "B-"]

    beta_spectrum_data = extract_data(PATH_TO_BETA_FILES, nuclide.name)
    beta_yields = compute_yields(beta_spectrum_data)
    for _, row in beta_yields['bins'].iterrows():
        energy = row['Energy_MeV']
        yield_fraction = row['Yield']
        if yield_fraction > 0:
            new_emission = RadiationEmission(
                radiation_type="B-",
                energy=energy,
                yield_fraction=yield_fraction
            )
            nuclide.emissions.append(new_emission)

    # Print some sanity check info

    total_yield = sum(em.yield_fraction for em in nuclide.emissions if em.radiation_type == "B-")
    mean_energy = sum(em.energy * em.yield_fraction for em in nuclide.emissions if em.radiation_type == "B-") / total_yield if total_yield > 0 else 0
    print(f"Nuclide: {nuclide.name} - Total beta yield: {total_yield:.6f} (should be 1.0), Mean beta energy: {mean_energy:.6f} MeV, Expected mean energy: {beta_yields['meanE_MeV']:.6f} MeV")

    return nuclide

if __name__ == "__main__":

    # Load the beta spectrum data
    beta_spectrum_data = extract_data(path_to_beta_files, 'Lu-177')

    df_spectrum = pd.DataFrame(beta_spectrum_data, columns=['E(MeV)', 'P(E)'])
    results = compute_yields(df_spectrum)

    # Print results
    print(f"Total yield (should be 1.0): {results['A']:.6f} per decay")
    print(f"Mean beta energy: {results['meanE_MeV']:.6f} MeV ({results['meanE_J']:.6e} J)")

    print(np.sum(results['bins']['Yield']))
    print(np.sum(results['bins']['BinArea']))
    print(np.sum(results['bins']['Yield'] * results['bins']['Energy_MeV']))



    plt.show()

