"""
This is a demontration script
to show a simple use case of Ashen

Say that you want to retreive the
energy radiated from a couple of radionuclides.

Some simple, like 177-Lu, and some more complex, like
224-Ra

"""

from ashen.ashen_utils import (
    make_decay_chain_db,
    load_icrp_107,
    energy_in_decay_chain,
    get_daughters,
)

# Usecase 1: Simple decay chain
# Load decay chains and emission data
emission_energy = load_icrp_107()
db = make_decay_chain_db(emission_data=emission_energy)

E = energy_in_decay_chain(
    db=db,
    nuc_name="Lu-177",
)

# This will return the energy of the decay chain of Lu-177
# since it is a simple decay, it will return the energy of the decay
# chain of Lu-177, which is roughly 0.133 MeV for the beta emission
# and 0.147 MeV for the combined non penetrating particle and gamma emissions
# The return value is a dictionary with the energy of the decay chain
# and the energy of the decay chain in MeV
# The dictionary will look like this:
# {
#     {'Total Alpha Energy': 0.0,
#     'Total Beta Energy': 0.1332697613073299,
#      'Total Non-Penetrative Energy': 0.1479480421288555
# }

# Usecase 2: Complex decay chain
# We now try to get the energy of the decay chain of 224-Ra

E = energy_in_decay_chain(
    db=db,
    nuc_name="Ra-224",
)
# This will return the energy of the decay chain of Ra-224, and will return the energy of the decay chain of Ra-224.
# The fractional contribution of the energy of the decay chain of Ra-224 is handled
# by the function.

# An assumption of a fraction of 1.0 for the mother is assumed, this can be changed
# if desired.

# Use case 3 - See the daughters of a decay chain

# We can user the function get_daughters to see the daughters of a decay chain

daughters = get_daughters(
    db=db,
    nuclide_name="Ra-224",
)

print(daughters)

# Use case 4 - Inspection of the different radionuclide data

# The db is a dictionary with the data of the radionuclides
# The individual radionuclides are stored in the dictionary with the name of the radionuclide as the key
# and the object is the Radionuclide object, which has the following attributes:
# - name: The name of the radionuclide
# - halflife: The half life of the radionuclide with numerical value and unit
# - daughters: The daughters of the radionuclide, which is a list of Radionuclide objects
# - emissions: A list of the emissions of the radionuclide, which is a list of RadiationEmission objects
# - These RadiationEmission objects have the following attributes:
# - energy: The energy of the emission with numerical value in MeV
# - radiation_type: The type of the emission, which is a string following the convention of the ICRP
# - yield: The yield of the emission, which is a float between 0 and 1 per decay

radionuclide = db.get_decay_info("Ra-224")
print(radionuclide.name)
print(radionuclide.halflife)
print(radionuclide.daughters)
print(radionuclide.emissions)

# Note that only one sted of the decay chain is stored in each radionuclide, and 
# the daughters down the chain must be accessed the chain-functions


