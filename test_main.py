"""

Test script for the program - this should only 
be on developer versions

"""

from ashen_utils import (
    make_decay_chain_db,
    load_icrp_107,
    energy_in_decay_chain,
)

from typing import Dict

# Load decay chains and emission data
emission_energy: Dict[str, float] = load_icrp_107()
db = make_decay_chain_db(emission_data=emission_energy)

E = energy_in_decay_chain(
    db=db,
    nuc_name="Lu-177",
)

print(E)
