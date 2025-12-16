"""

Test script for the program - this should only 
be on developer versions

"""

import sys
import numpy as np

from ashen.ashen_utils import (
    make_decay_chain_db,
    load_icrp_107,
    energy_in_decay_chain,

)

from ashen.red_marrow_module  import (
    return_list_of_sites_by_nuclide,
)

from typing import Dict

# Test the functions from ashen_utils

# Load decay chains and emission data
emission_energy: Dict[str, float] = load_icrp_107()
db = make_decay_chain_db(emission_data=emission_energy)

f_18 = db.get_decay_info("F-18")
lu_177 = db.get_decay_info("Lu-177")
ac_225 = db.get_decay_info("Ac-225")
cf_252 = db.get_decay_info("Cf-252")

print(return_list_of_sites_by_nuclide(ac_225))

print(f_18.is_alpha, lu_177.is_alpha, ac_225.is_alpha, cf_252.is_alpha)
print(f_18.is_electron_emitter, lu_177.is_electron_emitter, ac_225.is_electron_emitter, cf_252.is_electron_emitter)
print(f_18.is_neutron_emitter, lu_177.is_neutron_emitter, ac_225.is_neutron_emitter, cf_252.is_neutron_emitter)
print(f_18.is_fission_emitter, lu_177.is_fission_emitter, ac_225.is_fission_emitter, cf_252.is_fission_emitter)

#for em in f_18.emissions:
#    print(em)
#
#for em in lu_177.emissions:
#    print(em)

list_of_emission_types = [em.radiation_type for em in cf_252.emissions]

if "A" in list_of_emission_types:
    print("Cf-252 is an alpha emitter")

print(np.unique(list_of_emission_types))


sys.exit()

E = energy_in_decay_chain(
    db=db,
    nuc_name="Ac-225",
)

E_2 = energy_in_decay_chain(
    db=db,
    nuc_name="Ac-225",
    rbe_alpha=5,
)

print(E)
print(E_2)
