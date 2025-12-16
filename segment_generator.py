"""
Script to generate ilustrations of which regions are included for 
calculation in the red marrow dosimetry GUI module
"""

import numpy as np
import matplotlib.pyplot as plt

import nibabel as nib
import vtk

main_file = "resources/segment_ilus/segment.nii.gz"

main_segments = nib.load(main_file)

print(main_segments.header)

# List the unique values in the image
data = main_segments.get_fdata()
unique_values = np.unique(data)

print("Unique values in the image:", unique_values)

lumbar_IDs = [27, 28, 29, 30, 31]
thoracic_IDs = [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43]
cervical_IDs = [44, 45, 46, 47, 48, 49, 50]
scapulae_IDs = [71, 72]
clavicles_IDs = [73, 74]
humerus_IDs = [69, 70]
femur_IDs = [75, 76]
hip_bone_IDs = [77, 78]
skull_IDs = [91]
rib_IDs = [92, 93, 94, 95, 96, 97, 98, 99, 100,
               101, 102, 103, 104, 105, 106, 107,
               108, 109, 110, 111, 112, 113, 114, 115,
               ]
sacrum_IDs = [25]
sternum_IDs = [116]

single_IDs_lookup = {
    "Lumbar Vertebrae": lumbar_IDs,
    "Thoracic Vertebrae": thoracic_IDs,
    "Cervical Vertebrae": cervical_IDs,
    "Scapulae": scapulae_IDs,
    "Clavicles": clavicles_IDs,
    "Proximal Humeri": humerus_IDs,
    "Proximal Femora": femur_IDs,
    "Os Coxae": hip_bone_IDs,
    "Skull": skull_IDs,
    "Ribs": rib_IDs,
    "Sacrum": sacrum_IDs,
    "Sternum": sternum_IDs
}

# Combine IDs for each organ into a single value

main_IDs = {
    "Lumbar Vertebrae": int(1),
    "Thoracic Vertebrae": int(2),
    "Cervical Vertebrae": int(3),
    "Scapulae": int(4),
    "Clavicles": int(5),
    "Proximal Humeri": int(6),
    "Proximal Femora": int(7),
    "Os Coxae": int(8),
    "Skull": int(9),
    "Ribs": int(10),
    "Sacrum": int(11),
    "Sternum": int(12)

}

# Make a copy of the segmented data image

modified_data = np.zeros_like(data, dtype=np.int16)

for organ, ids in single_IDs_lookup.items():
    main_id = main_IDs[organ]
    for id_val in ids:
        modified_data[data == id_val] = main_id
    print(f"Assigned main ID {main_id} to organ {organ} for segment IDs {ids}")

# Save the modified image
modified_img = nib.Nifti1Image(modified_data, affine=main_segments.affine, header=main_segments.header)

modified_file = "resources/segment_ilus/segment_modified.nii.gz"
nib.save(modified_img, modified_file)
print(f"Modified image saved to {modified_file}")




