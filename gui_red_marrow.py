import tkinter as tk
from tkinter import ttk
import json
import nibabel as nib
import numpy as np
import vtk
from vtk.util import numpy_support  
from PIL import Image, ImageTk
from ashen.ashen_utils import (
    convert_to_hours
)

from ashen.myelodose_backend import (
    calculate_absorbed_dose_from_input_data
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

DEFAULT_VALUES_13 = {
    "cranifacial_bones": 38,
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

DEFAULT_VALUES_7 = {
    "cervical_vertebrae": None,
    "femur_head": None,
    "femur_neck": None,
    "iliac_crest": None,
    "lumbar_vertebrae": None,
    "ribs": None,
    "parietal_bone": None,
}


CHECKBOX_TITLES = json_data.get("CHECKBOX_TITLES", {})
# Mapping from checkbox title to value to display
CHECKBOX_VALUES = json_data.get("CHECKBOX_VALUES", {})

DO_SAVE_JSON = True

# --- Checkbox window class ------------------------------------------------------


class CheckBoxWindow(tk.Toplevel):
    """A reusable dynamic checkbox window."""

    def __init__(self, parent, selected_item):
        super().__init__(parent)
        self.parent = parent
        self.title("Additional Options")
        self.geometry("300x300")

        titles = CHECKBOX_TITLES.get(selected_item, [])

        tk.Label(self, text=f"Options for '{selected_item}':",
                 font=("Arial", 12, "bold")).pack(pady=10)

        self.vars = []
        frame = tk.Frame(self)
        frame.pack(pady=10)

        def make_callback(idx):
            def callback(*args):
                if self.vars[idx][1].get():
                    # Check all above
                    for j in range(idx):
                        self.vars[j][1].set(1)
                else:
                    # Uncheck all below
                    for j in range(idx+1, len(self.vars)):
                        self.vars[j][1].set(0)
            return callback

        for idx, title in enumerate(titles):
            var = tk.IntVar(value=1)  # All checked by default
            row_frame = tk.Frame(frame)
            row_frame.pack(anchor="w", fill="x")
            chk = tk.Checkbutton(row_frame, text=title, variable=var)
            chk.pack(side="left")
            # Show mapped value or N/A
            value = CHECKBOX_VALUES.get(title, "N/A")

            half_life_in_hours = convert_to_hours(float(value[:-1]), value[-1])  # Assuming last char is unit

            if half_life_in_hours > 24*10:

                value_label = tk.Label(row_frame, text=str(value), fg="red")

            else:
                value_label = tk.Label(row_frame, text=str(value), fg="gray")
            value_label.pack(side="left", padx=10)
            self.vars.append((title, var))
            var.trace_add('write', make_callback(idx))

        if not titles:
            tk.Label(self, text="No relevant daughters").pack(pady=10)

        tk.Button(self, text="Close", command=self.destroy).pack(pady=20)

    def get_values(self):
        """Return dict of checkbox selections."""
        return {title: var.get() for title, var in self.vars}


# --- Main window ----------------------------------------------------------------


class DynamicFormApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Dynamic Form")
        self.geometry("800x800")
        self.current_daughters = []  # Store daughters for selected nuclide
        self.checkbox_window = None  # Ensure checkbox_window is always defined


        # Extra numerical input field
        num_frame = tk.Frame(self)
        num_frame.pack(fill="x", padx=5, pady=5)
        tk.Label(num_frame, text="Alpha RBE value:").pack(side="left")
        self.alpha_RBE_var = tk.StringVar()
        self.alpha_RBE_entry = tk.Entry(num_frame, textvariable=self.alpha_RBE_var)
        self.alpha_RBE_entry.pack(side="left", padx=5)

        # Dropdown
        ttk.Label(self, text="Select or input radionuclide:").pack(anchor="w")
        self.combo = ttk.Combobox(self, values=PREDEFINED_VALUES)
        self.combo.bind("<<ComboboxSelected>>", self.on_nuclide_selected)
        self.combo.pack(fill="x", padx=5, pady=5)


        # Radio buttons
        self.radio_choice = tk.StringVar(value="RM")
        rad_frame = tk.Frame(self)
        rad_frame.pack(pady=5)
        tk.Radiobutton(rad_frame, text="Red marrow",
                       variable=self.radio_choice, value="RM").pack(side="left", padx=10)
        tk.Radiobutton(rad_frame, text="Trabecular bone surface",
                       variable=self.radio_choice, value="TBS").pack(side="left", padx=10)


        # Button → open second window
        tk.Button(self, text="Select daughters",
              command=self.open_checkbox_window).pack(pady=10)

        # Main frame for layout
        self.main_frame = tk.Frame(self)
        self.main_frame.pack(fill="both", expand=True)

        # Form container (left)
        self.form_container = tk.Frame(self.main_frame)
        self.form_container.grid(row=0, column=0, sticky="nsew")

        # VTK offscreen rendering to PNG, then display in Tkinter
        #import vtkmodules.all as vtk
        self.vtk_frame = tk.Frame(self.main_frame)
        self.vtk_frame.grid(row=0, column=1, sticky="nsew")

        #reader = vtk.vtkNIFTIImageReader()
        #reader.SetFileName("resources/segment_ilus/segment_modified.nii.gz")
        #reader.Update()

        #mc = vtk.vtkMarchingCubes()
        #mc.SetInputConnection(reader.GetOutputPort())
        #mc.SetValue(0, 0.5)
        #mc.Update()

        #mapper = vtk.vtkPolyDataMapper()
        #mapper.SetInputConnection(mc.GetOutputPort())
        #actor = vtk.vtkActor()
        #actor.SetMapper(mapper)

        #renderer = vtk.vtkRenderer()
        #renderer.AddActor(actor)
        #renderer.SetBackground(1, 1, 1)

        #render_window = vtk.vtkRenderWindow()
        #render_window.SetOffScreenRendering(1)
        #render_window.AddRenderer(renderer)
        #render_window.SetSize(400, 400)
        #render_window.Render()

        ## Capture the image to a PNG file
        #window_to_image = vtk.vtkWindowToImageFilter()
        #window_to_image.SetInput(render_window)
        #window_to_image.Update()

        #writer = vtk.vtkPNGWriter()
        png_path = "resources/vtk_rendered_surface.png"
        #writer.SetFileName(png_path)
        #writer.SetInputConnection(window_to_image.GetOutputPort())
        #writer.Write()

        # Display the PNG in Tkinter
        img = Image.open(png_path)

        # Get dimensions of the image
        dim = img.size

        # Resize and keep ratio

        img = img.resize((300, int(300 * dim[1] / dim[0])))

        self.current_image = ImageTk.PhotoImage(img)
        self.image_label = tk.Label(self.vtk_frame, image=self.current_image)
        self.image_label.pack(fill="both", expand=True)


        # Create both forms
        self.form_13, self.widgets_13 = self.create_form(13, ELECTRON_SITES)
        self.form_7, self.widgets_7 = self.create_form(7, ALPHA_SITES)


        self.active_widgets = self.widgets_13
        #self.update_image()

          # Buttons (create only once)
        btn_frame = tk.Frame(self)
        btn_frame.pack(pady=10)
        tk.Button(btn_frame, text="Apply ICRP CFs To All",
            command=self.apply_defaults_to_all).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Save To File",
            command=self.save_to_file).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Run Calculation",
            command=self.run_calculation).pack(side="left", padx=5)
    def save_to_file(self):
        """Save current input data to a file."""
        inputs = self.collect_all_values()
        with open("saved_input.json", "w") as f:
            json.dump(inputs, f, indent=4)
        print("Input data saved to saved_input.json")

        # Handle second window instance
        self.checkbox_window = None
    
    def on_nuclide_selected(self, event=None):
        selected = self.combo.get()
        self.current_daughters = CHECKBOX_TITLES.get(selected, [])
        print(f"Daughters for {selected}: {self.current_daughters}")
        self.update_form(event)

    # --- form creation -----------------------------------------------------

    def create_form(self, n, field_names=None):
        frame = tk.Frame(self.form_container)
        widgets = []
        dropdown_vals = [str(v) for v in range(10, 110, 10)]

        for i in range(n):
            row = i + 1
            if field_names and i < len(field_names):
                label_text = field_names[i]
            else:
                label_text = f"Field {row}:"

            tk.Label(frame, text=label_text).grid(row=i, column=0, sticky="w")

            entry = tk.Entry(frame)
            entry.grid(row=i, column=1, sticky="ew", padx=5)
            #entry.bind("<KeyRelease>", lambda event: self.update_image())

            combo = ttk.Combobox(frame, values=dropdown_vals, width=6)
            combo.grid(row=i, column=2, sticky="ew", padx=5)

            btn = tk.Button(frame, text="Use ICRP CF",
                            command=lambda c=combo, n=label_text: self.apply_default_combo(c, n))
            btn.grid(row=i, column=3, padx=5)

            widgets.append({'name': label_text, 'entry': entry, 'combo': combo, 'btn': btn})

        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(2, weight=1)
        return frame, widgets

    # --- dynamic form switching -------------------------------------------

    def show_form(self, form):
        for w in self.form_container.winfo_children():
            w.pack_forget()
        form.pack(fill="both", expand=True)

    def update_form(self, event=None):
        val = self.combo.get()
        if val in SEVEN_FIELD_VALUES:
            self.active_widgets = self.widgets_7
            self.show_form(self.form_7)
        else:
            self.active_widgets = self.widgets_13
            self.show_form(self.form_13)

    # --- defaults ----------------------------------------------------------

    def get_default_dict(self):
        return DEFAULT_VALUES_7 if self.active_widgets is self.widgets_7 else DEFAULT_VALUES_13

    def apply_default(self, entry, row):
        # Deprecated: No longer used for entry fields
        pass

    def apply_default_combo(self, combo, name):
        defaults = self.get_default_dict()
        value = defaults.get(name, "")
        combo.set(value)

    def apply_defaults_to_all(self):
        defaults = self.get_default_dict()
        for w in self.active_widgets:
            value = defaults.get(w['name'], "")
            w['combo'].set(value)

    # --- checkbox window ---------------------------------------------------

    def open_checkbox_window(self):
        if self.checkbox_window is not None and tk.Toplevel.winfo_exists(self.checkbox_window):
            return  # Already open

        selected = self.combo.get()
        self.checkbox_window = CheckBoxWindow(self, selected)

    # --- data collection ---------------------------------------------------

    def collect_all_values(self):
        """Returns ALL user inputs as a dictionary."""
        data = {}

        # Extra numerical value
        try:
            data["alpha_RBE_value"] = float(self.alpha_RBE_var.get())
        except ValueError:
            data["alpha_RBE_value"] = None

        # Selected item
        data["radionuclide"] = self.combo.get()

        # Radio choice
        data["source_tissue"] = self.radio_choice.get()

        # Dynamic form entries
        form_data = []
        for w in self.active_widgets:
            form_data.append({
                "name": w['name'],
                "MBqhrs_per_ml": w['entry'].get(),
                "CF": w['combo'].get(),
            })
        data["fields"] = form_data


        # Use the stored daughters list
        all_daughters_for_selected = self.current_daughters

        # Checkbox window values
        if self.checkbox_window is not None:
            data["daughters"] = self.checkbox_window.get_values()

        else:

            if all_daughters_for_selected:
                # If there are relevant daughters but window not opened, assume none selected
                data["daughters"] = {title: 1 for title in all_daughters_for_selected}

            else:
                data["daughters"] = {}

        return data

    # --- calculation -------------------------------------------------------

    def postprocess_results(self, results):
        """Post-process and display results."""
        print("\n=== Calculation Results ===")
        for k, v in results.items():
            print(k, ":", v)

        # First check if selected radionuclide and daughters are not alpha emitters

        nuclide = results["radionuclide"]
        daughters = results["daughters"]
        full_chain = [nuclide] + list(daughters.keys())

        if nuclide not in SEVEN_FIELD_VALUES and not any(d not in SEVEN_FIELD_VALUES for d in daughters):
            print("\nNote: Selected radionuclide is not an alpha emitter. No alpha-calculation needed.")
            results["only_electron_calculation"] = True
            return results
        
        else:
            results["only_electron_calculation"] = False

            # Now check if sites are compatible with alpha calculation

            if "iliac_crest" in [f['name'] for f in results['fields'] if f['MBqhrs_per_ml']]:
                print("Warning: 'iliac_crest' site is extra problematic")
                results['incompatible_sites'] = ['iliac_crest']
                results['sites_compatible'] = False
                return results

            sites = [f['name'] for f in results['fields'] if f['MBqhrs_per_ml']]

            incompatible_sites = [s for s in sites if s not in SITES_BOTH_ALPHA_ELECTRON]

            if len(incompatible_sites) > 0:
                print("\nWarning: The following selected sites are incompatible with alpha and electron calculations ")
                for s in incompatible_sites:
                    print(f" - {s}")
                results['incompatible_sites'] = incompatible_sites
                results['sites_compatible'] = False
                return results

            else:
                results['sites_compatible'] = True
                return results

    def run_calculation(self):
        inputs = self.collect_all_values()

        print("\n=== Calculation Input ===")
        for k, v in inputs.items():
            print(k, ":", v)

        print("\nStarting to process input...\n")

        # Remove the fields that have no value
        print("=== Processed Input ===")
        processed_fields = [f for f in inputs["fields"] if f["MBqhrs_per_ml"]]
        inputs["fields"] = processed_fields

        # Check that a radionuclide was selected
        if not inputs["radionuclide"]:
            print("Error: No radionuclide selected.")
            return

        inputs = self.postprocess_results(inputs)

        calc_result = calculate_absorbed_dose_from_input_data(inputs)

        if DO_SAVE_JSON:

        # Save all input to a JSON file for further processing and debugging

            with open("calculation_input.json", "w") as f:
                json.dump(inputs, f, indent=4)

            calc_result.save_to_json("calculation_results_test.json")

            print("\nCalculation input and results saved to JSON files.")



# --- run program ---------------------------------------------------------

if __name__ == "__main__":
    app = DynamicFormApp()
    app.mainloop()

