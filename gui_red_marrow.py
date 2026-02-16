class ToolTip:
    """Create a tooltip for a given widget"""
    def __init__(self, widget, text_func_or_str):
        self.widget = widget
        self.text_func_or_str = text_func_or_str
        self.tipwindow = None
        widget.bind("<Enter>", self.show_tip)
        widget.bind("<Leave>", self.hide_tip)

    def get_text(self):
        if callable(self.text_func_or_str):
            return self.text_func_or_str()
        return self.text_func_or_str

    def show_tip(self, event=None):
        if self.tipwindow or not self.get_text():
            return
        x, y, _, cy = self.widget.bbox("insert") if hasattr(self.widget, "bbox") else (0,0,0,0)
        x = x + self.widget.winfo_rootx() + 25
        y = y + cy + self.widget.winfo_rooty() + 25
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(tw, text=self.get_text(), background="#ffffe0", relief="solid", borderwidth=1)
        label.pack()

    def hide_tip(self, event=None):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

import tkinter as tk
from tkinter import ttk

import ttkbootstrap as ttk
from ttkbootstrap.constants import *

import json
import nibabel as nib
import numpy as np
import vtk

from vtk.util import numpy_support  
from PIL import Image, ImageTk
from collections import defaultdict

import pandas as pd

from ashen.ashen_utils import (
    convert_to_hours
)

from ashen.myelodose_backend import (
    calculate_absorbed_dose_from_input_data,
    CombinedCalculationResults,
    post_process_back_end_inputs,
    make_plot_figure,
    check_for_warnings,
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
        "craniofacial_bones",
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
    "craniofacial_bones": 38,
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
    def show_warning_popup(self):
        """
        Show a custom Toplevel window with warnings and Proceed/Cancel buttons.
        Returns True if user chooses to proceed, False if cancelled.
        """
        if not hasattr(self, 'warnings') or not self.warnings:
            return True
        warning_win = tk.Toplevel(self)
        warning_win.title("Warning")
        warning_win.geometry("400x250")
        warning_win.transient(self)
        warning_win.grab_set()
        tk.Label(warning_win, text="Warnings detected:", font=("TkDefaultFont", 12, "bold"), fg="red").pack(pady=(15,5))
        msg_frame = tk.Frame(warning_win)
        msg_frame.pack(fill="both", expand=True, padx=10)
        msg_box = tk.Text(msg_frame, wrap="word", height=8, width=45)
        msg_box.insert("1.0", "\n".join(self.warnings))
        msg_box.config(state="disabled", bg=warning_win.cget("bg"))
        msg_box.pack(fill="both", expand=True)
        btn_frame = tk.Frame(warning_win)
        btn_frame.pack(pady=10)
        result = {'proceed': False}
        def proceed():
            result['proceed'] = True
            warning_win.destroy()
        def cancel():
            result['proceed'] = False
            warning_win.destroy()
        tk.Button(btn_frame, text="Proceed", command=proceed, width=12).pack(side="left", padx=12)
        tk.Button(btn_frame, text="Cancel", command=cancel, width=12).pack(side="left", padx=12)
        self.wait_window(warning_win)
        return result['proceed']

    def __init__(self):
        super().__init__()
        self.title("Myelodose Beta")
        self.geometry("900x900")
        self.current_daughters = []  # Store daughters for selected nuclide
        self.checkbox_window = None  # Ensure checkbox_window is always defined

        self.warnings = []  # Store warnings

        # Extra numerical input field
        num_frame = tk.Frame(self)
        num_frame.pack(fill="x", padx=5, pady=5)
        tk.Label(num_frame, text="Alpha RBE value:").pack(side="left")
        self.alpha_RBE_var = tk.StringVar(value="1.0")
        self.alpha_RBE_entry = tk.Entry(num_frame, textvariable=self.alpha_RBE_var)
        self.alpha_RBE_entry.pack(side="left", padx=5)

        # Example: Add a tooltip to the Alpha RBE entry
        def rbe_tooltip_text():
            val = self.alpha_RBE_var.get()
            try:
                v = float(val)
                if v > 1.5:
                    return f"Warning: High RBE value ({v})"
                elif v < 0.5:
                    return f"Warning: Low RBE value ({v})"
                else:
                    return f"Current RBE: {v} (normal range)"
            except Exception:
                return "Enter a numeric RBE value."
        ToolTip(self.alpha_RBE_entry, rbe_tooltip_text)

        # Dropdown
        ttk.Label(self, text="Select or input radionuclide:").pack(anchor="w")
        self.combo = ttk.Combobox(self, values=PREDEFINED_VALUES)
        self.combo.bind("<<ComboboxSelected>>", self.on_nuclide_selected)
        self.combo.pack(fill="x", padx=5, pady=5)

        # Example: Tooltip for radionuclide combobox, static text
        ToolTip(self.combo, "Select a radionuclide. Choices affect available fields.")

        # Radio buttons
        self.radio_choice = tk.StringVar(value="RM")
        rad_frame = tk.Frame(self)
        rad_frame.pack(pady=5)
        tk.Radiobutton(rad_frame, text="Red marrow",
                       variable=self.radio_choice, value="RM").pack(side="left", padx=10)
        tk.Radiobutton(rad_frame, text="Trabecular bone surface",
                       variable=self.radio_choice, value="TBS").pack(side="left", padx=10)



        # Suppress warnings radio button
        self.supress_warnings = tk.StringVar(value="no")
        warn_frame = tk.Frame(self)
        warn_frame.pack(pady=5)
        tk.Label(warn_frame, text="Suppress warnings:").pack(side="left")
        tk.Radiobutton(warn_frame, text="Yes", variable=self.supress_warnings, value="yes").pack(side="left")
        tk.Radiobutton(warn_frame, text="No", variable=self.supress_warnings, value="no").pack(side="left")

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

        # Input unit selection (radio buttons) will be placed inside the form area
        self.input_unit = tk.StringVar(value="MBqhrs_per_ml")  # Default unit

        # Create both forms (radio buttons will be inside)
        self.form_13, self.widgets_13, self.input_label_13 = self.create_form(13, ELECTRON_SITES)
        self.form_7, self.widgets_7, self.input_label_7 = self.create_form(7, ALPHA_SITES)

        self.active_widgets = self.widgets_13

        # Buttons (create only once)
        btn_frame = tk.Frame(self)
        btn_frame.pack(pady=10)
        tk.Button(btn_frame, text="Apply ICRP CFs To All",
            command=self.apply_defaults_to_all).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Save To File",
            command=self.save_to_file).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Run Calculation",
            command=self.run_calculation).pack(side="left", padx=5)
        tk.Button(btn_frame, text="Clear Input",
            command=self.clear_all_inputs).pack(side="left", padx=5)

        # Add Results button for demonstration
        tk.Button(btn_frame, text="Show Results Window",
            command=self.show_results_window).pack(side="left", padx=5)

    def clear_all_inputs(self):
        # Clear alpha RBE entry
        self.alpha_RBE_var.set("")
        # Reset radio buttons
        self.radio_choice.set("RM")
        # Clear all form entries and combos
        for widgets in [self.widgets_13, self.widgets_7]:
            for w in widgets:
                w['entry'].delete(0, 'end')
                w['combo'].set("")

    def show_results_window(self):
        # Show warning popup before proceeding
        #if hasattr(self, 'warnings') and self.warnings:
        #    if not self.show_warning_popup():
        #        print("Operation cancelled due to warnings.")
        #        return
        results = self.run_calculation()
        if results is None:
            return
        rows = results.prepare_rows_for_plotting()
        ResultsWindow(self, rows, results)

    def save_to_file(self):
        # Show warning popup before proceeding
        if hasattr(self, 'warnings') and self.warnings:
            if not self.show_warning_popup():
                print("Save cancelled due to warnings.")
                return
        inputs = self.collect_all_values()
        with open("saved_input.json", "w") as f:
            json.dump(inputs, f, indent=4)
        print("Input data saved to saved_input.json")
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

        # Input unit radio buttons above the labels
        unit_frame = tk.Frame(frame)
        unit_frame.grid(row=0, column=0, columnspan=4, sticky="ew", pady=(0, 2))
        tk.Label(unit_frame, text="Input unit:").pack(side="left", padx=(0, 5))
        tk.Radiobutton(unit_frame, text="Total Activity (total_act)", variable=self.input_unit, value="MBqhrs").pack(side="left")
        tk.Radiobutton(unit_frame, text="Concentration (conc_act)", variable=self.input_unit, value="MBqhrs_per_ml").pack(side="left")

        # Add headers above input and combo columns, update label based on unit
        input_label = tk.Label(frame, text=self.get_input_label_text(), font=("TkDefaultFont", 10, "bold"))
        input_label.grid(row=1, column=1, sticky="ew", padx=5, pady=(0, 2))
        tk.Label(frame, text="Cellularity factor", font=("TkDefaultFont", 10, "bold")).grid(row=1, column=2, sticky="ew", padx=5, pady=(0, 2))

        for i in range(n):
            row = i + 2
            if field_names and i < len(field_names):
                label_text = field_names[i]
            else:
                label_text = f"Field {row - 1}:"

            tk.Label(frame, text=label_text).grid(row=row, column=0, sticky="w")

            entry = tk.Entry(frame)
            entry.grid(row=row, column=1, sticky="ew", padx=5)

            combo = ttk.Combobox(frame, values=dropdown_vals, width=6)
            combo.grid(row=row, column=2, sticky="ew", padx=5)

            btn = tk.Button(frame, text="Use ICRP CF",
                            command=lambda c=combo, n=label_text: self.apply_default_combo(c, n))
            btn.grid(row=row, column=3, padx=5)

            widgets.append({'name': label_text, 'entry': entry, 'combo': combo, 'btn': btn})

        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(2, weight=1)
        return frame, widgets, input_label

    def get_input_label_text(self):
        unit = self.input_unit.get()
        if unit == "total_act":
            return "MBq*hrs (total activity)"
        else:
            return "MBq*hrs/ml (concentration)"

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
            # Hide ICRP CF buttons for 7-field nuclides
            for w in self.widgets_7:
                w['btn'].grid_remove()
        else:
            self.active_widgets = self.widgets_13

            # Update both input labels when unit changes
            def update_input_labels(*args):
                self.input_label_13.config(text=self.get_input_label_text())
                self.input_label_7.config(text=self.get_input_label_text())
            self.input_unit.trace_add('write', update_input_labels)
            self.show_form(self.form_13)
            # Show ICRP CF buttons for 13-field nuclides
            for w in self.widgets_13:
                w['btn'].grid()

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

        # Input unit
        data["input_unit"] = self.input_unit.get()

        # Dynamic form entries
        form_data = []
        for w in self.active_widgets:
            form_data.append({
                "name": w['name'],
                "value": w['entry'].get(),
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

    def run_calculation(self):

        # Check if there are reasons for warning
        
        # DEBUG

        inputs = self.collect_all_values()
        inputs = post_process_back_end_inputs(inputs)
        calc_result = calculate_absorbed_dose_from_input_data(inputs)

        # Function that takes in input and returns True if there are warnings to show, False if not

        if DO_SAVE_JSON:
            with open("debug_calc_result.json", "w") as f:
                json.dump(self.collect_all_values(), f, indent=4)

        self.warnings = check_for_warnings(calc_result)

        # Show warning popup before proceeding
        if hasattr(self, 'warnings') and self.warnings:
            if not self.show_warning_popup():
                print("Calculation cancelled due to warnings.")
                return None

        if DO_SAVE_JSON:
            with open("calculation_input.json", "w") as f:
                json.dump(inputs, f, indent=4)
            calc_result.save_to_json("new_output_results.json")
            print("\nCalculation input and results saved to JSON files.")
        return calc_result

# --- Window for displaying the results --------------------------------------------



class ResultsWindow(tk.Toplevel):
    def __init__(self, parent, rows, results):
        super().__init__(parent)
        self.title("Absorbed dose results")
        self.geometry("700x500")

        main = ttk.Frame(self)
        main.pack(fill="both", expand=True, padx=10, pady=10)

        main.columnconfigure(0, weight=1)
        main.rowconfigure(0, weight=1)

        self._build_table(main, rows)
        self._build_totals(rows)

        # Set a label explaining the asterisk
        if any(r.get("surrogate_electron_site_used", False) for r in rows) and any(r.get("electron_unity_saf_used", False) for r in rows):
            tk.Label(self, text="* Unity AF used for surrogate electron site calculation", font=("TkDefaultFont", 8, "italic")).pack()
        elif any(r.get("electron_unity_saf_used", False) for r in rows):
            tk.Label(self, text="* Surrogate electron site used for dose calculation", font=("TkDefaultFont", 8, "italic")).pack()
        else:
            tk.Label(self, text="", font=("TkDefaultFont", 8, "italic")).pack()
        

        self.plot_window = None
        plot_btn = tk.Button(self, text="Show Plot", command=lambda: self.show_plot(rows, results))
        plot_btn.pack(pady=10)


        # Make a button to export results to a report

        export_btn = tk.Button(self, text="Export to Report", command=self.export_to_report)
        export_btn.pack(pady=10)

    def show_plot(self, rows, results):
        # Only one PlotWindow at a time
        if self.plot_window is not None and self.plot_window.winfo_exists():
            self.plot_window.lift()
            return
        self.plot_window = PlotWindow(self, rows, results)

    def export_to_report(self):
        # Placeholder for export functionality
        print("Exporting results to report...")

        report_path = "absorbed_dose_report.txt"

        with open(report_path, "w") as f:
            f.write("Absorbed Dose Report\n")
            f.write("====================\n\n")
            f.write("Radionuclide\tSite\tDose Alpha (Gy)\tDose Electron (Gy)\n")
            for row in self.tree.get_children():
                values = self.tree.item(row)['values']
                f.write(f"{values[0]}\t{values[1]}\t{values[2]}\t{values[3]}\n")
        print(f"Report saved to {report_path}")

    def _copy_tree_selection(self, event=None):
        tree = self.tree
        selection = tree.selection()

        if not selection:
            return

        rows = []
        for item in selection:
            values = tree.item(item, "values")
            rows.append("\t".join(str(v) for v in values))

        text = "\n".join(rows)

        tree.clipboard_clear()
        tree.clipboard_append(text)

    def _build_table(self, parent, rows):
        columns = (
            "radionuclide",
            'site',
            "dose_alpha",
            "dose_electron",
        )

        tree = ttk.Treeview(parent, columns=columns, show="headings")
        tree.grid(row=0, column=0, sticky="nsew", padx=(0, 10))
    
        headings = {
            "radionuclide": "Radionuclide",
            "site": "Site",
            "dose_alpha": "Absorbed dose α (Gy)",
            "dose_electron": "Absorbed dose e⁻ (Gy)",
        }

        for col in columns:
            tree.heading(col, text=headings[col])
            tree.column(col, anchor="center")

        for r in rows:

            if r.get("surrogate_electron_site_used", False):
                r["site_display"] = f"{r['site']}*"
            else:
                r["site_display"] = r["site"]

            tree.insert(
                "",
                "end",
                values=(
                    r["radionuclide"],
                    r["site_display"],
                    f"{r['dose_alpha']:.4g}",
                    f"{r['dose_electron']:.4g}",
                ),
            )

        self.tree = tree

        tree.bind("<Control-c>", self._copy_tree_selection)

    def _fmt_sig(self, value, sig=3):
        if value == 0:
            return "0"
        return f"{value:.{sig}g}"

    def _build_totals(self, rows):
        """
        Docstring for _build_totals
        
        :param self: The instance of the class
        :param rows: The data rows containing absorbed dose 
        """

        totals = defaultdict(lambda: {
            "alpha": 0.0,
            "electron": 0.0,
            "surrogate": False})

        for r in rows:
            site = r["site"]
            totals[site]["alpha"] += r["dose_alpha"]
            totals[site]["electron"] += r["dose_electron"]

            if r.get("surrogate_electron_site_used", False):
                totals[site]["surrogate"] = True

        frame = ttk.Frame(self)
        frame.pack(fill="x", pady=(10, 0))
        
        # Headers

        headers = ("Site", "Total Absorbed dose α (Gy)", "Total Absorbed dose e⁻ (Gy)", "Total Absorbed dose (Gy)")

        for col, text in enumerate(headers):
            ttk.Label(
                frame,
                text=text,
                font=("TkDefaultFont", 10, "bold"),
            ).grid(row=0, column=col, sticky="w", padx=(0,15))

        ttk.Separator(frame, orient="horizontal").grid(
            row=1, column=0, columnspan=len(headers), sticky="ew", pady=(2,4))
        
        self._totals_copy_lines = [
            "\t".join(headers)
        ]

        for row_idx, (site, dose_data) in enumerate(totals.items(), start=2):
            alpha = dose_data["alpha"]
            electron = dose_data["electron"]
            total = alpha + electron
            surrogate_mark = "*" if dose_data["surrogate"] else ""

            ttk.Label(
                frame,
                text=f"{site}{surrogate_mark}",
            ).grid(row=row_idx, column=0, sticky="w", padx=(0,15))

            ttk.Label(
                frame,
                text=self._fmt_sig(alpha),
            ).grid(row=row_idx, column=1, sticky="w", padx=(0,15))

            ttk.Label(
                frame,
                text=self._fmt_sig(electron),
            ).grid(row=row_idx, column=2, sticky="w", padx=(0,15))

            ttk.Label(
                frame,
                text=self._fmt_sig(total),
            ).grid(row=row_idx, column=3, sticky="w", padx=(0,15))

            self._totals_copy_lines.append(
                f"{site}{surrogate_mark}\t{self._fmt_sig(alpha)}\t{self._fmt_sig(electron)}\t{self._fmt_sig(total)}"
            )


        #total_alpha = sum(r["dose_alpha"] for r in rows)
        #total_electron = sum(r["dose_electron"] for r in rows)

        #frame = ttk.Frame(self)
        #frame.pack(fill="x", pady=(10, 0))

        #ttk.Label(
        #    frame,
        #    text=f"Total absorbed dose  |  α: {total_alpha:.3g} Gy   e⁻: {total_electron:.3g} Gy",
        #    font=("TkDefaultFont", 10, "bold"),
        #).pack(anchor="center")

class PlotWindow(tk.Toplevel):
    def __init__(self, parent, rows, results):
        super().__init__(parent)
        self.title("Absorbed dose plot")
        self.geometry("600x500")
        self._build_plot(rows, results)

        # Add a radio button to toggle RBE-adjusted alpha doses
        # if there are alpha doses present
        if any(r["dose_alpha"] > 0 for r in rows):

            self.plot_rbe_adjusted_alpha = tk.BooleanVar(value=False)
            rbe_check = tk.Checkbutton(self, text="RBE-adjusted alpha doses",
                                       variable=self.plot_rbe_adjusted_alpha,
                                       command=lambda: self._update_plot(rows, results))
            rbe_check.pack(pady=5)

    def _build_plot(self, rows, results):
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.figure import Figure

        fig = make_plot_figure(results)

        #fig = Figure(figsize=(5, 4))
        #ax = fig.add_subplot(111)

        #nuclides = [r["radionuclide"] for r in rows]
        #dose_electron = [r["dose_electron"] for r in rows]
        #dose_alpha = [r["dose_alpha"] for r in rows]

        #ax.bar(nuclides, dose_electron, label="Electrons")
        #ax.bar(nuclides, dose_alpha, bottom=dose_electron, label="Alpha")

        #ax.set_ylabel("Absorbed dose (Gy)")
        #ax.set_title("Absorbed dose per radionuclide")
        #ax.legend()
        #ax.tick_params(axis="x", rotation=45)

        #fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master=self)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas = canvas

    def _update_plot(self, rows, results):
        fig = make_plot_figure(
            results,
            plot_rbe_adjusted_alpha=self.plot_rbe_adjusted_alpha.get()
        )
        self.canvas.figure = fig
        self.canvas.draw()

    # ---------------- TOTALS ----------------

    def _build_totals(self, rows):
        total_alpha = sum(r["dose_alpha"] for r in rows)
        total_electron = sum(r["dose_electron"] for r in rows)

        frame = ttk.Frame(self)
        frame.pack(fill="x", pady=(10, 0))

        ttk.Label(
            frame,
            text=f"Total absorbed dose  |  α: {total_alpha:.3g} Gy   e⁻: {total_electron:.3g} Gy",
            font=("TkDefaultFont", 10, "bold"),
        ).pack(anchor="center")

# --- run program ---------------------------------------------------------

if __name__ == "__main__":
    app = DynamicFormApp()
    app.mainloop()

