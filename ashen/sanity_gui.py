"""

GUI for doing simple calculations with 
absorbed dose based on some simple assumptions
and a simple model.

Only self dose is considered, no other sources of dose are included.
Also, all energy is locally absorbed.

"""

import tkinter as tk
from tkinter import ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
from scipy.integrate import quad

from ashen.ashen_utils import (
    make_decay_chain_db,
    load_icrp_107,
    energy_in_decay_chain,
    mono_exp,
    f_2p,
)

DEBUG_MODE = False


def load_db():
    """
    Load the decay chain database if not already loaded.
    """
    global db
    if db is None:
        # Load decay chains and emission data
        emission_energy = load_icrp_107()
        db = make_decay_chain_db(emission_data=emission_energy)


def calculate_absorbed_dose_per_type(func, params, t_syn, db, mass, nuclide_name=None,
                                     rbe_alpha=1.0):
    """
    Calculate the absorbed dose for a given radionuclide and parameters.

    This is the function inside the sanity check GUI - in the 
    future this should be moved to the ashen package

    Parameters:
    - func: Function to calculate activity (e.g., mono_exp or f_2p).
    - params: Dictionary of parameters for the activity function.
    - t_syn: Time points for the activity calculation.
    - db: Decay chain database.
    - mass: Mass of the organ in kg.
    - nuclide_name: Name of the radionuclide (optional, for energy calculation).
    - rbe_alpha: Relative Biological Effectiveness for alpha radiation (default is 1.0).
    Returns:
    - result: Dictionary containing the total integrated activity, energies, and absorbed doses.
    """

    a = func(t_syn, **params)

    # Integrate the activity to get the absorbed dose

    I = quad(
        lambda t: func(t, **params), a=0, b=np.inf)

    TIA = I[0]  # Total Integrated Activity in MBq*h

    E_in_chain = energy_in_decay_chain(db, nuclide_name, rbe_alpha=rbe_alpha)

    alpha_E = E_in_chain["Total Alpha Energy"]  # in MeV
    beta_E = E_in_chain["Total Beta Energy"]  # in MeV
    total_E = E_in_chain["Total Non-Penetrative Energy"]  # in MeV

    # TODO: Not have the J per MeV as a magic number here
    alpha_E_J = alpha_E * 1.60218e-13  # Convert MeV to J
    beta_E_J = beta_E * 1.60218e-13  # Convert MeV to J
    total_E_J = total_E * 1.60218e-13  # Convert MeV to J

    total_destintegrations = TIA*1e6*3600  # Convert MBq*h to disintegrations

    total_E_alpha = total_destintegrations * alpha_E_J  # in J
    total_E_beta = total_destintegrations * beta_E_J  # in J
    total_E_non_penetrative = total_destintegrations * total_E_J  # in J

    total_AD_alpha = total_E_alpha / mass  # in Gy
    total_AD_beta = total_E_beta / mass  # in Gy
    total_AD_non_penetrative = total_E_non_penetrative / mass  # in Gy

    # Make a dictionary to return
    result = {
        "Total Integrated Activity (MBq*h)": TIA,
        "Total Alpha Energy (MeV)": alpha_E,
        "Total Beta Energy (MeV)": beta_E,
        "Total Non-Penetrative Energy (MeV)": total_E,
        "Absorbed Dose Alpha (Gy)": total_AD_alpha,
        "Absorbed Dose Beta (Gy)": total_AD_beta,
        "Absorbed Dose Non-Penetrative (Gy)": total_AD_non_penetrative,
        "RBE Alpha": rbe_alpha,
    }

    return result


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sanity check for absorbed dose calculations")
        self.geometry("1000x600")

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=3)

        emission_energy = load_icrp_107()
        self.db = make_decay_chain_db(emission_data=emission_energy)

        # Name of radionuclides

        self.nuclides = self.db.get_all_nuclide_names()

        self.create_widgets()

    def create_widgets(self):
        # Left panel
        self.input_frame = ttk.Frame(self, padding="10")
        self.input_frame.grid(row=0, column=0, sticky="nsew")

        left_panel_element_padding = 10

        # Right panel
        self.output_frame = ttk.Frame(self, padding="10")
        self.output_frame.grid(row=0, column=1, sticky="nsew")

        # Input fields
        ttk.Label(self.input_frame, text="A0 in mono-exp").grid(row=0,
                                                                column=0, sticky="w", pady=left_panel_element_padding)
        self.input_A0 = ttk.Entry(self.input_frame)
        self.input_A0.insert(0, "1")  # default
        self.input_A0.grid(row=0, column=1, pady=left_panel_element_padding)

        ttk.Label(self.input_frame, text="Biological half life").grid(
            row=1, column=0, sticky="w", pady=left_panel_element_padding)
        self.input_t_half_bio = ttk.Entry(self.input_frame)
        self.input_t_half_bio.insert(0, "30")  # default
        self.input_t_half_bio.grid(
            row=1, column=1, pady=left_panel_element_padding)

        ttk.Label(self.input_frame, text="Select radionuclide").grid(
            row=2, column=0, sticky="w", pady=left_panel_element_padding)
        # self.nuclide_var = tk.StringVar(value=self.nuclides[0])  # default to first nuclide
        self.nuclide_var = tk.StringVar(
            value="Lu-177")  # default to first nuclide
        self.nuclide_dropdown = ttk.Combobox(
            self.input_frame, textvariable=self.nuclide_var, values=self.nuclides, state="normal"
        )
        self.nuclide_dropdown.grid(
            row=2, column=1, sticky="ew", pady=left_panel_element_padding)

        # Button to insert administered activity
        ttk.Label(self.input_frame, text="Administered activity (MBq)").grid(
            row=3, column=0, sticky="w", pady=left_panel_element_padding)
        self.input_administered_activity = ttk.Entry(self.input_frame)
        self.input_administered_activity.insert(0, "1000")  # default
        self.input_administered_activity.grid(
            row=3, column=1, pady=left_panel_element_padding)

        # Button to insert fraction of administered activity in peak
        ttk.Label(self.input_frame, text="Fraction of administered activity in peak").grid(
            row=4, column=0, sticky="w", pady=left_panel_element_padding)
        self.input_fraction_peak = ttk.Entry(self.input_frame)
        self.input_fraction_peak.insert(0, "0.5")
        self.input_fraction_peak.grid(
            row=4, column=1, pady=left_panel_element_padding)

        # Button to insert mass of the object
        ttk.Label(self.input_frame, text="Mass of the organ (Kg)").grid(
            row=5, column=0, sticky="w", pady=left_panel_element_padding)
        self.input_mass = ttk.Entry(self.input_frame)
        self.input_mass.insert(0, "0.001")  # default
        self.input_mass.grid(row=5, column=1, pady=left_panel_element_padding)

        # Add radio buttons for half-life type
        self.halflife_type_var = tk.StringVar(value="bio")
        radio_frame = ttk.Frame(self.input_frame)
        radio_frame.grid(row=6, column=0, columnspan=2,
                         sticky="w", pady=left_panel_element_padding)
        ttk.Label(radio_frame, text="Half-life input type:").pack(side=tk.LEFT)
        ttk.Radiobutton(radio_frame, text="Biological",
                        variable=self.halflife_type_var, value="bio").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(radio_frame, text="Effective", variable=self.halflife_type_var,
                        value="effective").pack(side=tk.LEFT, padx=5)

        # Add button for effective half-life
        ttk.Label(self.input_frame, text="Effective half-life (h)").grid(row=7,
                                                                         column=0, sticky="w", pady=left_panel_element_padding)
        self.input_effective_half_life = ttk.Entry(self.input_frame)
        self.input_effective_half_life.insert(0, "10")
        self.input_effective_half_life.grid(
            row=7, column=1, pady=left_panel_element_padding)

        # Add checkbox for setting A0 manually or from administered activity and fraction
        ttk.Label(self.input_frame, text="Calculate A0 from administered activity and fraction").grid(
            row=8, column=0, columnspan=2, sticky="w", pady=left_panel_element_padding)
        self.calculate_A0_var = tk.BooleanVar(value=True)  # Default to True
        self.calculate_A0_checkbox = ttk.Checkbutton(
            self.input_frame,
            text="Calculate A0",
            variable=self.calculate_A0_var,
            command=lambda: self.input_A0.config(
                state="disabled" if self.calculate_A0_var.get() else "normal")
        )
        self.calculate_A0_checkbox.grid(
            row=9, column=0, columnspan=2, sticky="w", pady=left_panel_element_padding)

        # Input the RBE-value

        ttk.Label(self.input_frame, text="RBE for alpha radiation").grid(
            row=10, column=0, sticky="w", pady=left_panel_element_padding)
        self.input_rbe_alpha = ttk.Entry(self.input_frame)
        self.input_rbe_alpha.insert(0, "1.0")
        self.input_rbe_alpha.grid(
            row=10, column=1, pady=left_panel_element_padding)

        # Add bindings to trigger recalculate with <FocusOut> event

        self.input_A0.bind("<FocusOut>", lambda e: self.calculate())
        self.input_t_half_bio.bind("<FocusOut>", lambda e: self.calculate())
        self.input_administered_activity.bind(
            "<FocusOut>", lambda e: self.calculate())
        self.input_fraction_peak.bind("<FocusOut>", lambda e: self.calculate())
        self.input_mass.bind("<FocusOut>", lambda e: self.calculate())
        self.input_effective_half_life.bind(
            "<FocusOut>", lambda e: self.calculate())
        self.nuclide_dropdown.bind(
            "<<ComboboxSelected>>", lambda e: self.calculate())
        self.input_rbe_alpha.bind("<FocusOut>", lambda e: self.calculate())

        self.halflife_type_var.trace_add(
            "write", lambda *args: self.calculate())
        self.calculate_A0_var.trace_add(
            "write", lambda *args: self.calculate())

        ttk.Button(self.input_frame, text="Calculate", command=self.calculate).grid(
            row=12, column=0, columnspan=2, pady=10
        )

        # Matplotlib figure setup (2 plots stacked)
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.ax1 = self.figure.add_subplot(211)
        self.ax2 = self.figure.add_subplot(212)

        self.canvas = FigureCanvasTkAgg(self.figure, master=self.output_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Text output area
        self.text_output = tk.Text(self.output_frame, height=6)
        self.text_output.pack(fill="x", expand=False, pady=10)

        self.calculate()

    def calculate(self):

        try:
            t_half_bio = float(self.input_t_half_bio.get())
            administered_activity = float(
                self.input_administered_activity.get())
            fraction_peak = float(self.input_fraction_peak.get())
            mass = float(self.input_mass.get())
            effective_half_life = float(self.input_effective_half_life.get())
            rbe_alpha = float(self.input_rbe_alpha.get())

            if self.calculate_A0_var.get():
                # Calculate A0 from administered activity and fraction
                A0 = administered_activity * fraction_peak
                self.input_A0.delete(0, tk.END)
                self.input_A0.insert(0, str(A0))
            else:
                # Use the value from the input field
                A0 = float(self.input_A0.get())
            if self.halflife_type_var.get() == "effective":
                # Override and use effective half-life
                print(f"Using effective half-life: {effective_half_life}")
            elif self.halflife_type_var.get() == "bio":
                # Use biological half-life directly
                print(f"Using biological half-life: {t_half_bio}")
        except ValueError:
            self.text_output.delete("1.0", tk.END)
            self.text_output.insert(
                tk.END, "Invalid input! Please enter valid values.\n")
            return

        t_syn = np.linspace(0, 300, 1000)

        single_nuclide = self.db.get_decay_info(self.nuclide_var.get())
        lam_phys = single_nuclide.calculate_decay_constant()

        if self.halflife_type_var.get() == "effective":
            lam_eff = np.log(2) / effective_half_life
            a = f_2p(
                t_syn,
                A0=A0,
                lam_bio=None,
                lam_phys=lam_phys,
                use_effective_half_life=True,
                lam_eff=lam_eff

            )

            calc_lam_bio = lam_eff - lam_phys
            calc_t_half_bio = np.log(2) / calc_lam_bio

            if calc_lam_bio < 0:
                self.text_output.delete("1.0", tk.END)
                self.text_output.insert(
                    tk.END, "Effective half-life is shorter than physical half-life, biological half-life cannot be negative.\n")
                return

        else:
            lam_bio = np.log(2) / t_half_bio
            a = f_2p(t_syn, A0=A0, lam_bio=lam_bio, lam_phys=lam_phys)

            lam_eff = lam_bio + lam_phys
            calc_t_half_eff = np.log(2) / lam_eff

        a_only_phys = f_2p(
            t_syn, A0=A0, lam_phys=lam_phys, only_phys=True)

        # Set the X-axis limit to a sensible upper value
        # when a is one percent of A0
        # TODO: Get back to this...
        # cutoff_a = A0*0.01
        # a_sub = a[50:]
        # t_max = t_syn[np.where(a_sub > cutoff_a)[0]]
        # print(f"t_max: {t_max}")
        t_max = 10*24

        # Calculate the absorbed dose
        result = calculate_absorbed_dose_per_type(
            func=f_2p,
            params={
                "A0": A0,
                "lam_bio": lam_bio if self.halflife_type_var.get() == "bio" else None,
                "lam_phys": lam_phys,
                "use_effective_half_life": self.halflife_type_var.get() == "effective",
                "lam_eff": lam_eff if self.halflife_type_var.get() == "effective" else None
            },
            t_syn=t_syn,
            mass=mass,
            db=self.db,
            nuclide_name=self.nuclide_var.get(),
            rbe_alpha=rbe_alpha
        )

        plot_results = {
            "AD_alpha": result["Absorbed Dose Alpha (Gy)"],
            "AD_beta": result["Absorbed Dose Beta (Gy)"],
            "AD_total": result["Absorbed Dose Non-Penetrative (Gy)"],
        }

        # Update plots
        self.ax1.clear()
        self.ax1.plot(t_syn, a, label='Mono-exponential decay')
        self.ax1.plot(t_syn, a_only_phys,
                      label='Physical decay only', linestyle='--')
        self.ax1.set_title("Activity vs Time")
        self.ax1.set_xlabel("Time (h)")
        self.ax1.set_ylabel("Activity (MBq)")
        self.ax1.set_xlim(-1, t_max)
        self.ax1.legend()

        self.ax2.clear()
        self.ax2.bar(['Alpha', 'Beta', 'Total'],
                     [plot_results["AD_alpha"], plot_results["AD_beta"],
                         plot_results["AD_total"]],
                     color=['#66c2a5', '#fc8d62', '#8da0cb'])
        self.ax2.set_ylabel("Absorbed Dose (Gy)")

        self.canvas.draw()
        plt.tight_layout()

        # Update text output
        self.text_output.delete("1.0", tk.END)
        # self.text_output.insert(tk.END, f"Plotted sine and cosine with frequency={freq}, amplitude={amp}\n")
        self.text_output.insert(
            tk.END, f"The selected radionuclide is: {self.nuclide_var.get()}\n")

        if self.halflife_type_var.get() == "effective":
            self.text_output.insert(
                tk.END, f"Using effective half-life: {effective_half_life} h, biological half life is: {round(calc_t_half_bio,2)} h\n")
            self.text_output.insert(
                tk.END, f"Aborbed Dose Alpha: {round(plot_results['AD_alpha'], 3)} Gy\n")
            self.text_output.insert(
                tk.END, f"Aborbed Dose Beta: {round(plot_results['AD_beta'], 3)} Gy\n")
            self.text_output.insert(
                tk.END, f"Aborbed Dose Total: {round(plot_results['AD_total'], 3)} Gy\n")
        else:
            self.text_output.insert(
                tk.END, f"Using biological half-life: {t_half_bio} h, effective half life is: {round(calc_t_half_eff, 1)} h\n")
        
        self.text_output.insert(
            tk.END, f"Total Integrated Activity: {result['Total Integrated Activity (MBq*h)']:.2f} MBq*h\n")
        self.text_output.insert(
            tk.END, f"Aborbed Dose Alpha: {round(plot_results['AD_alpha'], 3)} Gy\n")
        self.text_output.insert(
            tk.END, f"Aborbed Dose Beta: {round(plot_results['AD_beta'], 3)} Gy\n")
        self.text_output.insert(
            tk.END, f"Aborbed Dose Total: {round(plot_results['AD_total'], 3)} Gy\n")


if __name__ == "__main__":

    if DEBUG_MODE:

        emission_energy = load_icrp_107()
        db = make_decay_chain_db(emission_data=emission_energy)
        nuclide_name = "Lu-177"

        print("Running in debug mode")

        A0 = 1
        lam_phys = np.log(2) / 120
        lam_bio = np.log(2) / 30
        lam_eff = lam_bio + lam_phys
        func_to_pass = f_2p
        mass = 0.001  # in Kg

        result_test = calculate_absorbed_dose_per_type(
            func=func_to_pass,
            params={
                "A0": A0,
                "lam_bio": lam_bio,
                "lam_phys": lam_phys,
                "use_effective_half_life": False,
                "lam_eff": lam_eff
            },
            t_syn=np.linspace(0, 100, 100),
            mass=mass,
            db=db,
            nuclide_name=nuclide_name
        )

        print("Test result:")
        for key, value in result_test.items():
            print(f"{key}: {value}")

    else:

        app = App()
        app.mainloop()
