"""
GUI for the red marrow module

"""

import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from red_marrow_module import calculate_absorbed_dose_to_red_marrow, AbsorbedDoseResult

# --- Dummy site data ---
site_data = {"Lumbar": 234.144, "Femur": 312.0, "Pelvis": 450.2}

class DoseApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Single site estimator")
        self.geometry("1000x600")
        self.configure(bg="#f7f7f7")
        
        # Main layout: 2 columns
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=2)
        
        self.create_input_frame()
        self.create_output_frame()

    # --- LEFT PANEL: Inputs ---
    def create_input_frame(self):
        frame = ttk.Frame(self, padding=20)
        frame.grid(row=0, column=0, sticky="nsew")

        ttk.Label(frame, text="Input Parameters", font=("Helvetica", 14, "bold")).grid(row=0, column=0, columnspan=2, pady=(0, 10))

        # Site
        ttk.Label(frame, text="Site:").grid(row=1, column=0, sticky="w")
        self.site_var = tk.StringVar(value="Lumbar")
        ttk.Combobox(frame, textvariable=self.site_var, values=list(site_data.keys())).grid(row=1, column=1, sticky="ew")

        # Cellularity factor
        ttk.Label(frame, text="Cellularity factor (CF):").grid(row=2, column=0, sticky="w")
        self.cf_var = tk.StringVar(value="70")
        ttk.Combobox(frame, textvariable=self.cf_var, values=["10", "20", "30", "40", "50", "60", "70", "80", "90", "100"]).grid(row=2, column=1, sticky="ew")

        # Source
        ttk.Label(frame, text="Source region:").grid(row=3, column=0, sticky="w")
        self.source_var = tk.StringVar(value="RM")
        ttk.Combobox(frame, textvariable=self.source_var, values=["RM", "TBS"]).grid(row=3, column=1, sticky="ew")

        # Target
        ttk.Label(frame, text="Target region:").grid(row=4, column=0, sticky="w")
        self.target_var = tk.StringVar(value="RM")
        ttk.Combobox(frame, textvariable=self.target_var, values=["RM"]).grid(row=4, column=1, sticky="ew")

        # Nuclide
        ttk.Label(frame, text="Nuclide:").grid(row=5, column=0, sticky="w")
        self.nuclide_var = tk.StringVar(value="Ac-225")
        ttk.Entry(frame, textvariable=self.nuclide_var).grid(row=5, column=1, sticky="ew")

        # Volume
        ttk.Label(frame, text="Spongiosa volume (cm³):").grid(row=6, column=0, sticky="w")
        self.volume_var = tk.DoubleVar(value=site_data["Lumbar"])
        ttk.Entry(frame, textvariable=self.volume_var).grid(row=6, column=1, sticky="ew")

        # Activity
        ttk.Label(frame, text="Cumulative activity (MBq·s):").grid(row=7, column=0, sticky="w")
        self.activity_var = tk.DoubleVar(value=1.0)
        ttk.Entry(frame, textvariable=self.activity_var).grid(row=7, column=1, sticky="ew")

        # Add a checkbox for doing single nuclide if needed
        self.single_nuclide_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frame, text="Single nuclide calculation", variable=self.single_nuclide_var).grid(row=8, column=0, columnspan=2, sticky="w")


        # Calculate button
        ttk.Button(frame, text="Calculate", command=self.update_output).grid(row=9, column=0, columnspan=2, pady=20)

        for i in range(10):
            frame.rowconfigure(i, pad=5)
        frame.columnconfigure(1, weight=1)

    # --- RIGHT PANEL: Results ---
    def create_output_frame(self):
        frame = ttk.Frame(self, padding=20)
        frame.grid(row=0, column=1, sticky="nsew")

        ttk.Label(frame, text="Results", font=("Helvetica", 14, "bold")).pack(anchor="w")

        self.result_label = ttk.Label(frame, text="Absorbed Dose: — Gy", font=("Helvetica", 12))
        self.result_label.pack(anchor="w", pady=10)

        # Matplotlib figure placeholder
        self.fig, self.ax = plt.subplots(figsize=(5, 4))
        self.ax.set_title("Dose Visualization")
        self.ax.set_xlabel("Parameter")
        self.ax.set_ylabel("Value")
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    # --- Business logic placeholder ---
    def update_output(self):
        cf = float(self.cf_var.get())
        vol = float(self.volume_var.get())
        act = float(self.activity_var.get())

        source = self.source_var.get()
        target = self.target_var.get()

        nuclide = self.nuclide_var.get()

        dose_result = calculate_absorbed_dose_to_red_marrow(
            input_site = self.site_var.get(),
            input_cf = cf,
            input_spongiosa_volume = vol,
            input_cumulative_activity = act,
            source_name = source,
            target_name = target,
            nuclide_name = nuclide
        )

        single_dose = dose_result.get_total_ad()

        # Replace this with your real dose calculation

        self.result_label.config(text=f"Absorbed Dose: {single_dose:.3f} Gy")

        # Example plot update (you can replace with your data)
        self.ax.clear()
        self.ax.set_title("Dose Visualization")
        self.ax.plot([0, 1, 2, 3], [single_dose * i for i in range(4)], marker="o")
        self.canvas.draw()

# --- Run app ---
if __name__ == "__main__":
    app = DoseApp()
    app.mainloop()
