"""

Test code for publication related to red marrow corrections.

"""

def plot_skeletal_site_data_electrons(sites):

    fig = plt.figure()

    dir_to_plots = "RM_manuscript"

    for site in sites:

        electron_safs = SKELETAL_SITE_DATA_ELECTRONS[site]["electrons"]["red_marrow"]["icrp"]

        icrp_mass = ircrp_masses[site]

        e = electron_safs.keys()
        saf_values = electron_safs.values()

        saf_array = np.array(list(saf_values)) * icrp_mass

        plt.plot(e, saf_array, marker='o', label=f"{site} - electrons")

    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel("Electron Energy (MeV)")

    plt.legend()

    plt.savefig(f"{dir_to_plots}/skeletal_site_safs_electrons.png", dpi=300)

def test_energy_inpolation(saf_data = SKELETAL_SITE_DATA_ELECTRONS):

    sites = list(saf_data.keys())

    dir_to_plots = "RM_manuscript"

    cols = sns.color_palette("husl", len(sites))

    E_to_test = np.linspace(0.001, 10, 1000)

    for interp_technique in ["spline", "linear", "loglog"]:

        print(f"Testing interpolation technique: {interp_technique}")

        interpolations = {}

        for site in sites:
            saf_values = saf_data[site]["electrons"]["red_marrow"]["icrp"]
            interpolated_values = []
            for E in E_to_test:
                phi = interpolate_phi(energy=E,
                                      saf_data=saf_values,
                                      interpolation_technique="spline")
                interpolated_values.append(phi)

            interpolated_values = np.array(interpolated_values)
            interpolated_values = interpolated_values * ircrp_masses[site]

            interpolations[site] = interpolated_values

        fig = plt.figure()
        fig.suptitle(f"Interpolation technique: {interp_technique}")
        for site in sites:
            plt.plot(E_to_test, interpolations[site], label=site, color=cols[sites.index(site)])
            saf_values = saf_data[site]["electrons"]["red_marrow"]["icrp"]
            e = np.array(list(saf_values.keys()))
            saf_vals = np.array(list(saf_values.values())) * ircrp_masses[site]
            plt.scatter(e, saf_vals, color=cols[sites.index(site)], marker='o', s=75)

        plt.xscale('log')

        plt.xlabel("Electron Energy (MeV)")
        plt.ylabel("Interpolated phi value")
        plt.legend()

        plt.savefig(f"{dir_to_plots}/interpolation_{interp_technique}.png", dpi=300)


    plt.show()




    print(sites)

    return 0