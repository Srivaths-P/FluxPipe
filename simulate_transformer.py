import json
import os
import sys

import gdsfactory as gf
import matplotlib.pyplot as plt
import numpy as np
from openems import openems

import transformer_pcell

# Unit Constants
MM = 1e-3
UM = 1e-6


def setup_simulation(params: dict, sim_name: str = "xf_sim") -> openems.OpenEMS:
    """Configures and initializes an OpenEMS simulation for a transformer.

    Args:
        params: Dictionary of transformer geometry and simulation settings.
        sim_name: Base name for simulation files.

    Returns:
        An openems.OpenEMS object configured for the simulation.
    """
    f_start_ghz = params["f1_ghz"]
    f_stop_ghz = params["f2_ghz"]
    N1 = params["N1"]
    N2 = params["N2"]
    w_um = params["w_um"]
    s_um = params["s_um"]
    r1_pri_in_um = params["r1_primary_inner_um"]
    eps_r_sub = params["eps_r_substrate"]
    h_sub_um = params["h_substrate_um"]

    metal_th_um = params.get("metal_thickness_um", 1.0)
    sub_tand = params.get("substrate_tand", 0.002)
    box_pad_factor = params.get("sim_box_padding_factor", 1.2)
    mesh_res_factor = params.get("mesh_resolution_factor", 4)
    end_crit = params.get("EndCriteria", 1e-4)
    nrts = params.get("NrTS", 3000000)
    fsteps = params.get("fsteps_output", 101)
    add_common_ground_via = params.get("add_common_ground_via", True)

    # Generate GDS device using the parametric cell
    gf.clear_cache()
    xf_device = transformer_pcell.transformer(
        N1=N1, N2=N2, w=w_um, s=s_um, r1_pri_in=r1_pri_in_um, eps_r=eps_r_sub, h_sub=h_sub_um
    )
    print("[INFO] Plotting GDS layout for visualization...")
    xf_device.plot(return_fig=True).show()

    # Get device dimensions to size the simulation box
    r_pri_out_edge = xf_device.info.get("r_pri_out_edge_um", r1_pri_in_um + w_um / 2)
    r_sec_in_cl = xf_device.info.get("r1_sec_in_cl_um", r_pri_out_edge + s_um + w_um / 2)
    r_sec_out_edge = xf_device.info.get("r_sec_out_edge_um", r_sec_in_cl + w_um / 2)
    max_dev_r_um = max(r_pri_out_edge, r_sec_out_edge)
    if N1 == 0 and N2 == 0:
        max_dev_r_um = r1_pri_in_um + w_um

    # Initialize OpenEMS Simulation
    em = openems.OpenEMS(
        sim_name,
        EndCriteria=end_crit,
        fmin=f_start_ghz * 1e9,
        fmax=f_stop_ghz * 1e9,
        fsteps=fsteps,
        NrTS=nrts,
    )

    # Define Materials
    copper = openems.Metal(em, "copper")
    ground_plane_metal = openems.Metal(em, "ground_plane_metal")
    sub_mat = openems.Dielectric(
        em, "substrate", eps_r=eps_r_sub, tand=sub_tand, fc=f_stop_ghz * 1e9
    )

    # Define Geometry Primitives from parameters
    h_sub_m = h_sub_um * UM
    metal_th_m = metal_th_um * UM
    sub_xy_span = 2 * max_dev_r_um * box_pad_factor * UM
    sub_start_xy = -sub_xy_span / 2
    sub_stop_xy = sub_xy_span / 2

    openems.Box(
        sub_mat, 1, start=[sub_start_xy, sub_start_xy, 0], stop=[sub_stop_xy, sub_stop_xy, h_sub_m]
    )
    gnd_plane_thickness_m = metal_th_m
    openems.Box(
        ground_plane_metal,
        priority=10,
        start=[sub_start_xy, sub_start_xy, 0 - gnd_plane_thickness_m],
        stop=[sub_stop_xy, sub_stop_xy, 0],
    )
    print(f"[INFO] Added ground plane from z={-gnd_plane_thickness_m:.3e} to z=0.")

    # Create transformer traces from GDS polygons
    gds_polys_by_spec = xf_device.get_polygons(by="tuple")
    metal_z0_on_sub = h_sub_m
    metal_z1_on_sub = h_sub_m + metal_th_m
    dbu_to_um = 0.001  # GDS database units to micrometers

    copper_gds_layers = [(1, 0), (2, 0), transformer_pcell.COMMON_ARM_LAYER_TUPLE]
    for spec, polys_list in gds_polys_by_spec.items():
        if spec not in copper_gds_layers or not polys_list:
            continue
        for polygon in polys_list:
            poly_um_np = None
            if hasattr(polygon, "points") and isinstance(polygon.points, np.ndarray):
                poly_um_np = polygon.points * dbu_to_um
            elif hasattr(polygon, "each_point_hull"):
                pts_list = [[p.x * dbu_to_um, p.y * dbu_to_um] for p in polygon.each_point_hull()]
                if pts_list:
                    poly_um_np = np.array(pts_list)

            if poly_um_np is None or poly_um_np.size < 6:
                print(f"[WARN] Invalid or empty polygon for layer {spec}. Skipping.")
                continue

            openems.Polygon(
                copper,
                points=poly_um_np * UM,
                elevation=[metal_z0_on_sub, metal_z1_on_sub],
                normal_direction="z",
                priority=10,
            )

    # Add a cylindrical via to ground the common port
    if add_common_ground_via and (N1 > 0 or N2 > 0):
        via_radius_um = params.get("common_via_radius_um", w_um)
        via_radius_m = via_radius_um * UM
        openems.Cylinder(
            copper,
            priority=11,
            start=[0, 0, metal_z0_on_sub],
            stop=[0, 0, 0],
            radius=via_radius_m,
        )
        print(f"[INFO] Added common ground via at (0,0) with radius {via_radius_m:.3e} m.")
        em.mesh.AddLine("x", [-via_radius_m, via_radius_m])
        em.mesh.AddLine("y", [-via_radius_m, via_radius_m])

    # Define Boundary Conditions and Mesh
    em.boundaries = ["PEC"] * 6
    min_feat_um = min(w_um, s_um) if w_um > 0 and s_um > 0 else max(w_um, s_um, 1.0)
    base_res_m = (min_feat_um / mesh_res_factor) * UM if mesh_res_factor > 0 else min_feat_um * UM

    # Define Ports from GDS port objects
    port_mesh_x, port_mesh_y = [], []
    port_feed_len_m = 3 * base_res_m if base_res_m > 0 else 3 * UM
    gds_ports_map = {}
    if N1 > 0 and "o1" in xf_device.ports:
        gds_ports_map["P1"] = xf_device.ports["o1"]
    if N2 > 0 and "o2" in xf_device.ports:
        gds_ports_map["P2"] = xf_device.ports["o2"]

    processed_ports_info = []
    for port_label, gds_port in gds_ports_map.items():
        center_m = np.array(gds_port.center) * UM
        width_m = gds_port.width * UM
        angle_deg = gds_port.orientation
        dx, dy = np.cos(np.deg2rad(angle_deg)), np.sin(np.deg2rad(angle_deg))
        px, py = center_m[0], center_m[1]
        p0, p1, d_char = [0, 0, 0], [0, 0, 0], "x"

        if abs(dx) > abs(dy):  # X-directed port
            d_char = "x"
            y_cs = [py - width_m / 2, py + width_m / 2]
            port_mesh_y.extend(y_cs)
            x_cs = [px - port_feed_len_m, px] if dx > 0 else [px, px + port_feed_len_m]
            p0, p1 = [x_cs[0], y_cs[0], metal_z0_on_sub], [x_cs[1], y_cs[1], metal_z1_on_sub]
            port_mesh_x.extend(x_cs)
        else:  # Y-directed port
            d_char = "y"
            x_cs = [px - width_m / 2, px + width_m / 2]
            port_mesh_x.extend(x_cs)
            y_cs = [py - port_feed_len_m, py] if dy > 0 else [py, py + port_feed_len_m]
            p0, p1 = [x_cs[0], y_cs[0], metal_z0_on_sub], [x_cs[1], y_cs[1], metal_z1_on_sub]
            port_mesh_y.extend(y_cs)
        processed_ports_info.append({"start": p0, "stop": p1, "direction": d_char})

    mesh_xy_ext = max_dev_r_um * UM
    em.mesh.AddLine("x", sorted(list(set([-mesh_xy_ext, 0, mesh_xy_ext] + port_mesh_x))))
    em.mesh.AddLine("y", sorted(list(set([-mesh_xy_ext, 0, mesh_xy_ext] + port_mesh_y))))

    air_pad_z_m = sub_xy_span / 4
    z_lines = sorted(
        list(
            {
                -gnd_plane_thickness_m,
                0,
                h_sub_m,
                metal_z1_on_sub,
                metal_z1_on_sub + air_pad_z_m,
                -air_pad_z_m,
            }
        )
    )
    em.mesh.AddLine("z", z_lines)
    em.resolution = [base_res_m] * 3

    for p_info in processed_ports_info:
        openems.Port(em, start=p_info["start"], stop=p_info["stop"], direction=p_info["direction"], z=50)

    return em


def run_and_extract(em: openems.OpenEMS, num_threads: int = 8) -> tuple:
    """Runs the simulation and extracts S-parameters.

    Args:
        em: The configured OpenEMS simulation object.
        num_threads: The number of threads to use for the simulation.

    Returns:
        A tuple containing (frequencies_array, s_parameters_dictionary).
    """
    s_raw_vec = em.run_openems(options="solve", numThreads=num_threads, show_plot=False)

    if s_raw_vec is None:
        print("[ERROR] Simulation did not return S-parameters. Exiting.")
        sys.exit(1)

    freqs = getattr(em, "frequencies", None)
    if freqs is None:
        print("[WARN] Frequency vector not found in em.frequencies. Reconstructing.")
        freqs = np.linspace(em.fmin, em.fmax, em.fsteps)
        if not s_raw_vec or not em.ports or len(s_raw_vec[0]) != len(freqs):
            print("[ERROR] Frequency length mismatch after reconstruction. Exiting.")
            sys.exit(1)

    s = {}
    n_ports = len(em.ports)

    s["S11"] = s_raw_vec[0] if n_ports >= 1 else np.zeros_like(freqs, dtype=complex)
    s["S21"] = s_raw_vec[1] if n_ports >= 2 else np.zeros_like(freqs, dtype=complex)
    s["S12"] = s["S21"]  # Assume reciprocal passive network
    s["S22"] = np.zeros_like(freqs, dtype=complex)  # Placeholder

    if n_ports < 2:
        s["S21"] = np.zeros_like(freqs, dtype=complex)
        s["S12"] = np.zeros_like(freqs, dtype=complex)

    for key in ["S11", "S12", "S21", "S22"]:
        if key not in s or s[key].size != freqs.size:
            s[key] = np.zeros_like(freqs, dtype=complex)

    return freqs, s


def save_results_csv(freqs: np.ndarray, s_params: dict, out_dir: str, name_base: str):
    """Saves S-parameters (magnitude and phase) to a CSV file.

    Args:
        freqs: Array of frequency points.
        s_params: Dictionary of S-parameter data.
        out_dir: The output directory.
        name_base: The base name for the output file.
    """
    if freqs is None or freqs.size == 0 or not s_params:
        print("[WARN] No valid data provided to save to CSV.")
        return

    filepath = os.path.join(out_dir, f"{name_base}_sparameters.csv")
    headers = ["Frequency (Hz)"]
    columns = [freqs]

    for key in ["S11", "S21", "S12", "S22"]:
        s_data = s_params.get(key)
        if s_data is not None and s_data.size == freqs.size:
            columns.extend([np.abs(s_data), np.angle(s_data)])
            headers.extend([f"|{key}|", f"Phase({key}) (rad)"])
        else:
            columns.extend([np.zeros_like(freqs)] * 2)
            headers.extend([f"|{key}| (no data)", f"Phase({key}) (rad) (no data)"])

    data_matrix = np.array(columns).T
    try:
        np.savetxt(filepath, data_matrix, delimiter=",", header=",".join(headers), comments="")
        print(f"[INFO] S-parameters saved to: {filepath}")
    except Exception as e:
        print(f"[ERROR] Could not save CSV file to {filepath}: {e}")


def plot_results(freqs: np.ndarray, s_params: dict, out_dir: str, name_base: str):
    """Generates and saves plots of S-parameter magnitudes.

    Args:
        freqs: Array of frequency points.
        s_params: Dictionary of S-parameter data.
        out_dir: The output directory.
        name_base: The base name for the output file.
    """
    if freqs is None or freqs.size == 0 or not s_params:
        print("[WARN] No valid data provided for plotting.")
        return

    f_ghz = freqs / 1e9
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(f"S-Parameters: {name_base}", fontsize=15)

    s11, s22 = s_params.get("S11"), s_params.get("S22")
    if s11 is not None and s11.size == freqs.size:
        axes[0].plot(f_ghz, 20 * np.log10(np.abs(s11) + 1e-10), label="|S11| (dB)")
    if s22 is not None and np.any(np.abs(s22) > 1e-9):
        axes[0].plot(f_ghz, 20 * np.log10(np.abs(s22) + 1e-10), label="|S22| (dB)", linestyle="--")

    axes[0].set_title("Return Loss", fontsize=13)
    axes[0].set_xlabel("Frequency (GHz)")
    axes[0].set_ylabel("Magnitude (dB)")
    axes[0].set_ylim(-80, 5)
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend()

    s21, s12 = s_params.get("S21"), s_params.get("S12")
    if s21 is not None and s21.size == freqs.size:
        axes[1].plot(f_ghz, 20 * np.log10(np.abs(s21) + 1e-10), label="|S21| (dB)")
    if s12 is not None and np.any(np.abs(s12) > 1e-9):
        axes[1].plot(f_ghz, 20 * np.log10(np.abs(s12) + 1e-10), label="|S12| (dB)", linestyle="--")

    axes[1].set_title("Transmission / Coupling", fontsize=13)
    axes[1].set_xlabel("Frequency (GHz)")
    axes[1].set_ylabel("Magnitude (dB)")
    axes[1].set_ylim(-150, 5)
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend()

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])

    for ext in ["png", "pdf"]:
        filepath = os.path.join(out_dir, f"{name_base}_sparameters.{ext}")
        try:
            plt.savefig(filepath)
            print(f"[INFO] Plot saved to: {filepath}")
        except Exception as e:
            print(f"[ERROR] Could not save plot to {filepath}: {e}")
    plt.close(fig)


if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    configs_dir = os.path.join(script_dir, "configs")
    output_dir = os.path.join(script_dir, "outputs")

    config_files_to_run = []
    if len(sys.argv) > 1:
        config_path_arg = sys.argv[1]
        if not os.path.exists(config_path_arg):
            sys.exit(f"[ERROR] Configuration file not found: {config_path_arg}")
        config_files_to_run.append(config_path_arg)
    else:
        print(f"[INFO] No specific config file provided. Searching in: {configs_dir}")
        if not os.path.isdir(configs_dir):
            sys.exit(f"[ERROR] 'configs' directory not found at: {configs_dir}")
        json_files = [os.path.join(configs_dir, f) for f in os.listdir(configs_dir) if f.endswith(".json")]
        if not json_files:
            sys.exit(f"[ERROR] No JSON configuration files found in {configs_dir}")
        config_files_to_run.extend(json_files)
        print(f"[INFO] Found {len(config_files_to_run)} config(s) to process.")

    for config_path in config_files_to_run:
        print(f"\n{'---' * 10}\n[INFO] Processing: {config_path}\n{'---' * 10}")
        try:
            with open(config_path, "r") as f:
                sim_params = json.load(f)
        except Exception as e:
            print(f"[ERROR] Failed to load or parse JSON from {config_path}: {e}")
            continue

        sim_base_name = os.path.splitext(os.path.basename(config_path))[0]
        specific_output_dir = os.path.join(output_dir, sim_base_name)

        for dir_to_create in [output_dir, specific_output_dir]:
            if not os.path.exists(dir_to_create):
                try:
                    os.makedirs(dir_to_create)
                except OSError as e:
                    print(f"[ERROR] Failed to create directory {dir_to_create}: {e}")
                    continue
        print(f"[INFO] Output for this run will be in: {specific_output_dir}")

        try:
            print("[INFO] Generating reference GDS file...")
            gf.clear_cache()
            ref_device = transformer_pcell.transformer(
                N1=sim_params["N1"],
                N2=sim_params["N2"],
                w=sim_params["w_um"],
                s=sim_params["s_um"],
                r1_pri_in=sim_params["r1_primary_inner_um"],
                eps_r=sim_params["eps_r_substrate"],
                h_sub=sim_params["h_substrate_um"],
            )
            gds_filepath = os.path.join(specific_output_dir, f"{sim_base_name}.gds")
            ref_device.write_gds(gds_filepath)
            print(f"[INFO] Reference GDS file saved: {gds_filepath}")
        except Exception as e:
            print(f"[ERROR] GDS generation failed for {sim_base_name}: {e}")
            continue

        print("[INFO] Setting up OpenEMS simulation...")
        ems_sim_name = f"sim_{sim_base_name}"
        ems_environment = None
        try:
            ems_environment = setup_simulation(sim_params, sim_name=ems_sim_name)
        except Exception as e:
            print(f"[ERROR] OpenEMS setup failed for {sim_base_name}: {e}")
            continue

        if ems_environment is None:
            print(f"[WARN] Skipping run for {sim_base_name} due to setup failure.")
            continue

        num_threads = sim_params.get("numThreads", 8)
        print(f"[INFO] Running simulation with {num_threads} threads...")
        freq_data, sparam_data = None, None
        try:
            freq_data, sparam_data = run_and_extract(ems_environment, num_threads=num_threads)
        except Exception as e:
            print(f"[ERROR] OpenEMS simulation run failed for {sim_base_name}: {e}")

        if freq_data is not None and sparam_data is not None and freq_data.size > 0:
            print("[INFO] Saving results and plotting...")
            save_results_csv(freq_data, sparam_data, specific_output_dir, sim_base_name)
            plot_results(freq_data, sparam_data, specific_output_dir, sim_base_name)
        else:
            print(f"[WARN] Skipping results processing for {sim_base_name} due to lack of data.")

    print(f"\n{'---' * 10}\n[INFO] All configurations processed. Script finished.\n{'---' * 10}")

