import gdsfactory as gf
import numpy as np

# A dedicated GDS layer for the common connection arms of the spirals.
COMMON_ARM_LAYER_TUPLE = (100, 0)


@gf.cell
def transformer(
    N1: int = 10,
    N2: int = 6,
    w: float = 5.0,
    s: float = 5.0,
    r1_pri_in: float = 50.0,
    eps_r: float = 11.68,
    h_sub: float = 500.0,
    inner_connect_width_factor: float = 1.0,
) -> gf.Component:
    """Generates a gdsfactory component for a planar RF spiral transformer.

    The transformer consists of two concentric spiral inductors. The inner
    connection arms for both spirals meet at the center (0,0) and are placed
    on a separate GDS layer to prevent geometric overlaps.

    Args:
        N1: Number of primary turns.
        N2: Number of secondary turns.
        w: Trace width (µm).
        s: Edge-to-edge spacing between turns (µm).
        r1_pri_in: Inner radius of the primary coil's first turn centerline (µm).
        eps_r: Substrate relative permittivity (for metadata).
        h_sub: Substrate thickness (µm) (for metadata).
        inner_connect_width_factor: Multiplier for the common port width at (0,0).

    Returns:
        A gf.Component containing the transformer geometry.
    """
    c = gf.Component(name=f"transformer_N1_{N1}_N2_{N2}_w_{w}_s_{s}")
    pitch = w + s  # Centerline-to-centerline distance between turns

    def create_spiral_and_arm(
        parent_component: gf.Component,
        r_in_cl: float,
        turns: int,
        coil_layer_tuple: tuple,
        arm_layer_tuple: tuple,
        trace_width: float,
        pts_per_turn: int = 100,
    ) -> tuple:
        """Helper to create one spiral and its inner connection arm."""
        r_out_cl_val = r_in_cl
        last_point = np.array([r_in_cl, 0.0])
        second_last_point = np.array([r_in_cl, 0.0])

        if turns > 0:
            # Create the Archimedean spiral segment
            th_max = 2 * np.pi * turns
            n_pts_curve = max(3, int(turns * pts_per_turn))
            theta_spiral = np.linspace(0, th_max, n_pts_curve)
            r_values = r_in_cl + (pitch / (2 * np.pi)) * theta_spiral
            curved_points = np.column_stack([r_values * np.cos(theta_spiral), r_values * np.sin(theta_spiral)])

            if curved_points.shape[0] > 0:
                spiral_path = gf.Path(curved_points)
                parent_component.add_ref(spiral_path.extrude(width=trace_width, layer=coil_layer_tuple))
                last_point = curved_points[-1]
                second_last_point = curved_points[-2] if n_pts_curve > 1 else curved_points[0]
                r_out_cl_val = r_in_cl + (pitch / (2 * np.pi)) * th_max

            # Create the straight inner arm connecting to the center
            if r_in_cl > 1e-9:
                arm_points = np.array([[0.0, 0.0], [r_in_cl, 0.0]])
                arm_path = gf.Path(arm_points)
                parent_component.add_ref(arm_path.extrude(width=trace_width, layer=arm_layer_tuple))

        return r_out_cl_val, last_point, second_last_point

    # Generate Primary Coil (N1)
    primary_coil_layer = (1, 0)
    r_pri_out_cl = r1_pri_in
    if N1 > 0:
        r_pri_out_cl, p1_last, p1_seclast = create_spiral_and_arm(
            c, r1_pri_in, N1, primary_coil_layer, COMMON_ARM_LAYER_TUPLE, w
        )
        if not np.allclose(p1_last, p1_seclast):
            angle_rad = np.arctan2(p1_last[1] - p1_seclast[1], p1_last[0] - p1_seclast[0]) + np.pi / 2
            c.add_port(
                name="o1",
                center=p1_last.tolist(),
                width=w,
                orientation=np.rad2deg(angle_rad),
                layer=primary_coil_layer,
            )

    # Generate Secondary Coil (N2)
    secondary_coil_layer = (2, 0)
    if N1 > 0:
        r1_sec_in_cl = (r_pri_out_cl + w / 2) + s + (w / 2)
    else:
        r1_sec_in_cl = r1_pri_in

    r_sec_out_cl = r1_sec_in_cl
    if N2 > 0:
        r_sec_out_cl, p2_last, p2_seclast = create_spiral_and_arm(
            c, r1_sec_in_cl, N2, secondary_coil_layer, COMMON_ARM_LAYER_TUPLE, w
        )
        if not np.allclose(p2_last, p2_seclast):
            angle_rad = np.arctan2(p2_last[1] - p2_seclast[1], p2_last[0] - p2_seclast[0]) + np.pi / 2
            c.add_port(
                name="o2",
                center=p2_last.tolist(),
                width=w,
                orientation=np.rad2deg(angle_rad),
                layer=secondary_coil_layer,
            )

    # Add the common port at the center
    c.add_port(
        name="common",
        center=(0, 0),
        width=w * inner_connect_width_factor,
        orientation=0,
        layer=COMMON_ARM_LAYER_TUPLE,
    )

    # Store key parameters and dimensions in the component's metadata
    c.info.update(
        {
            "N1": N1,
            "N2": N2,
            "w_um": w,
            "s_um": s,
            "r1_pri_in_um": r1_pri_in,
            "eps_r": eps_r,
            "h_sub_um": h_sub,
            "r_pri_out_edge_um": (r_pri_out_cl + w / 2) if N1 > 0 else (r1_pri_in + w / 2),
            "r1_sec_in_cl_um": r1_sec_in_cl,
            "r_sec_out_edge_um": (r_sec_out_cl + w / 2) if N2 > 0 else (r1_sec_in_cl + w / 2),
        }
    )
    return c


if __name__ == "__main__":
    # This block demonstrates the p-cell functionality when run directly.
    print("Generating and plotting example transformer designs...")
    gf.clear_cache()

    # Example 1: A 1-to-3 turn transformer
    transformer_1 = transformer(N1=1, N2=3, w=20, s=20, r1_pri_in=30)
    transformer_1.plot()

    # Example 2: A 2-to-5 turn transformer with different dimensions
    transformer_2 = transformer(N1=2, N2=5, w=5, s=10, r1_pri_in=15)
    transformer_2.plot()

    print("Plots displayed. To save, use component.write_gds('filename.gds')")
