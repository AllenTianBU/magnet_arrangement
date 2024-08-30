from multiprocessing.util import MAXFD

import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt
from plotly.graph_objs.layout.xaxis import Minor

#  Configuration
n_magnets = 16 #number of magnet per ring
n_rings = 8 #number of rings
cube_length = 10 #length of cube in mm
magnet_strength = 5492 #magnet strength in gauss

ring_diameter = 92 #ring diameter in mm
ring_length = 134 #ring length in mm
plot_length= 40 #diameter length field map in mm
dsv_diameter = 50 #diameter of the spherical volume in mm
eclipse_ratio = 0.983 #minor axis/major axis ratio

# Configuration for second layer of rain
ring_diameter_second_layer = 120
n_magnet_second_layer = 20
n_ring_second_layer = 2
ring_length_second_layer = 30

def field_compute(halbach):
    angles = np.linspace(0, 360, n_magnets, endpoint=False)

    # Protect code if n_rings = 1
    if n_rings != 1:
        ring_spacing = ring_length/(n_rings-1)
    else:
        ring_spacing = 1

    z_starting_pos = -ring_length/2
    for r in range(n_rings):
        #setting XYZ position per ring
        z_pos = z_starting_pos + r * ring_spacing
        y_pos = 0
        x_pos = ring_diameter / 2
        for ang in angles:
            # Calculating Eclipse
            a = ring_diameter / 2
            b = a*eclipse_ratio
            x_pos = (a*b)/(np.sqrt((b*np.cos(np.radians(ang)))**2 + (a*np.sin(np.radians(ang)))**2))

            # Calculate position to set ring diameter
            cube = magpy.magnet.Cuboid(
                dimension=(cube_length,cube_length,cube_length),
                magnetization=(magnet_strength/1.25663753*100,0,0), #converted using B = T/u
                position=(x_pos,y_pos,z_pos)
            )
            cube.rotate_from_angax(ang, 'z', anchor=0)
            cube.rotate_from_angax(ang, 'z')
            halbach.add(cube)
    return halbach

def add_second_layer_ring(halbach):
    angles = np.linspace(0, 360, n_magnet_second_layer, endpoint=False)
    ring_spacing = ring_length_second_layer / (n_ring_second_layer - 1)
    z_starting_pos = -ring_length_second_layer / 2
    for r in range(n_ring_second_layer):
        # setting XYZ position per ring
        z_pos = z_starting_pos + r * ring_spacing
        x_pos = ring_diameter_second_layer / 2
        y_pos = 0
        for a in angles:
            # Calculate position to set ring diameter
            cube = magpy.magnet.Cuboid(
                dimension=(cube_length, cube_length, cube_length),
                magnetization=(magnet_strength / 1.25663753 * 100, 0, 0),  # converted using B = T/u
                position=(x_pos, y_pos, z_pos)
            )
            cube.rotate_from_angax(a, 'z', anchor=0)
            cube.rotate_from_angax(a, 'z')
            halbach.add(cube)
    return halbach


def calculate_ppm(halbach, dsv_diameter):
    # Define the grid points over the DSV
    r = dsv_diameter / 2
    grid = np.mgrid[-r:r:20j, -r:r:20j, -r:r:20j].T.reshape(-1, 3)
    distances = np.linalg.norm(grid, axis=1)
    grid = grid[distances <= r]

    # Get magnetic field values at each grid point
    B = halbach.getB(grid)
    B_magnitude = np.linalg.norm(B, axis=1)

    # Calculate the average and maximum deviation
    B_avg = np.mean(B_magnitude)
    B_max_deviation = np.max(np.abs(B_magnitude - B_avg))

    # Calculate uniformity in ppm
    uniformity_ppm = (B_max_deviation / B_avg) * 1e6
    return uniformity_ppm, B_avg, B_max_deviation

def plot_halback(halbach):
    # Change 3-D graph view
    fig = halbach.show(backend='plotly', return_fig=True)

    # Update axis labels
    fig.update_layout(
        scene=dict(
            xaxis_title="X-axis (mm)",
            yaxis_title="Y-axis (mm)",
            zaxis_title="Z-axis (mm)"
        )
    )

    # Show the figure with updated axis labels
    fig.show()

# Compute and plot field on x-y grid
def plot_xy(halbach, print_result=False, arrow=False,plot_show =True):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    z_positions = [-plot_length/2, 0, plot_length/2]
    for z_pos, ax in zip(z_positions, axes):
        grid = np.mgrid[-plot_length/2:plot_length/2:100j, -plot_length/2:plot_length/2:100j, z_pos:z_pos:1j].T[0]
        X, Y, _ = np.moveaxis(grid, 2, 0)

        B = halbach.getB(grid)
        Bx, By, _ = np.moveaxis(B, 2, 0)
        Bamp = np.linalg.norm(B, axis=2)*1000
        pc = ax.contourf(X, Y, Bamp, levels=35, cmap="jet")
        if arrow == True:
            ax.streamplot(X, Y, Bx, By, color="k", density=1.5, linewidth=1)

        # Add colorbar
        fig.colorbar(pc, ax=ax, label="|B| mT")

        # Figure styling
        ax.set(
            title=f"Field at z = {z_pos} mm",
            xlabel="x-position",
            ylabel="y-position",
            aspect=1,
        )

    # Calculate and display uniformity in ppm
    uniformity_ppm, B_average, B_max_deviation = calculate_ppm(halbach, dsv_diameter)
    plt.suptitle(f"Average field: {B_average*1000:.2f}mT   Max deviation: {B_max_deviation*1000:.2f}mT   Homogeneity: {uniformity_ppm:.0f} ppm over {dsv_diameter} mm DSV", fontsize=16)
    plt.tight_layout()

    # Printing for excel spread sheet comparison
    if print_result:
        result = [f'{B_average*1000:.2f}', f'{uniformity_ppm:.0f}',f'{B_max_deviation*1000:.2f}', dsv_diameter, ring_diameter, ring_length, n_rings*n_magnets, n_magnets, n_rings, eclipse_ratio]
        formatted_data = '\t'.join(map(str, result))
        print (formatted_data)
    if plot_show == True:
        plt.show()

# Add array to keep the data, and mark the optimized one. Lastly set halbach to optimized one.
def optimize_eclipse(halbach, starting_eclipse_ratio,ending_eclipse_ratio,increments):
    while starting_eclipse_ratio <= ending_eclipse_ratio:
        eclipse_ratio = starting_eclipse_ratio
        temp_halbach = field_compute(halbach)
        plot_xy(temp_halbach,print_result=True,arrow=True,plot_show=False)
        starting_eclipse_ratio += increments

# For Optimization
'''
for x in range(983, 984, 1):
    eclipse_ratio = x/1000
    halbach = magpy.Collection()
    halbach = field_compute(halbach)
    plot_halback(halbach)
    plot_xy(halbach, print_result=True, arrow=True, plot_show=True)

'''
halbach = magpy.Collection()

halbach = field_compute(halbach)

#optimize_eclipse(halbach, starting_eclipse_ratio = 0.95,ending_eclipse_ratio = 1,increments = 0.01)

#halbach = add_second_layer_ring(halbach)

plot_halback(halbach)

plot_xy(halbach, print_result=True, arrow=True)
