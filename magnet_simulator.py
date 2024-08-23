import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt

#configuration
n_magnets = 16 #number of magnet per ring
n_rings = 16 #number of rings
cube_length = 10 #length of cube in mm
magnet_strength = 5492 #magnet strength in gauss

ring_diameter = 100 #ring diameter in mm
ring_length = 200 #ring length in mm
plot_length= 40 #diameter length field map in mm
dsv_diameter = 50 #diameter of the spherical volume in mm


def field_compute(show_graph=False):
    angles = np.linspace(0, 360, n_magnets, endpoint=False)
    halbach = magpy.Collection()

    #protect code if n_rings = 1
    if n_rings != 1:
        ring_spacing = ring_length/(n_rings-1)
    else:
        ring_spacing = 1

    z_starting_pos = -ring_length/2
    for r in range(n_rings):
        #setting XYZ position per ring
        z_pos = z_starting_pos + r * ring_spacing
        x_pos = ring_diameter/2
        y_pos = 0
        for a in angles:
            #Calculate position to set ring diameter
            cube = magpy.magnet.Cuboid(
                dimension=(cube_length,cube_length,cube_length),
                magnetization=(magnet_strength/1.25663753*100,0,0), #converted using B = T/u
                position=(x_pos,y_pos,z_pos)
            )
            cube.rotate_from_angax(a, 'z', anchor=0)
            cube.rotate_from_angax(a, 'z')
            halbach.add(cube)

    #change 3-D graph view
    fig = halbach.show(backend='plotly', return_fig=True)
    #update axis labels
    fig.update_layout(
        scene=dict(
            xaxis_title="X-axis (mm)",
            yaxis_title="Y-axis (mm)",
            zaxis_title="Z-axis (mm)"
        )
    )
    # Show the figure with updated axis labels

    if show_graph:
        fig.show()
    return halbach

def calculate_ppm(halbach, dsv_diameter):
    # Define the grid points over the DSV
    r = dsv_diameter / 2
    grid = np.mgrid[-r:r:35j, -r:r:35j, -r:r:35j].T.reshape(-1, 3)
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

# Compute and plot field on x-y grid
def plot_xy(halbach, print_result=False, arrow=False):
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
    #calculate and display uniformity in ppm
    uniformity_ppm, B_average, B_max_deviation = calculate_ppm(halbach, dsv_diameter)
    plt.suptitle(f"Average field: {B_average*1000:.2f}mT   Max deviation: {B_max_deviation*1000:.2f}mT   Homogeneity: {uniformity_ppm:.0f} ppm over {dsv_diameter} mm DSV", fontsize=16)
    plt.tight_layout()

    # Printing for excel spread sheet comparison
    if print_result:
        result = [f'{B_average*1000:.2f}', f'{B_max_deviation*1000:.2f}', f'{uniformity_ppm:.0f}', dsv_diameter, ring_diameter, ring_length]
        formatted_data = '\t'.join(map(str, result))
        print (formatted_data)
    plt.show()


# Repeated simulator for excel
'''
starting_dsv_diameter = 50
ending_dsv_diameter = 50
for sdd in range(starting_dsv_diameter, ending_dsv_diameter+10, 10):
    halbach_collection = field_compute(show_graph=True)
    plot_xy(halbach_collection,print_result=True,arrow=True)
    ring_diameter += 20
    ring_length += 20
    plot_length += 10
    dsv_diameter += 10
'''

halbach_collection = field_compute(show_graph=True)
plot_xy(halbach_collection, print_result=True, arrow=True)

