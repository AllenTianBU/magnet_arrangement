import numpy as np
import magpylib as magpy
import matplotlib.pyplot as plt
import magsimulator as magsim
import random


# Generate halbach array
def generate_magnet(halbach):
    # Protect code if n_rings = 1
    if n_rings != 1:
        ring_spacing = ring_length / (n_rings - 1)
    else:
        ring_spacing = 1

    # Define starting z position from the bottom, so center is at (0,0,0)
    z_starting_pos = -ring_length / 2
    z_starting_pos = 0
    # Loop to generate ring by ring
    for r in range(n_rings):
        # Setting diameter for each ring for optimization
        change = 0
        ring_diameter_array = [
            ring_diameter ,
            ring_diameter ,
            ring_diameter ,
            ring_diameter ,
            ring_diameter + 5.3+change,
            ring_diameter + 5.4+change,
            ring_diameter + 4.9+change,
            ring_diameter + 5.3+change,
            ring_diameter + 5.3+change,
            ring_diameter + 4.9+change,
            ring_diameter + 5.4+change,
            ring_diameter + 5.3+change,
            ring_diameter ,
            ring_diameter ,
            ring_diameter ,
            ring_diameter ,
        ]

        # Setting XYZ position per ring
        z_pos = z_starting_pos + r * ring_spacing
        y_pos = 0

        # Define the angles between each magnet in a circle
        angles = np.linspace(0, 360, n_magnets, endpoint=False)

        for ang in angles:
            # Calculate Eclipse
            a = ring_diameter_array[r] / 2
            b = a * eclipse_ratio
            x_pos = (a * b) / (np.sqrt((b * np.cos(np.radians(ang))) ** 2 + (a * np.sin(np.radians(ang))) ** 2))
            #random_strength = magnet_strength*random.uniform(0.99,1.01)
            random_strength = magnet_strength
            random_strength = random_strength / 1.25663753 * 100
            # Defining the position and rotation of magnet
            cube = magpy.magnet.Cuboid(
                dimension=(cube_length, cube_length, cube_length),
                magnetization=(random_strength,0,0),  # converted using B = T/u
                position=(x_pos, y_pos, z_pos)
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
        # Setting XYZ position per ring
        z_pos = z_starting_pos + r * ring_spacing
        x_pos = ring_diameter_second_layer / 2
        y_pos = 0
        for ang in angles:
            # Calculating Eclipse
            a = ring_diameter_second_layer / 2
            b = a * eclipse_ratio
            x_pos = (a * b) / (np.sqrt((b * np.cos(np.radians(ang))) ** 2 + (a * np.sin(np.radians(ang))) ** 2))

            # Calculate position to set ring diameter
            cube = magpy.magnet.Cuboid(
                dimension=(cube_length, cube_length, cube_length),
                magnetization=(magnet_strength / 1.25663753 * 100, 0, 0),  # converted using B = T/u
                position=(x_pos, y_pos, z_pos)
            )
            cube.rotate_from_angax(ang, 'z', anchor=0)
            cube.rotate_from_angax(ang, 'z')
            halbach.add(cube)
    return halbach


def calculate_ppm(halbach, dsv):
    # Define the grid points over the DSV
    r = dsv / 2
    grid = np.mgrid[-r:r:40j, -r:r:40j, -r:r:40j].T.reshape(-1, 3)

    # Spherical shape is defined
    distances = np.linalg.norm(grid, axis=1)
    grid = grid[distances <= r]

    # Get magnetic field values at each grid point
    B = halbach.getB(grid)
    B_magnitude = np.linalg.norm(B, axis=1)

    # Calculate the mean and maximum deviation
    B_avg = np.mean(B_magnitude)
    B_min = np.min(B_magnitude)
    B_max = np.max(B_magnitude)
    B_max_deviation = B_max - B_min

    # Calculate uniformity in ppm
    uniformity_ppm = (B_max_deviation / B_avg) * 1e6
    return uniformity_ppm, B_avg, B_max_deviation


def plot_halback(halbach):
    # Defining graph
    fig = halbach.show(backend='plotly', return_fig=True)
    # Update axis labels
    fig.update_layout(
        scene=dict(
            xaxis_title="X-axis (mm)",
            yaxis_title="Y-axis (mm)",
            zaxis_title="Z-axis (mm)"
        )
    )
    # Show the figure
    fig.show()


# Compute and plot field on x-y grid
def plot_result(halbach, print_result=False, arrow=False, plot_show=True):
    sensor = magsim.define_sensor_points_on_sphere(20, dsv/2, [0,0,0])
    homogeneity, meanfield = magsim.cost_func(sensor.getB(halbach))

    if plot_show == True:
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        # Set the z position view of the three graphs
        z_positions = [-plot_length / 2, 0, plot_length / 2]

        for z_pos, ax in zip(z_positions, axes):
            # Grid is defined
            grid = np.mgrid[-plot_length / 2:plot_length / 2:50j, -plot_length / 2:plot_length / 2:50j, z_pos:z_pos:1j].T[
                0]

            # Separate grid in to X and Y axis
            X, Y, _ = np.moveaxis(grid, 2, 0)
            B = halbach.getB(grid)

            # Magnitude of field vector is taken and converted to mT
            Bamp = np.linalg.norm(B, axis=2) * 1000

            # Define field map
            pc = ax.contourf(X, Y, Bamp, levels=35, cmap="jet")

            # Show field vector
            if arrow == True:
                Bx, By, _ = np.moveaxis(B, 2, 0)
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

        # Calculate and display uniformity in ppm, field strength, and max deviation on graph
        #uniformity_ppm, B_average, B_max_deviation = calculate_ppm(halbach, dsv)
        plt.suptitle(
            f"Average field: {meanfield*1000:.2f}mT Homogeneity: {homogeneity:.0f} ppm over {dsv} mm DSV",
            fontsize=16)
        plt.tight_layout()
        plt.show()

    # Print for excel spread sheet comparison. Have no effect on outcome.
    if print_result:
        result = [f'{meanfield * 1000}', f'{homogeneity:.0f}', dsv,
                  ring_diameter, ring_length, n_rings * n_magnets, n_magnets, n_rings, eclipse_ratio]
        formatted_data = '\t'.join(map(str, result))
        print(formatted_data)


# Set the eclipse ratio automatically per ring
def optimize_eclipse(starting_eclipse_ratio,ending_eclipse_ratio,increments):
    optimized_eclipse_array = []
    while starting_eclipse_ratio <= ending_eclipse_ratio:
        eclipse_ratio = starting_eclipse_ratio
        temp_halbach = magpy.Collection()
        temp_halbach = generate_magnet(temp_halbach)
        homogeneity_result = plot_result(temp_halbach,print_result=True,arrow=False,plot_show=False)
        optimized_eclipse_array.append([homogeneity_result, eclipse_ratio])
        starting_eclipse_ratio += increments
    return optimized_eclipse_array


if __name__ == '__main__':
    # Configuration
    n_magnets = 18  # Number of magnet per ring
    n_rings = 1  # Number of rings
    cube_length = 10  # Length of cube in mm
    #magnet_strength = 5754  # Magnet strength in gauss
    magnet_strength = 12467 # Measured

    ring_diameter = 93  # Ring diameter in mm
    ring_length = 180  # Ring length in mm
    plot_length = 45  # Length of field map plotted in mm
    dsv = 50  # Diameter of the spherical volume in mm
    eclipse_ratio = 1  # Minor axis/major axis ratio

    # Configuration for second layer of rain
    ring_diameter_second_layer = ring_diameter + 35
    n_magnet_second_layer = 6
    n_ring_second_layer = 2
    ring_length_second_layer = 90

    # For Optimization
    """
    for x in range(0, 100, 1):
        halbach = magpy.Collection()
        halbach = generate_magnet(halbach)
        plot_result(halbach, print_result=True, arrow=True, plot_show=False)
    """

    halbach = magpy.Collection()

    halbach = generate_magnet(halbach)

    plot_halback(halbach)

    plot_result(halbach, print_result=True, arrow=False)