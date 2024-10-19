import numpy as np
import matplotlib.pyplot as plt

# Parameters
r_axis_size = 13
z_axis_size = 150
r_len = 0.01
z_len = 1

# Size of grid:
dz = z_len / z_axis_size
dr = r_len / r_axis_size
epoch = 1000

# Grid initialization:
grid = np.zeros((2 * r_axis_size, z_axis_size))
grid[:,:] = 20

# Initial boundary conditions:
T_low = 20
T_high = 300
grid[:, 0] = T_low  # Left boundary (z = 0)
grid[:, -1] = T_high  # Right boundary (z = z_len)
grid[0,:30] = T_low   # Upper side of steel tube on 20K side
grid[-1,:30] = T_low  # Lower side of steel tube on 20K side
grid[0,-30:] = T_high # Upper side of steel tube on 300K side
grid[-1,-30:] = T_high # Lower side of steel tube on 300K side
# Inputting initial boundary conditions:
h_he = 50
k_aisi20 = 3.5
k_aisi300 = 16.2
# Heat diffusion simulation
for t in range(epoch):
    for i in range(1, 2 * r_axis_size-1):
        if i < r_axis_size:
            r_m = r_len - i * dr  # Adjust radial position to the center of the cell
        else:
            r_m = (i - r_axis_size) * dr
        for j in range(1, z_axis_size - 1):# Iterate over internal z-axis points
            sum_nom = 0
            # Denominator for the equation
            # Neighbors contributions
            # Handle out-of-bounds explicitly
            if i in [1, 23] and j in [148]:
                nom_1 = h_he * dr * grid[i + 1][j]
                nom_2 = h_he * dz * grid[i][j + 1]
                nom_3 = k_aisi300 * dr * dz * grid[i - 1][j]
                nom_4 = k_aisi300 * dz * dr * grid[i][j - 1]
                # Denominator term
                denom = h_he * dr + h_he * dz + 2 * k_aisi300 * dr * dz
            elif i in [1,23] and j in [1]:
                nom_1 = h_he * dr * grid[i + 1][j]
                nom_2 = h_he * dz * grid[i][j - 1]
                nom_3 = k_aisi20 * dr * dz * grid[i - 1][j]
                nom_4 = k_aisi20 * dz * dr * grid[i][j + 1]
                # Denominator term
                denom = h_he * dr + h_he * dz + 2 * k_aisi300 * dr * dz
            elif i in [2,24] and j in [1]:
                nom_1 = h_he * dr * grid[i - 1][j]
                nom_2 = h_he * dz * grid[i][j - 1]
                nom_3 = k_aisi20 * dr * dz * grid[i + 1][j]
                nom_4 = k_aisi20 * dz * dr * grid[i][j + 1]
                # Denominator term
                denom = h_he * dr + h_he * dz + 2 * k_aisi300 * dr * dz
            elif i in [2, 24] and j in [148]:
                nom_1 = h_he * dr * grid[i - 1][j]
                nom_2 = h_he * dz * grid[i][j + 1]
                nom_3 = k_aisi300 * dr * dz * grid[i - 1][j]
                nom_4 = k_aisi300 * dz * dr * grid[i][j - 1]
                # Denominator term
                denom = h_he * dr + h_he * dz + 2 * k_aisi300 * dr * dz
            elif i in [1, 23] and j in list(range(121,148)) + list(range(0,30)):
                nom_1 = h_he * dr * grid[i + 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi300 * dz * grid[i][j + 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi300 * dr * grid[i - 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = k_aisi300 * dz * grid[i][j - 1]  # k_aisi300 * dz * T_m,n-1
                denom = h_he * dr + k_aisi300 * dz + k_aisi300 * dr + k_aisi300 * dz
            elif i in [1, 23] and j in range(0,30):
                nom_1 = h_he * dr * grid[i + 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi20 * dz * grid[i][j + 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi20 * dr * grid[i - 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = k_aisi20 * dz * grid[i][j - 1]  # k_aisi300 * dz * T_m,n-1
                denom = h_he * dr + k_aisi20 * dz + k_aisi20 * dr + k_aisi20 * dz
            elif i in [2, 24] and j in list(range(121,148)) + list(range(0,30)):
                nom_1 = h_he * dr * grid[i - 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi300 * dz * grid[i][j + 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi300 * dr * grid[i + 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = k_aisi300 * dz * grid[i][j - 1]  # k_aisi300 * dz * T_m,n-1
                denom = h_he * dr + k_aisi300 * dz + k_aisi300 * dr + k_aisi300 * dz
            elif i in [2, 24] and j in range(0,30):
                nom_1 = h_he * dr * grid[i - 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi20 * dz * grid[i][j + 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi20 * dr * grid[i + 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = k_aisi20 * dz * grid[i][j - 1]  # k_aisi300 * dz * T_m,n-1
                denom = h_he * dr + k_aisi20 * dz + k_aisi20 * dr + k_aisi20 * dz
            elif i in [1, 23] and j in [120]:
                nom_1 = h_he * dr * grid[i - 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi300 * dz * grid[i][j + 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi300 * dr * grid[i + 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = 0
                denom = h_he * dr + k_aisi300 * dz + k_aisi300 * dr + k_aisi300 * dz
            elif i in [2, 24] and j in [120]:
                nom_1 = 0
                nom_2 = h_he * dz * grid[i][j + 1]
                nom_3 = k_aisi300 * dr * dz * grid[i - 1][j]
                nom_4 = k_aisi300 * dz * dr * grid[i][j - 1]
                denom = h_he * dr + h_he * dz + 2 * k_aisi300 * dr * dz
            elif i in [1, 23] and j in [30]:
                nom_1 = h_he * dr * grid[i + 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi20 * dz * grid[i][j - 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi20 * dr * grid[i - 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = 0
                denom = h_he * dr + k_aisi20 * dz + k_aisi20 * dr + k_aisi20 * dz
            elif i in [1, 23] and j in [30]:
                nom_1 = h_he * dr * grid[i - 1][j]  # h_he * dr * T_m+1,n
                nom_2 = k_aisi20 * dz * grid[i][j - 1]  # k_aisi300 * dz * T_m,n+1
                nom_3 = k_aisi20 * dr * grid[i + 1][j]  # k_aisi300 * dr * T_m-1,n
                nom_4 = 0
                denom = h_he * dr + k_aisi20 * dz + k_aisi20 * dr + k_aisi20 * dz
            else:
                try:
                    nom_1 = (r_m + dr) * grid[i + 1, j] / dr ** 2
                except IndexError:
                    nom_1 =  (r_m + dr) * grid[i, j] / dr ** 2
                try:
                    nom_2 = (r_m - dr) * grid[i - 1, j] / dr ** 2
                except IndexError:
                    nom_2 = (r_m - dr) * grid[i, j] / dr ** 2
                try:
                    nom_3 = grid[i, j + 1] / dz ** 2
                except IndexError:
                    nom_3 = grid[i, j ] / dz ** 2
                try:
                    nom_4 = grid[i, j - 1] / dz ** 2
                except IndexError:
                    nom_4 = grid[i, j ] / dz ** 2
                denom = 2 * r_m / dr ** 2 + 2 / dz ** 2
                # Summing up the contributions
            sum_nom = nom_1 + nom_2 + nom_3 + nom_4
            grid[i,j] = sum_nom / denom


# Visualization

fig, ax = plt.subplots()

# Display the grid matrix with imshow
cax = ax.imshow(grid, cmap='viridis', aspect='auto')

# Add a colorbar
fig.colorbar(cax, label='Temperature')

# Set title and axis labels
ax.set_title('Heat Diffusion in a Cylindrical Tube')
ax.set_xlabel('Z-axis (Length)')
ax.set_ylabel('R-axis (Radial)')

# Set custom y-ticks for the range -10 mm to 10 mm
yticks = np.linspace(0, grid.shape[0] - 1, 5)  # Positions based on the grid's size
yticklabels = np.linspace(-10, 10, 5)  # Custom labels from -10 mm to 10 mm

# Apply the custom y-ticks and labels
ax.set_yticks(yticks)
ax.set_yticklabels([f'{label:.1f} mm' for label in yticklabels])

# Set custom x-ticks for the range 0 to 1 meter
xticks = np.linspace(0, grid.shape[1] - 1, 5)  # Positions based on the grid's size
xticklabels = np.linspace(0, 1, 5)  # Custom labels from 0 to 1 meter

# Apply the custom x-ticks and labels
ax.set_xticks(xticks)
ax.set_xticklabels([f'{label:.2f} m' for label in xticklabels])
plt.show()