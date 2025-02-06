import matplotlib

matplotlib.use('Qt5Agg')  # Define o backend como Qt5Agg

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_unit_cell(filename):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    index = 0
    while not lines[index][0].isdigit() and not lines[index][0] == '-':
        index += 1

    lattice_constant = float(lines[index].split()[0])
    index += 1

    # Unit cell vectors (a1, a2, a3)
    a1 = np.array([float(x) for x in lines[index].split()])
    a2 = np.array([float(x) for x in lines[index + 1].split()])
    a3 = np.array([float(x) for x in lines[index + 2].split()])
    index += 3

    # Número de elementos
    num_elements = int(lines[index])
    index += 1

    atom_data = {}
    for _ in range(num_elements):
        element_name = lines[index].strip()
        index += 1
        atom_count = int(lines[index].split()[0])
        index += 1
        positions = []

        for _ in range(atom_count):
            if index < len(lines):
                pos = np.array([float(x) for x in lines[index].split()[:3]])
                positions.append(pos)
                index += 1

        atom_data[element_name] = np.array(positions)

    return lattice_constant, a1, a2, a3, atom_data


def replicate_atoms(atom_data, a1, a2, a3, nx=3, ny=3, nz=3):
    """Replicate atoms in the unit cell."""
    replicated_data = {}
    for element, positions in atom_data.items():
        replicated_positions = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    translation = i * a1 + j * a2 + k * a3
                    for atom in positions:
                        replicated_positions.append(atom + translation)
        replicated_data[element] = np.array(replicated_positions)
    return replicated_data


def plot_unit_cell(atom_data, a1, a2, a3, view_orientation='001'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colors = {'H': 'lightgray', 'He': 'cyan', 'Li': 'red', 'Be': 'darkgreen', 'B': 'brown', 'C': 'black', 'N': 'blue',
              'O': 'blue', 'F': 'green', 'Ne': 'cyan', 'Na': 'blue', 'Mg': 'green', 'Al': 'gray', 'Si': 'yellow',
              'P': 'orange',
              'S': 'yellow', 'Cl': 'green', 'Ar': 'cyan', 'K': 'purple', 'Ca': 'darkgray', 'Sc': 'pink', 'Ti': 'silver',
              'V': 'gray', 'Cr': 'gray', 'Mn': 'purple', 'Fe': 'brown', 'Co': 'blue', 'Ni': 'green', 'Cu': 'orange',
              'Zn': 'lightblue', 'Ga': 'red', 'Ge': 'darkgray', 'As': 'blue', 'Se': 'brown', 'Br': 'red', 'Kr': 'cyan'}

    for element, positions in atom_data.items():
        color = colors.get(element, 'k')
        ax.scatter(*positions.T, color=color, s=100, label=element)

        for i in range(len(positions)):
            for j in range(len(positions)):
                if i != j:
                    diff = positions[j] - positions[i]
                    if np.allclose(diff[:2], [0, 0]) and diff[2] > 0:
                        ax.plot([positions[i][0], positions[j][0]],
                                [positions[i][1], positions[j][1]],
                                [positions[i][2], positions[j][2]],
                                color='gray', lw=0.5)

    if view_orientation == '001':
        ax.view_init(elev=90, azim=0)
    elif view_orientation == '111':
        ax.view_init(elev=30, azim=45)
    elif view_orientation == '110':
        ax.view_init(elev=60, azim=45)
    else:
        ax.view_init(elev=20, azim=30)

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    plt.legend()

    def set_axes_fixed(ax, length=15):
        """Define os limites fixos dos eixos de 0 a 'length' Å."""
        ax.set_xlim(-0.5, length)
        ax.set_ylim(-0.5, length)
        ax.set_zlim(-0.5, 0.75)

    set_axes_fixed(ax, length=1.2)

    plt.show()


# Main program
filename = "Ga2O3.txt"
lattice_constant, a1, a2, a3, atom_data = read_unit_cell(filename)

replicated_positions = replicate_atoms(atom_data, a1, a2, a3, nx=1, ny=1, nz=1)

view_orientation = '001'
plot_unit_cell(replicated_positions, a1, a2, a3, view_orientation=view_orientation)