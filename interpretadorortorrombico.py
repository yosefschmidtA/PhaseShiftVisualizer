import matplotlib

matplotlib.use('Qt5Agg')  # Define o backend como Qt5Agg

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def read_unit_cell(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Lattice constant
    lattice_constant = float(lines[1].strip().split()[0])

    # Unit cell vectors (a1, a2, a3)
    a1 = np.array([float(x) for x in lines[2].strip().split()])
    a2 = np.array([float(x) for x in lines[3].strip().split()])
    a3 = np.array([float(x) for x in lines[4].strip().split()])

    # Read number of elements
    index = 5
    num_elements = int(lines[index].strip())
    index += 1

    atom_data = {}

    for _ in range(num_elements):
        element_info = lines[index].strip().split()
        element_name = element_info[1]  # Get element symbol
        index += 1
        atom_count = int(lines[index].strip().split()[0])
        index += 1
        positions = []

        for _ in range(atom_count):
            pos = np.array([float(x) for x in lines[index].strip().split()[:3]])
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
                    # Apply translation for each replica
                    translation = i * a1 + j * a2 + k * a3
                    for atom in positions:
                        replicated_positions.append(atom + translation)
        replicated_data[element] = np.array(replicated_positions)
    return replicated_data

def plot_unit_cell(atom_data, a1, a2, a3, view_orientation='001'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Color mapping for elements
    colors = {
        'H': 'lightgray', 'He': 'cyan', 'Li': 'red', 'Be': 'darkgreen', 'B': 'brown', 'C': 'black', 'N': 'blue',
        'O': 'red', 'F': 'green', 'Ne': 'cyan', 'Na': 'blue', 'Mg': 'green', 'Al': 'gray', 'Si': 'yellow', 'P': 'orange',
        'S': 'yellow', 'Cl': 'green', 'Ar': 'cyan', 'K': 'purple', 'Ca': 'darkgray', 'Sc': 'pink', 'Ti': 'silver',
        'V': 'gray', 'Cr': 'gray', 'Mn': 'purple', 'Fe': 'brown', 'Co': 'blue', 'Ni': 'green', 'Cu': 'orange',
        'Zn': 'lightblue', 'Ga': 'red', 'Ge': 'darkgray', 'As': 'blue', 'Se': 'brown', 'Br': 'red', 'Kr': 'cyan'
    }

    # Plot atomic positions as spheres
    for element, positions in atom_data.items():
        color = colors.get(element, 'k')  # Default to black if element not listed
        ax.scatter(*positions.T, color=color, s=100, label=element)

        # Connect atoms only along the vertical direction (z-axis)
        for i in range(len(positions)):
            for j in range(len(positions)):
                if i != j:
                    diff = positions[j] - positions[i]
                    if np.allclose(diff[:2], [0, 0]) and diff[2] > 0:  # Only connect lower to upper
                        ax.plot([positions[i][0], positions[j][0]],
                                [positions[i][1], positions[j][1]],
                                [positions[i][2], positions[j][2]],
                                color='gray', lw=0.5)

    # Set camera view based on the specified orientation
    if view_orientation == '001':
        ax.view_init(elev=90, azim=0)
    elif view_orientation == '111':
        ax.view_init(elev=30, azim=45)
    elif view_orientation == '110':
        ax.view_init(elev=60, azim=45)
    else:
        ax.view_init(elev=20, azim=30)  # Default view if not recognized

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Z-axis')
    plt.legend()
    plt.show()

# Main program
filename = "arquivomodelo.txt"  # Replace with your file name
lattice_constant, a1, a2, a3, atom_data = read_unit_cell(filename)

# Replicate atoms to form a larger structure
replicated_positions = replicate_atoms(atom_data, a1, a2, a3, nx=2, ny=2, nz=1)

view_orientation = '001'
plot_unit_cell(replicated_positions, a1, a2, a3, view_orientation=view_orientation)


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib

matplotlib.use('Qt5Agg')

def read_unit_cell(filename):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    index = 0
    while not lines[index][0].isdigit() and not lines[index][0] == '-':
        index += 1

    lattice_constant = float(lines[index].split()[0])
    index += 1

    a1 = np.array([float(x) for x in lines[index].split()])
    a2 = np.array([float(x) for x in lines[index + 1].split()])
    a3 = np.array([float(x) for x in lines[index + 2].split()])
    index += 3

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


def plot_unit_cell(a1, a2, a3, atom_data):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colors = {'Ba': 'red', 'Ti': 'blue','O': 'green'}

    # Desenhando a célula unitária
    origin = np.array([0, 0, 0])
    vertices = [
        origin, origin + a1, origin + a2, origin + a1 + a2,
                origin + a3, origin + a1 + a3, origin + a2 + a3, origin + a1 + a2 + a3
    ]
    faces = [
        [vertices[0], vertices[1], vertices[5], vertices[4]],
        [vertices[1], vertices[3], vertices[7], vertices[5]],
        [vertices[3], vertices[2], vertices[6], vertices[7]],
        [vertices[2], vertices[0], vertices[4], vertices[6]],
        [vertices[4], vertices[5], vertices[7], vertices[6]],
        [vertices[0], vertices[1], vertices[3], vertices[2]]
    ]
    ax.add_collection3d(Poly3DCollection(faces, alpha=0.3, edgecolor='k'))

    # Plotando os átomos
    for element, positions in atom_data.items():
        color = colors.get(element, 'black')
        ax.scatter(*positions.T, color=color, s=100, label=element)
        for pos in positions:
            ax.text(pos[0], pos[1], pos[2], element, color=color, fontsize=10)

    # Ajustes dos eixos
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([np.linalg.norm(a1), np.linalg.norm(a2), np.linalg.norm(a3)])
    plt.legend()
    plt.show()


# Executar o código
filename = "arquivomodelo.txt"
lattice_constant, a1, a2, a3, atom_data = read_unit_cell(filename)
plot_unit_cell(a1, a2, a3, atom_data)
