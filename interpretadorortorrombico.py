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

    lattice_constant = float(lines[index].split()[0]) * 0.5292  # Conversão para Å
    index += 1

    a1 = np.array([float(x) for x in lines[index].split()]) * lattice_constant
    a2 = np.array([float(x) for x in lines[index + 1].split()]) * lattice_constant
    a3 = np.array([float(x) for x in lines[index + 2].split()]) * lattice_constant
    index += 3

    # Ajuste do vetor a3 para estrutura monoclínica (β = 103.7°)
    beta = np.radians(103.7)
    a3 = np.array([a3[0] * np.cos(beta), a3[1], a3[2] * np.sin(beta)])

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
                pos = np.array([float(x) for x in lines[index].split()[:3]]) * lattice_constant
                positions.append(pos)
                index += 1

        atom_data[element_name] = np.array(positions)

    return a1, a2, a3, atom_data


def replicate_atoms(atom_data, a1, a2, a3, nx=3, ny=3, nz=3):
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


def plot_unit_cell(atom_data, a1, a2, a3):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    colors = {'Ga': 'red', 'O': 'blue'}

    for element, positions in atom_data.items():
        color = colors.get(element, 'k')
        ax.scatter(*positions.T, color=color, s=100, label=element)

    ax.set_xlabel('X-axis (Å)')
    ax.set_ylabel('Y-axis (Å)')
    ax.set_zlabel('Z-axis (Å)')
    ax.view_init(elev=20, azim=30)
    plt.legend()
    plt.show()


# Main program
filename = "Ga2O3.txt"  # Arquivo contendo os dados do Ga2O3
a1, a2, a3, atom_data = read_unit_cell(filename)

replicated_positions = replicate_atoms(atom_data, a1, a2, a3, nx=2, ny=2, nz=1)
plot_unit_cell(replicated_positions, a1, a2, a3)
