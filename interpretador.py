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

    # Atomic positions
    atom_count = int(lines[7].strip().split()[0])
    atomic_positions = []
    for i in range(8, 8 + atom_count):
        pos = np.array([float(x) for x in lines[i].strip().split()])
        atomic_positions.append(pos)

    return lattice_constant, a1, a2, a3, atomic_positions

def replicate_atoms(atomic_positions, a1, a2, a3, nx=3, ny=3, nz=3):
    """Replicate atoms in the unit cell."""
    replicated_positions = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                # Apply translation for each replica
                translation = i * a1 + j * a2 + k * a3
                for atom in atomic_positions:
                    replicated_positions.append(atom + translation)
    return np.array(replicated_positions)

def plot_unit_cell(atomic_positions, a1, a2, a3, distance_threshold=1.5, view_orientation='001'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot atomic positions as spheres
    ax.scatter(*np.array(atomic_positions).T, color='r', s=100, label='Atoms')  # Red for atoms

    # Draw lines between atoms that are close enough (within distance threshold)
    for i in range(len(atomic_positions)):
        for j in range(i + 1, len(atomic_positions)):
            dist = np.linalg.norm(atomic_positions[i] - atomic_positions[j])
            if dist < distance_threshold:
                ax.plot([atomic_positions[i][0], atomic_positions[j][0]],
                        [atomic_positions[i][1], atomic_positions[j][1]],
                        [atomic_positions[i][2], atomic_positions[j][2]], color='b', lw=0.5)

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
    plt.show()  # Ensure interactivity

# Main program
filename = "bcc_001.txt"  # Replace with your file name
lattice_constant, a1, a2, a3, atomic_positions = read_unit_cell(filename)

# Replicate atoms to form a larger structure
replicated_positions = replicate_atoms(atomic_positions, a1, a2, a3, nx=2, ny=2, nz=1)


view_orientation = '111'
plot_unit_cell(replicated_positions, a1, a2, a3, distance_threshold=1.5, view_orientation=view_orientation)
