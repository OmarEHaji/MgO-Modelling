import numpy as np

def parse_gulp_file_coordinates(file_path):
    """
    Parse the GULP input file to extract atomic coordinates (Mg core, O core, O shel)
    and ignore non-coordinate lines.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_lines = []
    footer_lines = lines[-13:]  # Last 13 lines to be replaced

    for line in lines:
        # Split the line into components
        parts = line.split()
        if len(parts) > 4 and all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in parts[2:5]):
            atom_lines.append(line)

    return atom_lines, footer_lines


def filter_atom_type_strict(atom_lines, atom_type):
    """
    Strictly filter atom lines by a specific type (e.g., 'Mg core', 'O core', 'O shel').
    """
    return [line for line in atom_lines if atom_type == f"{line.split()[0]} {line.split()[1]}"]


def filter_atoms_in_octant(atom_lines, x_min, x_max, y_min, y_max, z_min, z_max):
    """
    Filter atoms to include only those within the specified octant boundaries.
    """
    filtered_atoms = []
    for line in atom_lines:
        coords = list(map(float, line.split()[2:5]))
        if x_min <= coords[0] < x_max and y_min <= coords[1] < y_max and z_min <= coords[2] < z_max:
            filtered_atoms.append(line)
    return filtered_atoms


def find_most_central_atom(atoms):
    """
    Find the most central atom based on coordinates.
    """
    coords = [list(map(float, atom.split()[2:5])) for atom in atoms]
    distances = [np.linalg.norm(coord) for coord in coords]
    central_index = np.argmin(distances)
    return atoms[central_index]


def find_closest_atom(reference_atom, atoms):
    """
    Find the closest atom to a reference atom.
    """
    ref_coords = np.array(list(map(float, reference_atom.split()[2:5])))
    atom_coords = [list(map(float, atom.split()[2:5])) for atom in atoms]
    distances = [np.linalg.norm(ref_coords - np.array(coord)) for coord in atom_coords]
    closest_index = np.argmin(distances)
    return atoms[closest_index]


def calculate_distance(atom1, atom2):
    """
    Calculate the distance between two atoms based on their coordinates.
    """
    coords1 = np.array(list(map(float, atom1.split()[2:5])))
    coords2 = np.array(list(map(float, atom2.split()[2:5])))
    return np.linalg.norm(coords1 - coords2)


def create_defect_files(full_atom_lines, footer_lines, octant_mg_atoms, octant_o_core_atoms, octant_o_shel_atoms, max_iterations=5):
    """
    Create defect files by iteratively removing the most central Mg core atom
    and its closest O core and O shel atoms, ensuring the full supercell is output.
    """
    for iteration in range(1, max_iterations + 1):
        if not octant_mg_atoms or not octant_o_core_atoms or not octant_o_shel_atoms:
            print(f"Stopping at iteration {iteration}, not enough atoms left in the octant.")
            break

        # Find the most central Mg core atom in the octant
        most_central_mg = find_most_central_atom(octant_mg_atoms)

        # Find the closest O core atom
        closest_o_core = find_closest_atom(most_central_mg, octant_o_core_atoms)

        # Find the corresponding O shel atom with the same coordinates as the O core
        closest_o_shel = None
        for atom in octant_o_shel_atoms:
            if np.allclose(list(map(float, atom.split()[2:5])), list(map(float, closest_o_core.split()[2:5]))):
                closest_o_shel = atom
                break

        if closest_o_shel is None:
            print(f"Error: Could not find matching O shel for O core {closest_o_core}.")
            break

        # Debugging: Print the atoms being removed
        print(f"Iteration {iteration}: Removing Mg {most_central_mg}, O core {closest_o_core}, O shel {closest_o_shel}")

        # Calculate the distance between the Mg and the O core
        distance = calculate_distance(most_central_mg, closest_o_core)

        # Remove only the selected atoms from the full atom list
        updated_atoms = []
        for atom in full_atom_lines:
            atom_type = f"{atom.split()[0]} {atom.split()[1]}"
            coords = list(map(float, atom.split()[2:5]))

            if (atom_type == "Mg core" and coords == list(map(float, most_central_mg.split()[2:5]))) or \
               (atom_type == "O core" and coords == list(map(float, closest_o_core.split()[2:5]))) or \
               (atom_type == "O shel" and coords == list(map(float, closest_o_shel.split()[2:5]))):
                continue  # Skip the selected atoms
            updated_atoms.append(atom)

        # Add the required header
        header = f"""# 
# Keywords:
# 
energy  
# 
# Options:
# 
title
"Distance between Mg and O removed: {distance:.4f} Å"
end
cell
  29.47839554  29.47839892  29.47839554  90.00000000  90.00000000  90.00000000
cartesian region  1
"""

        # Update the footer to include the library reference
        updated_footer = ["lib exercise.lib\n"]

        # Save the new file
        output_file = f"defect_{iteration}.gin"
        with open(output_file, "w") as file:
            file.write(header)
            file.writelines(updated_atoms)
            file.writelines(updated_footer)
        
        print(f"Defect file created: {output_file}")

        # Remove the selected atoms from the octant-specific lists
        octant_mg_atoms.remove(most_central_mg)
        octant_o_core_atoms.remove(closest_o_core)
        octant_o_shel_atoms.remove(closest_o_shel)


# Main execution
file_path = "input777.gin"  # Replace with the actual path to your GULP input file

# Parse the file to get atom data and footer
atom_lines_coords, footer_lines_coords = parse_gulp_file_coordinates(file_path)

# Define supercell boundaries
L = 29.47839554  # Replace with the actual size of your supercell
x_min, x_max = L / 2, L  # Adjust for the octant you want to explore
y_min, y_max = L / 2, L
z_min, z_max = L / 2, L

# Filter atoms in the specified octant
octant_atom_lines = filter_atoms_in_octant(atom_lines_coords, x_min, x_max, y_min, y_max, z_min, z_max)

# Filter atom types within the selected octant
octant_mg_atoms = filter_atom_type_strict(octant_atom_lines, "Mg core")
octant_o_core_atoms = filter_atom_type_strict(octant_atom_lines, "O core")
octant_o_shel_atoms = filter_atom_type_strict(octant_atom_lines, "O shel")

# Debugging: Print the counts of filtered atoms
print(f"Total atoms in octant: {len(octant_atom_lines)}")
print(f"Mg core atoms: {len(octant_mg_atoms)}")
print(f"O core atoms: {len(octant_o_core_atoms)}")
print(f"O shel atoms: {len(octant_o_shel_atoms)}")

# Check if there are enough atoms to proceed
if not octant_mg_atoms or not octant_o_core_atoms or not octant_o_shel_atoms:
    print("No atoms found in the octant or insufficient atoms for defect creation.")
else:
    # Create defect files
    create_defect_files(atom_lines_coords, footer_lines_coords, octant_mg_atoms, octant_o_core_atoms, octant_o_shel_atoms, max_iterations=20)

