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

def create_defect_files(atom_lines, footer_lines, mg_atoms, o_core_atoms, o_shel_atoms, max_iterations=5):
    """
    Create defect files by iteratively removing the most central Mg core atom
    and its closest O core and O shel atoms.
    """
    for iteration in range(1, max_iterations + 1):
        if not mg_atoms or not o_core_atoms or not o_shel_atoms:
            print(f"Stopping at iteration {iteration}, not enough atoms left.")
            break

        # Find the most central Mg core atom
        most_central_mg = find_most_central_atom(mg_atoms)

        # Find the closest O core and O shel atoms
        closest_o_core = find_closest_atom(most_central_mg, o_core_atoms)
        closest_o_shel = find_closest_atom(most_central_mg, o_shel_atoms)

        # Calculate the distance between the Mg and the O core
        distance = calculate_distance(most_central_mg, closest_o_core)

        # Remove the selected atoms from the lists
        mg_atoms.remove(most_central_mg)
        o_core_atoms.remove(closest_o_core)
        o_shel_atoms.remove(closest_o_shel)

        # Create the updated atom list for the defect file
        updated_atoms = [
            atom
            for atom in atom_lines
            if atom not in [most_central_mg, closest_o_core, closest_o_shel]
        ]

        # Add the required header
        header = f"""# 
# Keywords:
# 
opti conp  
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

# Main execution
file_path = "input777.gin"  # Replace with the actual path to your GULP input file

# Parse the file to get atom data and footer
atom_lines_coords, footer_lines_coords = parse_gulp_file_coordinates(file_path)

# Define supercell boundaries
L = 30.0  # Replace with the actual size of your supercell
x_min, x_max = L / 2, L  # Adjust for the octant you want to explore
y_min, y_max = L / 2, L
z_min, z_max = L / 2, L

# Filter atoms in the specified octant
filtered_atom_lines = filter_atoms_in_octant(atom_lines_coords, x_min, x_max, y_min, y_max, z_min, z_max)

# Filter atom types within the selected octant
mg_core_atoms = filter_atom_type_strict(filtered_atom_lines, "Mg core")
o_core_atoms = filter_atom_type_strict(filtered_atom_lines, "O core")
o_shel_atoms = filter_atom_type_strict(filtered_atom_lines, "O shel")

# Create defect files
create_defect_files(filtered_atom_lines, footer_lines_coords, mg_core_atoms, o_core_atoms, o_shel_atoms, max_iterations=20)

