# Full updated script with corrected core-shell indexing based on the GULP file structure.

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

def find_most_central_atom(atoms, supercell_size):
    """
    Find the most central atom based on coordinates relative to the true center of the supercell.
    """
    center = np.array([supercell_size / 2] * 3)
    coords = [list(map(float, atom.split()[2:5])) for atom in atoms]
    distances = [np.linalg.norm(np.array(coord) - center) for coord in coords]
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

def filter_atoms_with_indices(atom_lines, x_min, x_max, y_min, y_max, z_min, z_max):
    """
    Filter atoms within a specific cube and return both the filtered atoms and their original indices.
    """
    filtered_atoms = []
    indices = []
    for i, line in enumerate(atom_lines):
        coords = list(map(float, line.split()[2:5]))
        if x_min <= coords[0] < x_max and y_min <= coords[1] < y_max and z_min <= coords[2] < z_max:
            filtered_atoms.append(line)
            indices.append(i)
    return filtered_atoms, indices

def create_defect_files(file_path, max_iterations=5, delta=6.0, supercell_size=29.47839554):
    """
    Create defect files by defining a smaller cube around the most central Mg atom
    and iteratively removing the Mg and its closest O core and O shel atoms.
    """
    full_atom_lines, footer_lines = parse_gulp_file_coordinates(file_path)

    mg_atoms = filter_atom_type_strict(full_atom_lines, "Mg core")
    o_core_atoms = filter_atom_type_strict(full_atom_lines, "O core")
    o_shel_atoms = filter_atom_type_strict(full_atom_lines, "O shel")

    most_central_mg = find_most_central_atom(mg_atoms, supercell_size)
    mg_coords = np.array(list(map(float, most_central_mg.split()[2:5])))

    print(f"Most central Mg core: {most_central_mg.strip()} with coordinates: {mg_coords}")

    for iteration in range(1, max_iterations + 1):
        # Define the smaller cube boundaries
        x_min, x_max = mg_coords[0], min(mg_coords[0] + delta, supercell_size)
        y_min, y_max = mg_coords[1], min(mg_coords[1] + delta, supercell_size)
        z_min, z_max = mg_coords[2], min(mg_coords[2] + delta, supercell_size)

        # Filter atoms within the cube and maintain index alignment
        cube_o_core_atoms, core_indices = filter_atoms_with_indices(o_core_atoms, x_min, x_max, y_min, y_max, z_min, z_max)
        cube_o_shel_atoms = [o_shel_atoms[i] for i in core_indices]

        if not cube_o_core_atoms or not cube_o_shel_atoms:
            print(f"Stopping at iteration {iteration}: Not enough atoms left in the cube.")
            break

        # Find the closest O core
        closest_o_core = find_closest_atom(most_central_mg, cube_o_core_atoms)

        # Find the corresponding O shel by index
        core_index = cube_o_core_atoms.index(closest_o_core)
        closest_o_shel = cube_o_shel_atoms[core_index]

        distance = calculate_distance(most_central_mg, closest_o_core)

        print(f"Iteration {iteration}:")
        print(f"  Mg core: {most_central_mg.strip()}")
        print(f"  Closest O core: {closest_o_core.strip()}")
        print(f"  Closest O shel: {closest_o_shel.strip()}")
        print(f"  Distance between Mg and O core: {distance:.4f} Å")

        # Remove the selected atoms from the global lists
        o_core_atoms.pop(core_indices[core_index])
        o_shel_atoms.pop(core_indices[core_index])
        mg_atoms.remove(most_central_mg)

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
  {supercell_size}  {supercell_size}  {supercell_size}  90.00000000  90.00000000  90.00000000
cartesian region  1
"""

        # Update the footer to include the library reference
        updated_footer = ["lib exercise.lib\n"]

        # Write the new defect file
        output_file = f"defect_{iteration}.gin"
        with open(output_file, "w") as file:
            file.write(header)
            file.writelines(mg_atoms + o_core_atoms + o_shel_atoms + footer_lines)
            file.writelines(updated_footer)

        print(f"Defect file created: {output_file}")

        # Update the atom pool and recalculate the most central Mg core
        if not mg_atoms:
            print("No more Mg core atoms left.")
            break
        most_central_mg = find_most_central_atom(mg_atoms, supercell_size)
        mg_coords = np.array(list(map(float, most_central_mg.split()[2:5])))

# Example usage
file_path = "input777.gin"  # Replace with the actual path to your GULP input file
create_defect_files(file_path, max_iterations=20, delta=6.0, supercell_size=29.47839554)


