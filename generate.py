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

def create_defect_files(file_path, max_iterations=5, delta=5.0, supercell_size=29.47839554):
    """
    Create defect files by defining a smaller cube around the most central Mg atom
    and iteratively removing the Mg and its closest O core and O shel atoms.
    """
    # Parse the file to get atom data and footer
    full_atom_lines, footer_lines = parse_gulp_file_coordinates(file_path)

    # Filter Mg core atoms
    mg_atoms = filter_atom_type_strict(full_atom_lines, "Mg core")
    o_core_atoms = filter_atom_type_strict(full_atom_lines, "O core")
    o_shel_atoms = filter_atom_type_strict(full_atom_lines, "O shel")

    # Find the most central Mg core atom
    most_central_mg = find_most_central_atom(mg_atoms)
    mg_coords = np.array(list(map(float, most_central_mg.split()[2:5])))

    # Define the smaller cube boundaries, ensuring it does not exceed the supercell limits
    x_min = mg_coords[0]
    x_max = min(mg_coords[0] + delta, supercell_size)
    y_min = mg_coords[1]
    y_max = min(mg_coords[1] + delta, supercell_size)
    z_min = mg_coords[2]
    z_max = min(mg_coords[2] + delta, supercell_size)

    # Filter atoms within the smaller cube
    def filter_atoms_in_cube(atom_lines):
        filtered_atoms = []
        for line in atom_lines:
            coords = list(map(float, line.split()[2:5]))
            if x_min <= coords[0] < x_max and y_min <= coords[1] < y_max and z_min <= coords[2] < z_max:
                filtered_atoms.append(line)
        return filtered_atoms

    cube_o_core_atoms = filter_atoms_in_cube(o_core_atoms)
    cube_o_shel_atoms = filter_atoms_in_cube(o_shel_atoms)

    # Iteratively create defect files
    for iteration in range(1, max_iterations + 1):
        if not cube_o_core_atoms or not cube_o_shel_atoms:
            print(f"Stopping at iteration {iteration}, not enough oxygen atoms left in the cube.")
            break

        # Find the closest O core and O shel atoms to the most central Mg
        closest_o_core = find_closest_atom(most_central_mg, cube_o_core_atoms)
        closest_o_shel = find_closest_atom(most_central_mg, cube_o_shel_atoms)

        # Calculate the distance between the Mg and the O core
        distance = calculate_distance(most_central_mg, closest_o_core)

        # Remove the Mg, O core, and O shel atoms for this defect file
        updated_atoms = [
            atom for atom in full_atom_lines
            if atom != most_central_mg and atom != closest_o_core and atom != closest_o_shel
        ]

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

        # Save the new file
        output_file = f"defect_{iteration}.gin"
        with open(output_file, "w") as file:
            file.write(header)
            file.writelines(updated_atoms)
            file.writelines(updated_footer)

        print(f"Defect file created: {output_file}")

        # Remove the selected oxygen atoms from the cube lists
        cube_o_core_atoms.remove(closest_o_core)
        cube_o_shel_atoms.remove(closest_o_shel)

# Example usage
file_path = "input777.gin"  # Replace with the actual path to your GULP input file
create_defect_files(file_path, max_iterations=20, delta=5.0, supercell_size=29.47839554)
