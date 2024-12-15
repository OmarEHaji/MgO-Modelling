import numpy as np

def parse_gulp_file_coordinates(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    atom_lines = []
    footer_lines = lines[-13:]  # Last 13 lines to be replaced

    for line in lines:
        parts = line.split()
        if len(parts) > 4 and all(part.replace('.', '', 1).isdigit() or part.lstrip('-').replace('.', '', 1).isdigit() for part in parts[2:5]):
            atom_lines.append(line)

    return atom_lines, footer_lines

def filter_atom_type_strict(atom_lines, atom_type):
    return [line for line in atom_lines if atom_type == f"{line.split()[0]} {line.split()[1]}"]

def find_most_central_atom(atoms, supercell_size):
    """
    Find the atom closest to the center of the full supercell.
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

def calculate_displacement_vector(atom1, atom2):
    """
    Calculate the displacement vector between two atoms.
    """
    coords1 = np.array(list(map(float, atom1.split()[2:5])))
    coords2 = np.array(list(map(float, atom2.split()[2:5])))
    return coords2 - coords1

def is_vector_unique(new_vector, previous_vectors, tol=1e-6):
    """
    Check if the new displacement vector is unique compared to previous vectors.
    """
    for vec in previous_vectors:
        if np.allclose(new_vector, vec, atol=tol):
            return False
    return True

def create_defect_files(file_path, max_iterations=5, supercell_size=29.47839554):
    """
    Create defect files by removing the central Mg core and the closest O core and O shel atoms,
    ensuring symmetry uniqueness of the defects.
    """
    full_atom_lines, footer_lines = parse_gulp_file_coordinates(file_path)

    mg_atoms = filter_atom_type_strict(full_atom_lines, "Mg core")
    o_core_atoms = filter_atom_type_strict(full_atom_lines, "O core")
    o_shel_atoms = filter_atom_type_strict(full_atom_lines, "O shel")

    # Find the central Mg core
    fixed_central_mg = find_most_central_atom(mg_atoms, supercell_size)
    mg_coords = np.array(list(map(float, fixed_central_mg.split()[2:5])))

    print(f"Fixed central Mg core: {fixed_central_mg.strip()} with coordinates: {mg_coords}")

    previous_vectors = []

    for iteration in range(1, max_iterations + 1):
        if not o_core_atoms or not o_shel_atoms:
            print(f"Stopping at iteration {iteration}: Not enough oxygen atoms left.")
            break

        # Find the closest O core to the fixed central Mg core
        closest_o_core = find_closest_atom(fixed_central_mg, o_core_atoms)
        core_coords = np.array(list(map(float, closest_o_core.split()[2:5])))

        # Calculate the displacement vector
        displacement_vector = calculate_displacement_vector(fixed_central_mg, closest_o_core)

        # Ensure the vector is unique
        if not is_vector_unique(displacement_vector, previous_vectors):
            print(f"Skipping O core at iteration {iteration}: Not symmetry unique.")
            o_core_atoms.remove(closest_o_core)
            continue

        # Find the corresponding O shel by index
        core_index = o_core_atoms.index(closest_o_core)
        closest_o_shel = o_shel_atoms[core_index]

        # Calculate the distance
        distance = calculate_distance(fixed_central_mg, closest_o_core)

        print(f"Iteration {iteration}:")
        print(f"  Fixed Mg core: {fixed_central_mg.strip()}")
        print(f"  Closest O core: {closest_o_core.strip()}")
        print(f"  Closest O shel: {closest_o_shel.strip()}")
        print(f"  Distance between Mg and O core: {distance:.4f} Å")

        # Add the displacement vector to the list of previous vectors
        previous_vectors.append(displacement_vector)

        # Remove the selected atoms
        o_core_atoms.pop(core_index)
        o_shel_atoms.pop(core_index)

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

        # Write the new defect file
        output_file = f"defect_{iteration}.gin"
        with open(output_file, "w") as file:
            file.write(header)
            file.writelines(mg_atoms + o_core_atoms + o_shel_atoms + footer_lines)

        print(f"Defect file created: {output_file}")

# Example usage
file_path = "input777.gin"
create_defect_files(file_path, max_iterations=20, supercell_size=29.47839554)

