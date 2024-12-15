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

def calculate_distance(atom1, atom2):
    """
    Calculate the distance between two atoms based on their coordinates.
    """
    coords1 = np.array(list(map(float, atom1.split()[2:5])))
    coords2 = np.array(list(map(float, atom2.split()[2:5])))
    return np.linalg.norm(coords1 - coords2)

def create_defect_files(file_path, supercell_size=29.47839554, max_distance=10.0):
    """
    Explore O cores near the central Mg first and move outward. Create defect files
    for unique distances, stopping once a distance exceeds the max_distance.
    """
    full_atom_lines, footer_lines = parse_gulp_file_coordinates(file_path)

    mg_atoms = filter_atom_type_strict(full_atom_lines, "Mg core")
    o_core_atoms = filter_atom_type_strict(full_atom_lines, "O core")
    o_shel_atoms = filter_atom_type_strict(full_atom_lines, "O shel")

    # Find the central Mg core
    fixed_central_mg = find_most_central_atom(mg_atoms, supercell_size)
    print(f"Fixed central Mg core: {fixed_central_mg.strip()}")

    # Calculate distances of all O cores to the central Mg
    distances = [
        (o_core, o_shel, calculate_distance(fixed_central_mg, o_core))
        for o_core, o_shel in zip(o_core_atoms, o_shel_atoms)
    ]

    # Sort the O cores by distance
    distances.sort(key=lambda x: x[2])  # Sort by distance (x[2])

    previous_distances = []  # Track unique distances
    iteration = 1

    for o_core, o_shel, distance in distances:
        # Stop if the distance exceeds the maximum allowed distance
        if distance > max_distance:
            print(f"Stopping: Distance exceeded {max_distance} Å at iteration {iteration}.")
            break

        # Skip if the distance is the same as a previous defect
        if any(np.isclose(distance, d, atol=1e-6) for d in previous_distances):
            print(f"Skipping O core: Distance {distance:.4f} Å already used.")
            continue

        # Add the distance to the list of previously used distances
        previous_distances.append(distance)

        print(f"Iteration {iteration}:")
        print(f"  Fixed Mg core: {fixed_central_mg.strip()}")
        print(f"  Closest O core: {o_core.strip()}")
        print(f"  Closest O shel: {o_shel.strip()}")
        print(f"  Distance between Mg and O core: {distance:.4f} Å")

        # Remove the fixed Mg core, selected O core, and O shel
        mg_atoms.remove(fixed_central_mg)  # Ensure central Mg core is removed
        o_core_atoms.remove(o_core)
        o_shel_atoms.remove(o_shel)

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

        iteration += 1

# Example usage
file_path = "input777.gin"
create_defect_files(file_path, supercell_size=29.47839554, max_distance=10.0)

