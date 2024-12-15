import time
#read the file in


#define the centre of the supercell

# the dimensions of the supercell are 12.63359910 in all directions

#define lattice parameters
a = b = c =   12.63359910 

centre_a = centre_b = centre_c = 12.63359910/2 

#define centre of mass coordinates 

centre = [centre_a, centre_b, centre_c]

file = open("input777.gin")

supercell = file.readlines()[13:4129]

parsed_data = []

for line in supercell:
    row = line.split()                 
    if row[1] == 'shel':
        atom = row [0] + ' shell'
        atom_coords = row[2:5]
        atom_charge = row[5]
        atom_spring_constant = row[6]
        atom_shell_mass = row[7]
    else:
        atom  = row[0] + ' core'                    
        atom_coords = row[2:5]   
        atom_charge = row[5]
        atom_spring_constant =row[6]
        atom_shell_mass = row[7]
    parsed_data.append({"type":atom , "coords" : atom_coords, "charge" : atom_charge,"spring": atom_spring_constant,"shell mass" : atom_shell_mass})


CENTRE_X = 14.73919776  # x,y and z coords of the most central Mg 
CENTRE_Y = 12.63360034 
CENTRE_Z = 14.73919776

centre_array = [CENTRE_X, CENTRE_Y , CENTRE_Z]
for atom in parsed_data:
    atom['coords'] = [float(coord) for coord in atom['coords']]

#for index, atom in enumerate(parsed_data):
 #   if atom['coords'] == centre_array:
  #      print("Central Mg found")
   #     time.sleep(2)
    #    print("Removing central Mg")




def locate_close_oxygen(data,mg_coords):
    closest = None
    min_distance = 1000 #pick a sufficiently large number

    for atom in data:
        if atom["type"] == 'O core':
            
