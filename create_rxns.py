import os
import numpy as np

def read_xyz(filename):
    """
    Reads an XYZ file.
    Returns a list of tuples (element, np.array([x, y, z])).
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    try:
        natoms = int(lines[0].strip())
    except ValueError:
        raise ValueError("The first line should be the number of atoms.")
    
    # Skip the comment (line 2) and then read atoms

    atoms = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        element = parts[0]
        x, y, z = map(float, parts[1:4])
        atoms.append((element, np.array([x, y, z])))
    return atoms

def write_xyz(filename, atoms, comment=""):
    """
    Writes the list of atoms (element, np.array([x,y,z])) to an XYZ file.
    """
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{comment}\n")
        for element, coord in atoms:
            f.write("{:s} {: .5f} {: .5f} {: .5f}\n".format(element, coord[0], coord[1], coord[2]))

def move_to_origin(atoms):
    """
    Recenters a molecule so that its geometric center is at the origin.
    """
    coords = np.array([coord for _, coord in atoms])
    center = coords.mean(axis=0)
    new_atoms = [(elem, coord - center) for elem, coord in atoms]
    return new_atoms

def combine_molecules(molecules, spacing_vector=np.array([10.0, 0.0, 0.0])):
    """
    Combines a list of molecules into one. Each molecule is offset by a multiple of
    the given spacing vector. Then the entire combined molecule is recentered.
    
    Parameters:
      molecules - list of molecule data (each a list of (element, np.array) tuples)
      spacing - tuple (dx, dy, dz) specifying the displacement between each molecule.
    """
    combined = []
    for i, mol in enumerate(molecules):
        offset = spacing_vector * i
        for elem, coord in mol:
            combined.append((elem, coord + offset))
    combined_centered = move_to_origin(combined)
    return combined_centered

def process_reactions(rxns, species_dir="./", write_to_dir="./"):
    """
    Processes a list of reactions defined by the rxns list.
    Each reaction is defined as: [reactants, products], where reactants and products
    are lists of species IDs.
    
    For each reaction, the reactant and product molecules are combined and saved 
    into an 'ini' file and a 'fin' file. The file naming convention is:
    
      rxn_{reactant_species}={product_species}_ini.xyz   and 
      rxn_{reactant_species}={product_species}_fin.xyz
    """
    for rxn in rxns:
        if len(rxn) == 3: 
            reactants, products, spacing_vector = rxn
        else: 
            reactants, products = rxn
            spacing_vector = np.array([10.0, 0.0, 0.0])

        # Process reactants
        reactant_mols = []
        for species in reactants:
            fname = f"{species_dir}mol_{species}.xyz"
            if not os.path.exists(fname):
                raise FileNotFoundError(f"File {fname} not found.")
            mol = read_xyz(fname)
            mol_centered = move_to_origin(mol)
            reactant_mols.append(mol_centered)
        combined_reactants = combine_molecules(reactant_mols, spacing_vector)
        
        # Process products
        product_mols = []
        for species in products:
            fname = f"{species_dir}mol_{species}.xyz"
            if not os.path.exists(fname):
                raise FileNotFoundError(f"File {fname} not found.")
            mol = read_xyz(fname)
            mol_centered = move_to_origin(mol)
            product_mols.append(mol_centered)
        combined_products = combine_molecules(product_mols, spacing_vector)
        
        # Create a string for the file name based on product species IDs (joined with '+')
        reactants_name_str = "+".join(str(s) for s in reactants)
        products_name_str = "+".join(str(s) for s in products)
        rxn_ini_filename = f"rxn_{reactants_name_str}={products_name_str}_ini.xyz"
        rxn_fin_filename = f"rxn_{reactants_name_str}={products_name_str}_fin.xyz"
        
        write_xyz(f"{write_to_dir}{rxn_ini_filename}", combined_reactants,
                  comment=f"Reaction {reactants_name_str}={products_name_str} initial: reactants {reactants}")
        write_xyz(f"{write_to_dir}{rxn_fin_filename}", combined_products,
                  comment=f"Reaction {reactants_name_str}={products_name_str} final: products {products}")
        
        print(f"Processed reaction {reactants_name_str}={products_name_str}")

if __name__ == "__main__":
    # Each reaction is of the form: [ [reactant species IDs], [product species IDs], (optional shift between reactants and products) ]
    rxns = [
        [[1], [2]],
        [[1], [3, 4]],
        [[1], [5, 6], np.array([-10, 0, 0])], 
        [[1], [8, 13], np.array([-10, 0, 0])], 
        [[1], [10, 11, 12]],
        [[1], [10, 25]],
        [[2], [3, 4], np.array([-10, 10, 0])],
        [[2], [10, 16]], 
        [[2], [11, 18]], 
        [[2], [17]],
        [[2], [20, 21]],
        [[3], [7, 10]],
        [[3], [20, 23]],
        [[4], [7, 22]],
        [[4], [24]],
        [[5], [7, 8]],
        [[5], [9]],
        [[13], [14]],
        [[13], [15]],
        [[25], [4, 7]],
        [[25], [11, 12]],
        [[25], [26]], 
        [[25], [27]]
    ]
    process_reactions(rxns, species_dir="./molecules/", write_to_dir="./reactions/")
