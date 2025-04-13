import numpy as np
from pyscf import gto, dft, geomopt

mol = gto.M(
    atom = '''
    C    0.0000   0.0000   0.0000
    O    1.2000   0.0000   0.0000
    H   -0.5400   0.9400   0.0000
    H   -0.5400  -0.9400   0.0000
    ''',
    basis = '6-31+G(d)',
    unit = 'Angstrom',
    charge = 0,
    spin = 0,
)

mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf.kernel()

# Optimize the geometry using a BFGS algorithm.
mol_opt = geomopt.bfgs_solver(mol, mf)

# Write the optimized geometry to an XYZ file.
coords = mol_opt.atom_coords()
with open("optimized_species.xyz", "w") as f:
    f.write(f"{len(coords)}\n")
    f.write("Optimized geometry using B3LYP/6-31+G(d) with BFGS optimization\n")
    for i, coord in enumerate(coords):
        sym = mol_opt.atom_symbol(i)
        f.write(f"{sym} {coord[0]:.4f} {coord[1]:.4f} {coord[2]:.4f}\n")

print("Geometry optimization complete. Optimized geometry saved to optimized_species.xyz")
