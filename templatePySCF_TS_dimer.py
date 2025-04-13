import numpy as np
from pyscf import gto, dft, geomopt

mol_init = gto.M(
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

mf_init = dft.RKS(mol_init)
mf_init.xc = 'B3LYP'
mf_init.kernel()

mol_final = gto.M(
    atom = '''
    C    0.1000   0.1000   0.0000
    O    1.3000   0.1000   0.0000
    H   -0.4400   1.0400   0.0000
    H   -0.4400  -0.8400   0.0000
    ''',
    basis = '6-31+G(d)',
    unit = 'Angstrom',
    charge = 0,
    spin = 0,
)

mf_final = dft.RKS(mol_final)
mf_final.xc = 'B3LYP'
mf_final.kernel()

# Use the dimer method to search for the transition state.
# The dimer_solver function takes an initial geometry, its calculator (mf_init),
# and a guess for the TS (here, provided by mol_final).
# You may need to adjust the tolerance and maximum steps.
ts_mol = geomopt.dimer_solver(mol_init, mf_init, guess=mol_final, tol=1e-3, max_steps=100)

# Write the transition state geometry to an XYZ file.
coords = ts_mol.atom_coords()
with open("transition_state.xyz", "w") as f:
    f.write(f"{len(coords)}\n")
    f.write("Transition state geometry using B3LYP/sto-3g and the dimer method\n")
    for i, coord in enumerate(coords):
        sym = ts_mol.atom_symbol(i)
        f.write(f"{sym} {coord[0]:.4f} {coord[1]:.4f} {coord[2]:.4f}\n")

print("Transition state search complete. TS geometry saved to transition_state.xyz")
