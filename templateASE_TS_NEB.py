# neb_ts.py
from ase import io
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator, all_properties
from pyscf import gto, dft
import numpy as np

# Define the PySCF-based calculator (as before) for use in the NEB calculation.
class PySCFCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, xc='B3LYP', basis='6-31+G(d)', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.xc = xc
        self.basis = basis

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_properties):
        Calculator.calculate(self, atoms, properties, system_changes)
        # Create a PySCF molecule from ASE Atoms
        atom_list = []
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            atom_list.append(f"{sym} {pos[0]} {pos[1]} {pos[2]}")
        mol = gto.M(atom="\n".join(atom_list), basis=self.basis, unit='Angstrom', charge=0, spin=0)
        mf = dft.RKS(mol)
        mf.xc = self.xc
        energy = mf.kernel()
        forces = -mf.grad()  # Negative gradient
        self.results['energy'] = energy
        self.results['forces'] = forces

# Load initial and final state structures (ASE Atoms objects)
initial = io.read("initial.xyz")
final = io.read("final.xyz")

# Number of intermediate images you want (excluding initial and final)
n_images = 10
images = [initial]
for i in range(n_images):
    image = initial.copy()
    images.append(image)
images.append(final)

# Attach the PySCF calculator to all images.
for image in images:
    image.calc = PySCFCalculator(xc='B3LYP', basis='6-31+G(d)')

# Set up the NEB
neb = NEB(images)
neb.interpolate()

# Optimize the NEB path using BFGS optimizer
optimizer = BFGS(neb, trajectory='neb.traj')
optimizer.run(fmax=0.05)

# Save each image after optimization
for i, image in enumerate(images):
    io.write(f'neb_image_{i}.xyz', image)
print("NEB optimization complete. NEB images saved as neb_image_*.xyz")
