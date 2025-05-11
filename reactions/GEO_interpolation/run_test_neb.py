# neb_ts_from_xyz.py
from ase import io
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator, all_properties
from pyscf import gto, dft
import numpy as np

# Define a PySCF-based calculator
class PySCFCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, xc='B3LYP', basis='6-31+G(d)', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.xc = xc
        self.basis = basis

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_properties):
        Calculator.calculate(self, atoms, properties, system_changes)
        atom_list = [
            f"{sym} {x:.10f} {y:.10f} {z:.10f}"
            for sym, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions())
        ]
        mol = gto.M(atom="\n".join(atom_list), basis=self.basis, unit='Angstrom', charge=0, spin=0)
        mf = dft.RKS(mol)
        mf.xc = self.xc
        energy = mf.kernel()
        forces = -mf.nuc_grad_method().kernel()
        self.results['energy'] = energy
        self.results['forces'] = forces

# === Load interpolated images ===
images = io.read('/home/jkl6090/reactions/interpolation/rxn_1=2_ini_fin_interpolation.xyz', index=':')

# === Attach PySCF calculator to all images ===
for image in images:
    image.calc = PySCFCalculator(xc='B3LYP', basis='6-31+G(d)')

# === Set up NEB ===
neb = NEB(images)

# === Optimize the NEB path ===
optimizer = BFGS(neb, trajectory='neb.traj')
optimizer.run(fmax=0.05)

# === Save optimized images ===
for i, image in enumerate(images):
    io.write(f'neb_image_{i}.xyz', image)

print("NEB optimization complete. NEB images saved as neb_image_*.xyz")

