from ase import io
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator, all_properties
from pyscf import gto, dft
import numpy as np
import os

class PySCFCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, xc='B3LYP', basis='6-31+G(d)', charge=0, spin=0, **kwargs):
        super().__init__(**kwargs)
        self.xc = xc
        self.basis = basis
        self.charge = charge
        self.spin = spin

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_properties):
        Calculator.calculate(self, atoms, properties, system_changes)

        atom_list = []
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            atom_list.append(f"{sym} {pos[0]} {pos[1]} {pos[2]}")

        mol = gto.Mole()
        mol.atom = atom_list
        mol.basis = self.basis
        mol.charge = self.charge
        mol.spin = self.spin
        mol.unit = 'Angstrom'
        mol.build()

        mf = dft.RKS(mol)
        mf.xc = self.xc
        mf.conv_check = False
        mf.kernel()

        # Total energy
        self.results['energy'] = mf.e_tot

        # Forces (negative gradient)
        gradients = mf.nuc_grad_method().kernel()
        self.results['forces'] = -gradients

# === Load the molecule from 1.xyz ===
atoms = io.read("mol_1.xyz")

# === Set the PySCF calculator ===
atoms.calc = PySCFCalculator()

# === Run BFGS optimization ===
optimizer = BFGS(atoms, trajectory='opt.traj', logfile='opt.log')
optimizer.run(fmax=0.05)

# === Save the optimized structure ===
io.write("1_opt.xyz", atoms)

