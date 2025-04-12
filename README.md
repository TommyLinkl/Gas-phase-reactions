# Gas-phase Reaction Network Example

This directory provides a starting point for building and testing our framework for accelerated transition state (TS) search and the autonomous expansion of reaction networks.

## Source

The reaction network is composed of small gas-phase reactants containing fewer than 4 carbon atoms and 3 oxygen atoms. The design is inspired by Figure 5c from the ChemRxiv paper:

> "Harnessing Machine Learning to Enhance Transition State Search with Interatomic Potentials and Generative Models"

## Directory Structure

- **Reaction Network Labelling:**  
  The molecular species are labelled in the provided PDF file [Reaction Network labelled.pdf](<Reaction Network labelled.pdf>). 
  
- **Molecular Geometries:**  
  The geometry files for the individual molecular species are located in the `molecules/` folder.

- **Reactions:**  
  The initial and final states of the reactions are stored in the `reactions/` folder. The Python script `create_rxns.py` was used to process the molecular geometries and construct these reaction states.

  > **Note: Some of the final state files may require post-processing (e.g., rotating certain chemical species) to resolve the issue of overlapping atoms during the initial interpolation in Minimum Energy Path calculations (such as NEB).**

## Next Steps

The next step is to test the computational framework using these inputs. 

1. Benchmarking Computational Costs:
   - Evaluate the computational expense of generating *ab initio* data for both the minimized reactant configurations and the transition states.
   - Use quantum chemistry packages (e.g., PySCF, Gaussian, etc.) to benchmark various basis sets, density functionals, CPU vs. GPU, and understand the cost and scaling of NEB or dimer calculations for systems of these sizes.

2. Start training MLIP using the data, including TS and surrounding configurations.

3. Build and test the workflow. Test whether we can recover the DFT energy barriers, and whether we can construct or expand the reaction network faster and more efficiently. Explore whether we can discover new reaction pathways. 