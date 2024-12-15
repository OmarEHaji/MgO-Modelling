================================================================================
MgO Defect Energy Analysis Using GULP
================================================================================

This repository contains scripts and input files to study the defect energy in a 
7x7x7 MgO supercell by removing Mg and O atom pairs. The project explores the 
relationship between defect energy and the distance between the removed Mg and 
O atoms.

--------------------------------------------------------------------------------
Project Overview
--------------------------------------------------------------------------------

This project investigates defect formation energies in an MgO supercell using 
GULP (General Utility Lattice Program). The script automates the creation of 
defect input files based on a non-defective supercell input file. These defect 
files are submitted to GULP for computation on a high-performance computing 
(HPC) cluster.

Key objectives:
- Remove Mg and O pairs from a supercell.
- Study defect energy as a function of distance between the removed Mg and O 
  atoms.
- Automate the creation of defect files for efficient workflow.

--------------------------------------------------------------------------------
Features
--------------------------------------------------------------------------------

- Automated Input File Generation: Creates defect files by iteratively removing 
  Mg and O pairs based on their distance.
- Distance-based Filtering: Ensures only unique distances are used for defect 
  generation.
- Customizable Parameters: Allows for adjustments to the supercell size and 
  maximum distance.
- The program can be altered to search for/ calculate interstitial sites 

--------------------------------------------------------------------------------
Workflow
--------------------------------------------------------------------------------

1. Parse the coordinates from the non-defective GULP input file.
2. Identify the most central Mg atom in the supercell.
3. Calculate distances to all O atoms.
4. Create defect files for unique Mg-O distances, adhering to a maximum cutoff 
   range which can be chosen by the user and also avoiding symmetrically
   equivalent defects to reduce uneccessary compute.
5. Submit the defect files to GULP for defect energy calculations.
6. The use of 1 HPC core is sufficient for this calculation as the supercell
   has already been optimised.

--------------------------------------------------------------------------------
Requirements
--------------------------------------------------------------------------------

- Python 3.x
- NumPy
- GULP -- General Utility Lattice Program 
- HPC cluster (for running GULP jobs)
--------------------------------------------------------------------------------

