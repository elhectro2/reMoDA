# -*- coding: utf-8 -*-

"""
Modules for checking if a feature of a trajectory of a protein, related to unfolding, has been currently analyzed.
If the trajectory has not been analyzed, the analysis is performed. If it has been analyzed, it is properly saved.

Each module has two main functions: One for performing the analysis and the other for reading the results.

Modules contained in this package are the following ones:

- energy: Calculation of the total energy of the protein during the simulation
- gyration: Calculation of the radius of gyration (compactness) of the protein during the simulation
- h_bonds_global: Calculation of the hydrogen bonds established by the protein during the simulation
- h_bonds_local: Calculation of the hydrogen bonds established by and with the residue of interest
- native_contacts_global: Calculation of the native contacts established by the protein during the simulation
- native_contacts_local: Calculation of the native contacts established by and with the residue of interest
- pressure: Calculation of the pressure of the system during the simulation
- rmsd: Calculation of the Root Mean Square Deviation compared to the minimized structure
- rmsf: Calculation of the Root Mean Square Fluctuation compared to the minimized structure, by residue
- sasa: Calculation of the Solvent Accessible Surface Area of the whole protein, residue by residue and in total
- ss: Calculation of the number of residues on a secondary structure of interest (helices for alpha proteins,
      beta sheets or strands for beta proteins) and of unstructured regions (coils)

Some modules require also some helper bash functions:
-TODO LIST HERE
"""

__version__ = "1.0.0"

