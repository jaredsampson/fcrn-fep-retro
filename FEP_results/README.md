# FEP_results

This directory contains the results of the free energy perturbation (FEP) calculations performed in this study. Original .fmp files from which the .csv files here are too large for Github but may be requested from the authors.

The .csv files were produced by running `protein_fep_groups.py` from the Schrödinger 2024-4 release on the original .fmp files with `-ph 6.0`.

The results are divided into the following subdirectories:

- `pka-calc`: OPLS4 and OPLS5 pKa calculations for titratable interface residues on Fc and FcRn.
- `original-retrospective`: Original retrospective FEP calculations performed in this study using the full dataset.
- `histidine-coupling`: Subset of mutations selected for histidine coupling analysis to Fc H310 and H435 pH sensors.
- `apo-solvent`: Results after running solvent legs with apo Fc models (wt and YTE) and merging those solvent legs with the original retrospective complex and fragment legs (using `merge_graphs_by_leg.py` from this repo).