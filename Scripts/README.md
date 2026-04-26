# Scripts

The following scripts were used in the preparation and analysis of data for the
study.

  - make_sparse_graph.py - Remove redundant edges from the FEP graph for the DHS variant with coupled H310/H435 titration.
  - merge_graph_avg.py - Combine multiple FEP graphs by averaging per-leg dG values for edges present in more than one graph. 
  - merge_graphs_by_leg.py - Combine leg-specific dG values from different graphs into a single graph. This was used to merge holo Fc-FcRn complex leg with apo Fc solvent leg dG values.
  - prepare_expt_data.py - Merge the experimental data from different sources.
  - rebase_variants.py - Adjust variant names for YTE variant to include the M252Y/S254T/T256E mutations.