# Input structures

This folder contains all prepared input structure files used in the manuscript.

## Primary Fc-FcRn complex all-atom models used for retrospective FEP calculation

Holo *wt* and YTE Fc-FcRn complex all-atom models used for retrospective FEP
calculations of Fc and FcRn mutations, prepared from 6WNA with protonation
states assigned using OPLS5-based pKas and pH 6.0.

- 6WNA_IgG1_WT_opls5_pH6.mae = for Fc and FcRn mutations in *wt* background
- 6WNA_IgG1_YTE_opls5_pH6.mae = for Fc mutations in YTE background

## Apo-conformation Fc all-atom models used for solvent leg-only FEP simulations

Apo-conformation *wt* and YTE Fc all-atom models used for solvent leg-only FEP
simulations of Fc mutations, prepared from the M252R variant low-pH structure
(9D06) with protonation states matching the holo models.

- 9D06_WT_opls5_ph6.mae = for apo solvent leg Fc mutations in *wt* background
- 9D06_YTE_opls5_ph6.mae = for apo solvent leg Fc mutations in YTE background

## L309D/Q311H Fc model used for DHS variant calculations

Holo Fc-FcRn complex all-atom model of the L309D/Q311H variant used for FEP
calculations investigating the coupling between DHS mutations and H310/H435,
prepared from 6WNA with protonation states as in the *wt* holo model, and
mutated residues H:D309 and H:H311 modeled as ASP and HIP, respectively.

- 6WNA_IgG1_DH_opls5_pH6.mae

## T307W models used for calculations reported in Supplementary Table 6

Holo- and apo-conformation all-atom models of the T307W Fc variants used for
FEP calculations reported in Supplementary Table 6, prepared from 6WNA and
9D06, respectively, with varying rotamers of the T307W side chain. Used for
investigation of T307W outlier when using apo-conformation model for FEP
solvent legs.

- holo_T307W_HIP310_HIP435_alt_rot1_min_6A.mae
- holo_T307W_HIP310_HIP435_alt_rot2_min_6A.mae
- holo_T307W_HIP310_HIP435_alt_rot3_min_6A.mae
- holo_T307W_HIP310_HIP435_alt_stab_frame.mae
- apo_T307W_HIP310_HIP435_alt_rot1_min_6A.mae
- apo_T307W_HIP310_HIP435_alt_rot2_min_6A.mae
- apo_T307W_HIP310_HIP435_alt_rot3_min_6A.mae
- apo_T307W_HIP310_HIP435_alt_stab_frame.mae
