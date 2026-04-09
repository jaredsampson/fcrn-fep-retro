Low frequency normal modes, describing collective atomic motions typically associated with the most
flexible regions, were calculated for all monomeric holo and apo Fc structures as well as for a dimeric YTE
Fc structure (prepared from PDB 4N0U via crystallographic symmetry) using a computationally efficient
normal mode analysis approach (A. Hoffmann, S. Grudinin, NOLB: Nonlinear Rigid Block Normal-Mode Analysis Method, Journal of Chemical Theory and
Computation 13 (5) (2017) 2123–2134. doi:10.1021/acs.jctc.7b00197), which was run using the command:

***./NOLB structure.pdb***

To visualize, open a reference structure and a structure corresponding to a given low frequency mode in pymol, e.g.:  
***pymol 4N0U_Fc_dimer.pdb 4N0U_Fc_dimer_nlb_1.pdb***

and hit "play button" in the right bottom corner to view the normal mode motion.

In all structures (both monomers and the dimer), the lowest-frequency mode corresponds to a hinge-like
“bending” motion, analogous to the symmetric in-plane bending vibrational mode of a water molecule. The
second-lowest-frequency mode represented an out-of-plane “twisting” motion. Linear combinations of these
two lowest-frequency modes could describe transitions between all observed structures (see subfolder "transitions") 
