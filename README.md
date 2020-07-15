# 2D Structure Designer
This Python code allows to generate 2D disordered ensembles of particles with tailored structure properties. 

This algorithm consists of two main parts:
• First, hard (non-overlapping) particles were added using a random sequential approach until the desired filling fraction was reached;
• Second, the difference between the targeted S(q) and the one of the structure was minimised. In detail, the positions of the particles were gradually changed following a gradient descending minimisation protocol.

This code has been used to understand the role of structural correlations in scattering optimisation and in the creation of isotropic structural colourations (see CITATION.txt).

cylinderEnsemble.py and powertools.py are supporting files containing ausiliary functions.
Disordered structures are generated by running MWE.py.

The input parameters are:
r: mean size (in µm)
rsigma: size polidispersity 
ff : filling fraction
bbox: size of the system (in µm)
k0: average distance between the particles (in k-space)
sigmaK: uncertainty in the k (k-error in the structure factor)
sigmaPhi: uncertainty in the phi (angular error in the structure factor)

Different target structure factors are already builtin, namely: Hexagonal, Rectangular, Elliptical, Layered

IMPORTANT: You need to optimise the cSep parameter every time you change the input parameters to ensure that the particles are packing accordingly to the desired structure factor whilst avoiding overlapping.
