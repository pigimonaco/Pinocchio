# Pinocchio

PINOCCHIO is a fast code to generate catalogues of cosmological dark matter halos with known mass, position, velocity and merger history.

It is able to reproduce, with very good accuracy, the hierarchical formation of dark matter halos from a realization of an initial (linear) density perturbation field, given on a 3D grid.

Its setup is similar to that of a conventional N-body simulation, but it is based on the powerful Lagrangian Perturbation Theory. It runs in just a small fraction of the computing time taken by an equivalent N-body simulation, producing promptly the merging histories of all halos in the catalog.

PINOCCHIO is distributed under a gnu-gpl v2 license.

The latest version, V4.1, is available in github. This version is written in C, presents much better scaling properties, and is highly recommended. Please read the DOCUMENTATION file for further details.

Important: version 4.1 contains an important bug fix with respect to v 4.0

Latest version: 4.1.3, yet another bug fix in the lightcone reconstruction, and improvements in the writing of snapshots.

Webpage:
http://adlibitum.oats.inaf.it/monaco/pinocchio.html

Papers:

1. The PINOCCHIO algorithm: pinpointing orbit-crossing collapsed hierarchical objects in a linear density field, Pierluigi Monaco, Tom Theuns & Giuliano Taffoni, 2002, MNRAS, 331, 587, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2002.MNRAS.331.587.pdf)
2. Predicting the Number, Spatial Distribution and Merging History of Dark Matter Haloes, Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Fabio Governato, Tom Quinn & Joachim Stadel, 2002, ApJ, 564, 8, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2002.ApJ.564.8.pdf)
3. PINOCCHIO and the Hierarchical Build-Up of Dark-Matter Haloes, Giuliano Taffoni, Pierluigi Monaco & Tom Theuns, 2002, MNRAS, 333, 623, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/taffoni.2002.MNRAS.333.623.pdf)
4. An accurate tool for the fast generation of dark matter halo catalogs, Pierluigi Monaco, Emiliano Sefusatti, Stefano Borgani, Martin Crocce, Pablo Fosalba, Ravi Sheth & Tom Theuns, 2013, MNRAS, 433, 2389, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2013.MNRAS.433.2389.pdf)
5. nIFTy cosmology: Galaxy/halo mock catalogue comparison project on clustering statistics, C.-H. Chuang, C. Zhao, F. Prada et al., 2015, MNRAS, 452, 686, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/chuang.2015.MNRAS.452.686.pdf)
6. Improving the prediction of dark matter halo clustering with higher orders of Lagrangian Perturbation Theory, E. Munari, P. Monaco, E. Sefusatti, E. Castorina, F.G. Mohammad, S. Anselmi & S. Borgani, submitted to MNRAS. arXiv:1605.04788, [PDF](http://arxiv.org/abs/1605.04788)
7. Simulating cosmologies beyond LambdaCDM with Pinocchio, L.A. Rizzo, F. Villaescusa-Navarro, P. Monaco, E. Munari, S. Borgani, E. Castorina & E. Sefusatti, 2016, submitted to JCAP. arXiv:1610.07624, [PDF](https://arxiv.org/abs/1610.07624)
8. Approximated methods for the generation of dark matter halo catalogs in the age of precision cosmology, P.Monaco, 2016, Galaxies, N. 4, 53, special issue: "Dark Matter: Large versus Small Scale Structures", ed. J. Gaite & A. Diaferio. [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2016.Galaxies.4.53.pdf)
