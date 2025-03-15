# Pinocchio V5.1

PINOCCHIO is a fast code to generate catalogues of cosmological dark matter halos with known mass, position, velocity and merger history. As a byproduct, it can also generate the matter density field and, upon post-processing, its lensing potential.

PINOCCHIO starts from the realisation of a linear density field on a regular grid, as in the generation of the initial conditions of a cosmological N-body simulation. It is based on Lagrangian Perturbation Theory (LPT), excursion set theory and ellipsoidal collapse. The code is made of two main parts, the computation of collapse times and LPT displacements for each particle, and the grouping of collapsed particles into haloes (`fragmentation'), with the costruction of halo merger histories and light-cone with continuous time sampling. 
Collapse times are computed by Gaussian-smoothing the linear density field on many smoothing radii, then computing the second derivatives of the initial potential with FFTs; these are used to compute the collapse redshift of each particle using ellipsoidal collapse. We define the inverse collapse time as $F=1+z_{\rm c}$, and store its highest value $F_{\rm max}$ for all smoothing radii. At the final smoothing radius $R=0$ (meaning that the variance of the linear density field is only limited by the Lagrangian grid), the LPT displacement fields are computed, amounting to four vectors for each particle (of the three 3LPT displacement fields we only compute the first two, the third rotational term being negligible). 
The second part of the code uses the collapse times and displacements to group particles into haloes, with an algorithm that mimics hierarchical clustering. Because collapse here is identified with orbit crossing, collapsed particles are not necessarily contained in haloes, they may be part of the filamentary network that joins haloes. So particles may be classified into uncollapsed (still in single-stream regime), filaments and halo particles. Having recognised the haloes without really running the simulation (we just performed a single 3LPT time-step when needed), we can see PINOCCHIO as a halo finder that works on the Lagrangian space of initial conditions, plus a 3LPT engine to place the haloes at the right position.

PINOCCHIO is distributed under a gnu-gpl v2 license.

[Documentation]https://github.com/pigimonaco/Pinocchio/blob/master/DOCUMENTATION
[Installation]https://github.com/pigimonaco/Pinocchio/blob/master/INSTALLATION

Github: https://github.com/pigimonaco/Pinocchio

Webpage: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

Papers:

1. The PINOCCHIO algorithm: pinpointing orbit-crossing collapsed hierarchical objects in a linear density field, Pierluigi Monaco, Tom Theuns & Giuliano Taffoni, 2002, MNRAS, 331, 587, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2002.MNRAS.331.587.pdf)
2. Predicting the Number, Spatial Distribution and Merging History of Dark Matter Haloes, Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Fabio Governato, Tom Quinn & Joachim Stadel, 2002, ApJ, 564, 8, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2002.ApJ.564.8.pdf)
3. PINOCCHIO and the Hierarchical Build-Up of Dark-Matter Haloes, Giuliano Taffoni, Pierluigi Monaco & Tom Theuns, 2002, MNRAS, 333, 623, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/taffoni.2002.MNRAS.333.623.pdf)
4. An accurate tool for the fast generation of dark matter halo catalogs, Pierluigi Monaco, Emiliano Sefusatti, Stefano Borgani, Martin Crocce, Pablo Fosalba, Ravi Sheth & Tom Theuns, 2013, MNRAS, 433, 2389, [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2013.MNRAS.433.2389.pdf)
5. Improving the prediction of dark matter halo clustering with higher orders of Lagrangian Perturbation Theory, E. Munari, P. Monaco, E. Sefusatti, E. Castorina, F.G. Mohammad, S. Anselmi & S., 2017, MNRAS, 465,4658 [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/munari.2017.MNRAS.465.4658.pdf)
6. Simulating cosmologies beyond LambdaCDM with Pinocchio, L.A. Rizzo, F. Villaescusa-Navarro, P. Monaco, E. Munari, S. Borgani, E. Castorina & E. Sefusatti, 2017, JCAP, 01/2017, 008 [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/rizzo.2017.JCAP.0117.008.pdf)
7. Approximated methods for the generation of dark matter halo catalogs in the age of precision cosmology, P.Monaco, 2016, Galaxies, N. 4, 53, special issue: "Dark Matter: Large versus Small Scale Structures", ed. J. Gaite & A. Diaferio. [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/monaco.2016.Galaxies.4.53.pdf)
8. Fast numerical method to generate halo catalogues in modified gravity (part I): second-order Lagrangian perturbation theory, C. Moretti, S. Mozzon, P. Monaco, E. Munari, M. Baldi, M. 2020, MNRAS, 493, 1153 [PDF](http://adlibitum.oats.inaf.it/monaco/Papers/moretti.2020.mnras.493.1153.pdf)
9. Euclid Preparation. Simulating thousands of Euclid spectroscopic skies, Euclid Consortium: P. Monaco, G. Parimbelli, Y. Elkhashab et al., in preparation.
