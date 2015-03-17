This document is formatted in 
[Python Markdown](https://pypi.python.org/pypi/Markdown).

"Worms" are a shorthand name for "Poisson wavelet multi-scale edges".

The fundamental worming algorithm is described fully in [Hornby et al., (1999)](http://dx.doi.org/10.1046/j.1365-246x.1999.00788.x).

Here, we present a brief summary of the highlights.

First and foremost, this is a _processing_ algorithm (in the spirit of Jon Claerbout's "processing vs. inversion" dichotomy). The inversion is _induced_ by a physical interpretation of the algorithm [Hornby et al. (2002)](http://www.cosis.net/abstracts/EGS02/05568/EGS02-A-05568.pdf) recognized after the initial theory paper.

Briefly, horizontal gradients of the Green's function of the Poisson equation -- when upward continued and amplitude corrected appropriately -- form both the analysing and synthesizing functions of a continuous (2D) vector-valued wavelet transform of any function derived from Poisson's equation ("potential fields"). In geophysics, gravity and magnetic fields (the latter suitably operated on in the usual fashion found in textbooks) are prominent examples of such fields. 

Multi-scale edges are the locations and scaled amplitude of local maxima in horizontal gradient of the field. When stitched together horizontally, they form a "worm", and when different levels of upward continuation are connected (i.e. different wavelet scales), such "worms" form "worm sheets".

So-far, worms form a viable (and physically/mathematically relevant) "skeletonization" of the source-free field (which extends everywhere above ground). [Mallat and Zhong, (1992)](http://dx.doi.org/10.1109/34.142909) show that the original signal (in our case, the the field before upward continuation) can be reconstructed to a close approximation from the information contained in the wavelet multi-scale edges alone. This means that the worms themselves contain enough information to re-build the original field (again within a close approximation).

Up to this point, we have a skeletonization, and a method for re-constructing the field from that skeletonization. (This is still "processing" in a Claerbout-ian sense.) This produces a nice, interpretable visualization in and of itself -- and it was used as such for the first 5+ years of worm deployment in the Australian mining industry. However, [Hornby et al, (2002)](http://www.cosis.net/abstracts/EGS02/05568/EGS02-A-05568.pdf) showed -- via a physical interpretation of the inverse wavelet transform derived in the earlier paper -- that the vector-valued wavelet in fact can alternatively be viewed as the field due to source-dipoles of appropriate strength and orientation. From this perspective, the maxima in horizontal gradients (i.e. the worms) are now the locations of the strongest (or highest concentration) of dipole sources oriented perpendicular to the worms. This leads to a natural physical interpretation of the worms being lateral boundaries in (approximately piecewise-constant) source regions. In other words, the positive side of the dipole is in the direction of the higher-valued (monopole) source region, while the negative side of the dipole is similarly in the direction of the lower-valued region. We have found an "edge" in the rock sources.

Moreover, again drawing on the result from Mallat and Zhong (1992) we see that these edges (and their associated monopole source regions) generate the observed field (not upward continued) in exactly the same fashion as before. The edges are in fact the same worm sheets as before, just draped underground now and re-interpreted as lateral source boundaries.

This leads to the conclusion that out of the infinite families of source distributions that can exactly generate the observed potential field, the underground worm sheets "induce" a source inversion which locates (some) lateral boundaries. 

Experience has shown that the while we lose a little field fidelity to the residual via the "close approximation" mentioned above, we gain a lot of physical intuition and interpretability from the worm visualisation. That being said, the inversion is clearly incorrect in detail, since it gets smoother with depth, unlike the real world.

A picture is worth 1000 words, so the following figure tries to illustrate many of the concepts described above:

![Worm Cartoon](./WormCartoon.png "Worm Cartoon")

