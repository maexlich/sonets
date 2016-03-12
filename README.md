# sonets
Generates Second Order Networks (SONETs) with prescribed second order motif frequencies

## Overview

C++ code that generates random directed unweighted graphs with prescribed second order statistics.
Additional structure (i.e., higher-order statistics) are approximately minimized using a near maximum entropy approach.

Current version of code can generate a single statistically homogeneous network or a two population network
where the statistics can vary as a function of the populations participating in each connection.

Second order statistics are prescribed based on frequencies of four types of two-edge motifs (subnetworks containing two edges):

* reciprocal motifs (two nodes reciprocally connected:  o <-> o)
* convergent motifs (two nodes each connecting onto a third node:  o -> o <- o)
* divergent motifs (two nodes each receving a connection from a third node: o <- o -> o)
* chain motifs (three nodes connected in sequence: o -> o -> o)

Reference: L. Zhao, B. Beverlin II, T Netoff, and D. Q. Nykamp.  Synchronization from second order network connectivity statistics.
*Frontiers in Computational Neuroscience*, 5:28, 2011.  [Publisher's web site](http://dx.doi.org/10.1007/s10827-011-0373-5)

## Requirements

* C++ compiler (such a [gcc](https://gcc.gnu.org/)
* [GNU Software Library (GSL)](http://www.gnu.org/software/gsl/) 
(If installing a gsl package on a system, such as on Ubuntu, where the development library is separated from the main library,
also install the development library.)
