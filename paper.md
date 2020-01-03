---
title: 'dispersioneering: A Matlab package for the dispersion relation of magnetized plasmas'
tags:
  - Matlab
  - plasma physics
  - dispersion relation
  - magnetized
  - hot plasma
  - cold plasma
  - kinetic effects
authors:
  - name: David L. Green
    orcid: 0000-0003-3107-1170
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Oak Ridge National Laboratory
   index: 1
date: 02 January 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

It is often of interest to plasma physicists to know which plasma waves can propagate in a magnetized plasma at a given frequency. The calculation of such reduces to solving for which $k_{\perp}$ (or $k_{\parallel}$) the 0-D time-harmonic wave equation operator has zero determinant. For a cold (i.e., no finite-temperature effects) magnetized plasma, simple analytic polynomial expressions for this problem are available in several texts [@stix1992waves; @brambilla1998kinetic; @swanson2003plasma] and $k_{\perp}(\omega)$ (or $k_{\parallel}(\omega)$) being solveable via the quadratic equation. However, for a hot plasma where finite-temperature effects must be included, the problem becomes one of non-linear root finding with a rather complex objective function. As such, dispersion solvers for hot magnetized plasmas are not readily available. This repository aims to make available this capability for Matlab.  

``dispersioneering`` is a Matlab package for calculating the general dispersion relation for hot or cold magnetized plasmas. It allows 1-D specification of plasma profiles of magnetic field $B(x)$ (T), density $n(x)$ ($m^{-3}$), and temperature $T(x)$ (eV) for an aribtrary number of species at a specified frequency. It calculates both hot and cold solutions, with the harmonic number an input parameter for the hot calculation. 

Figures can be included like this: ![Test figure.](figures/output1.png)

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures



# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
