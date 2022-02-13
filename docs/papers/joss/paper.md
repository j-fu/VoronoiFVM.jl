---
title: 'VoronoiFVM.jl: PDE System solver using the Voronoi finite volume method.'
tags:
  - Julia
  - Partial differential equations
  - Finite Volume method
authors:
  - name: Jürgen Fuhrmann^[Custom footnotes for e.g. denoting who the corresponding author is can be included like this.]
    orcid: 0000-0003-0872-7098
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    affiliation: 2
  - name: Author with no affiliation
    affiliation: 3
affiliations:
 - name: Weierstrass Institute for Applied Analysis and Stochastics
   index: 1
date: February 28, 2022
bibliography: paper.bib
---

# Summary

# Statement of need
Systems of partial differential equations are used to model coupled physical processes evolving in space and time. Except for some special cases, solutions to such equations cannot be described by simple mathematical formulas. Instead, one mostly relies on approximating them by large, but finite dimensional systems of equations which can be solved on a computer.

The finite volume method provides an approximation strategy especilly well suited for partial differential equations derived from conservation laws for physical quantities. It essentially consists in the subdivision of a computational domain into a finite number of representative elementary volumes, called control volumes, and the accounting for interacting  sources/sinks of quantities in each control volume, and the transfer of quantities betweeen neighboring
control volumes.


# The method

- General form of systems
- Finite volume picture

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Implementation

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References


