---
title: 'DRAFT: VoronoiFVM.jl: Finite Volume Solver for Coupled Nonlinear Partial Differential Equations'
tags:
  - Julia
  - Partial differential equations
  - Finite Volume method
  - Automatic differentiation
authors:
  - name: Jürgen Fuhrmann
    orcid: 0000-0003-4432-2434
    affiliation: 1
affiliations:
  - name: Weierstrass Insitute for Applied Analysis and Stochastics, Berlin, Germany
    index: 1
date: January 31, 2023
bibliography: paper.bib
---
# Summary

The Voronoi box based finite volume method is a special case of the two-point flux finite volume method
which used Voronoi cells as control volumes. The method is nowadays theoretically well understood,
and it has a wide range of possible applications, calling for flexible, efficient and user friendly implementations.
The Julia programming language provides efficient means to implement algorithms in the realms of scientific computing
and data science. It combines the possibility of user friendly problem descriptions with the ability to generate
fast code via  just-ahead-of-time compilation. Its interface oriented API design
growing package ecosystem provides ways to combine different algorithms, like e.g. automatic differentiation,
method of lines for time integration.
In this paper, we introduce VoronoiFVM.jl, a Julia package for solving nonlinear systems of
convection-diffusion-reaction equations  with the help of the Voronoi box based finite volume method.


# The method


Regard $n$ coupled PDEs in $d$ dimensional domain $\Omega\subset \mathbb R^d$ with boundary $\Gamma=\partial\Omega$ and time interval $(0,T)$:
$$
\begin{aligned}
  \partial_t \mathbf s(\mathbf u) + \nabla \cdot \vec{\mathbf j}(\mathbf u, \vec \nabla \mathbf u) + \mathbf r(\mathbf u) &= \mathbf f(\vec x,t)\quad\text{in}\; \Omega\\
  -\vec{\mathbf j}(\mathbf u, \vec \nabla \mathbf u)\cdot \vec n + \mathbf b(\mathbf u,\vec x,t) &=0 \quad \text{on}\; \Gamma
\end{aligned}
$$

  Here, $\mathbf u(\vec x,t)=(u_1(\vec x,t)\dots u_n(\vec x,t)): \Omega\times (0,T) \to \mathbb R^n$ represents the unknwons of the problem to be solved for. The storage term $\mathbf s: \mathbb R^n \to \mathbb R^n$ denotes the local amount of species. The reaction term $\mathbf r: \mathbb R^n \to \mathbb R^n$ gives the reaction between species  within $\Omega$. 
The vector function $\vec{\mathbf j}: \mathbb R^n \times \mathbb R^{nd} \to \mathbb R^{nd}$ represents the species flux. The function $\mathbf f: \Omega\times (0,T) \to \mathbb R^{n}$ provides the sources and sinks within the domain.
  
  The boundary term $\mathbf b: \mathbb R^n\times\partial\Omega\times (0,T)$ describes nonlinear Robin boundary conditions which  via corresponding choices of $\mathbf b$ can result in Neumann boundary conditions. Dirichlet boundary conditions will  be realized by the penalty method.

For simplicity, we omitted possible additional dependencies of $\mathbf s, \vec{\mathbf j}, \mathbf r$ on $\vec x$ and $t$.



|                                                  | $\mathbf s(u)$  | $\vec{\mathbf j}(u,\vec\nabla u)$        | $\mathbf r(u)$  |
| -----------------------------------------------: | :-------------: | :--------------------------------------: | :-------------: |
| Diffusion (Fick's law):                          | $u$             | $-D\vec \nabla u$                       | 0               |
| Convection-diffusion in velocity field $\vec v$: | $u$             | $-D\vec \nabla u + u \vec v$             | 0               |
| Unsaturated porous media flow:                   | $\theta(u)$     | $-k(\theta(u)) (\vec \nabla u -\vec g)$ | 0               |
| Charge equilibrium in electrolytes:              | $0$             | $-\varepsilon\vec \nabla u$              | $\sinh(u)$      |
|--------------------------------------------------|-----------------|------------------------------------------|-----------------|

Examples for constitutive relationships for scalar problems

  
  
  Table \ref{tab:scalar} gives a number of examples for scalar problems.
  Notable systems of equations which can be described in this context include multiphase flow in porous media,
  Reaction-diffusion systems and charge transport in semiconductors and electrolytes.
  
  
  
  
  Assume $\Omega\subset \mathbb{R}^d$ is a polygonal domain such that
  $\partial\Omega=\bigcup_{m\in\mathcal{G}} \Gamma_m$, where $\Gamma_m$ are planar such that $\vec{n}|_{\Gamma_m}=\vec{n}_m$.
  We subdivide $\Omega$ into into a finite number of  control volumes  such that $\bar\Omega= \bigcup_{k\in \mathcal N} \bar \omega_k$  and 
  -  $\omega_k$ are open  convex domains such that  $\omega_k\cap \omega_l = \emptyset$ if
    $\omega_k\neq \omega_l$ 
  -  The intersections of the closures of the control volumes  $\sigma_{kl}=\bar \omega_k \cap \bar \omega_l$ are either empty,    points or straight lines.
    If $|\sigma_{kl}|>0$ we say that $\omega_k$, $\omega_l\in$ are neighbours.
  -   $\vec{n}_{kl}\perp \sigma_{kl}$: normal of $\partial\omega_k$ at $\sigma_{kl}$
  -   $\mathcal{N}_k = \{l \in \mathcal{N}: |\sigma_{kl}|>0\}$: set of neighbours of $\omega_k$
  -   $\gamma_{k}=\partial\omega_k \cap \Gamma$: boundary part of $\partial\omega_k$


  
  $\Rightarrow$ $\partial\omega_k= \left(\bigcup_{l\in \mathcal{N}_k} \sigma_{kl}\right)  \bigcup\gamma_{k}$

  To each control volume $\omega_k$ assign a collocation
  point: $\vec{x}_k \in \bar \omega_k$ which
  for a given function $u:\Omega \to \mathbb{R}$ will allow to associate its value $u_k^m=u(\vec{x}_k,t^m)$
  as the value of the discrete unknown  at $\vec{x}_k$

  -  _Admissibility condition:_ if $l\in \mathcal N_k$ then the   line $\vec{x}_k\vec{x}_l$ is orthogonal to $\sigma_{kl}$..
      For two neigboring control volumes $\omega_k, \omega_l$ , this will allow to approximate
      $\vec\nabla u \cdot \vec{n}_{kl} \approx  \frac{u_l - u_k}{|\vec x_k - \vec x_l}$
  -   _Placement of boundary unknowns at the boundary:_ if
      $\omega_k$ is situated at the boundary, i.e. for 
      $|\partial \omega_k \cap \partial \Omega| >0$,
      then $\vec{x}_k \in \partial\Omega$. This condition allows to apply boundary condition in  a direct way.


In certain cases (rectangular meshes, 1D problems), such subdivision can be constructed in an explicit manner.
In general, the procedure to obtain this  construction starts with the construction of a boundary conforming Delaunay triangulation
of $\Omega$ with vertices $\vec{x}_k$ \cite{Si}. Such triangulations can be obtained by grid generators \cite{Triangle} and \cite{TetGen}.
With the triangulation given, $\omega_k$ can be constructed as restricted Voronoi cells $\omega_k$.
In 2D, corners of these restricted Voronoi cells are  either triangle  circumcenters or midpoints of boundary   edges. 
This approach fulfills the admissibility condition in a natural way and uses  the simplex nodes $\vec x_k$ as collocation points.
As a consequence, the degrees of freedom are the same as the ones of the P1 finite element method.
    

  Integrating over a space-time control volume gives
  $$
  \begin{aligned}
    \int_{t^{m-1}}^{t^{m}}\int_{\omega} \mathbf s(\mathbf u(\vec x,t))\,d\vec x,dt  +
    \int_{t^{m-1}}^{t^{m}}\left( \int_{\partial\omega} \vec{\mathbf j} \cdot\vec n\,ds
    +  \int_{\omega} \mathbf r(\mathbf u)\, d\vec x -  \int_{\omega} \mathbf f \, d\vec x\right)\,dt = 0.\\
\end{aligned}
$$
  Using Newton-Leibniz rule,  Gauss' theorem and replacing integrals over $\omega_k$ by one-point quadrature centered at $x_k$
  results in 
  
$$
    |\omega_k|\frac{\mathbf s(\mathbf u_k^{m})-\mathbf s(\mathbf u_k^{m-1})}{t^{m}-t^{m-1}}+ |\omega_k|\mathbf r(\mathbf u_k^m) - |\omega_k|\mathbf f_k 
+    \sum_{\omega_l  \text{neigbour of} \omega_k} \int_{\omega_k\cap\omega_l} \vec{\mathbf j} \cdot \vec n \; +  \int_{\gamma_k} \vec{\mathbf j} \cdot \vec n \; =0
$$

In next step, we replace integrals over boundries by another one-point quadrature. For the flux between neigboring control volumes,
we introduce the flux function $\mathbf g(u_k,u_l)$ which approximates $\vec{\mathbf j}\cdot(\vec x_k - \vec x_l)$. Furthermore, we
use the boundary condition.
$$
    |\omega_k|\frac{\mathbf s(\mathbf u_k^{m})-\mathbf s(\mathbf u_k^{m-1})}{t^{m}-t^{m-1}}+ |\omega_k|\mathbf r(\mathbf u_k^m) - |\omega_k|\mathbf f_k 
    +    \sum_{\omega_l  \text{neigbour of} \omega_k} \frac{|\omega_k\cap\omega_l|}{|\vec x_k-\vec x_l|} \mathbf  g(\mathbf u^m_k , \mathbf u^m_l) +
    |\gamma_k| \mathbf b(u_k^m) =0
$$

This is a system of nonlinear equations for $u_k^m$ realising the implicit Euler method.
In the scalar case, for the heat equation with Robin boundary conditions

In the scalar linear case
$$
\begin{aligned}
  u_t - \nabla D\cdot \nabla u  = 0 \; \text{in} \; \Omega \quad  \partial_{\vec n} u + \alpha u - \beta=  0 \; \text{on} \; \partial \Omega,
\end{aligned}
$$
one uses for approximation  $s(u)=u, r(u)=0, g(u_k, u_l) = D(u_k-u_l), f=0, b(u)=\alpha u - \beta$. For $D,\alpha >0$, this leads
to a linear system for $u_k^m$ with an M-matrix and correspomding results for positivity of the solution, mass conseervation.
Similar results can be obtained  for many other well defined situations. \cite{EGH}.


  -  Inherent local mass conservation
  -  Robust concepts for upwinding and stabilization $\Rightarrow$ avoid unphysical oscillations, negative concentrations 
  -  Provable consistency to thermodynamic principles - discrete 2nd law of thermodynamics 
  -  Convergence theory for nonlinear systems based on compactness methods


 Gajewski/Gärtner, ZAMM, doi 10.1002/zamm.19960760502, 1996\\
    F.,Computer Phys. Comm. doi 10.1016/j.cpc.2015.06.004, 2015\\
    Gaudeul, F., Numer. Math, doi 10.1007/s00211-022-01279-y 2022



With $\alpha = \frac1{\varepsilon}$ and $\beta=\frac1{\varepsilon}u_=0$, for small $\varepsilon$, one realizes the Dirichlet
boundary condition  via the penalty method. In particular, Dirichlet values are almost exactly attained at the corresponding boundary
points $x_k$. With proper implementation and choice of $\varepsilon$, they are attained up to floating point accuracy.

# Implementation

Core of the implementation is the need to solve the nonlinear system $A(u^m)= f^m$ described by \eqref{eq:impeul}. 
For the solution of the nonlinear system of equations, we use Newton's method, possibly with a simple damping strategy and
coupled to parameter embedding. This requires the calculation  of the residual $R(u)=A(u)-f$ and the Jacobian $A'(u)$
in each iteration step. Here, we provide some details of the implementation.


## API Idea
The API for the prolem description is inspired by the structure of the discretization:

-  Create grid (e.g. 2D grid from coordinate vectors $X$,$Y$)\\
  \verb|grid=simplexgrid(X,Y)|
-  Write functions describing physics of the problem
   -   Storage $\mathbf s: \mathbb R^n \to \mathbb R^n$: ``s!(f,u,node)= ...``
   -   Reaction $\mathbf r: \mathbb R^n \to \mathbb R^n$: ``r!(f,u,node)= ...``
   -   Flux between neighbor REVs $\mathbf g: \mathbb R^n \times \mathbb R^n \to \mathbb R^n$:   ``g!(f,u,edge)= ...``
   -   Boundary reaction $\mathbf b: \mathbb R^n \to \mathbb R^n$:     ``b!(f,u,node)= ...``
-  Create discrete system: ``system=VoronoiFVM.System(grid,flux=g!,storage=s!,reaction=r!, breaction=b!)``
-  Solve via built-in implicit Euler method  ``solve(system;times=[0,T])``


## General assembly strategy
For the assembly of $R(U)$ and $A'(u)$ it is not necessary to construct the Voronoi cells in an
explicit way.  It is sufficient to obtain know the connectivity of between the $\vec x_k$
given by the primary triangulation used to construct the Voronoi diagram, and the geometry data $|\omega_k|$,
$\frac{|\omega_k\cap\omega_l|}{|\vec x_k-\vec x_l|}$ and $|\gamma_k|$. These can be split into contributions
from each simplex of the triangulation, and thus the assembly loop is run over the triangles instead of
the Voronoi cells. This approach resembles the standard assembly approach for P1 finite elemenet methods.
It allows to utilize the ability of grid generators to describe boundaries between subdomains which align
with the simplex boundaries. As a result, both  $R(U)$ and $A'(u)$ are assembled from local contributions
from simplices, and the local contributions to  $A'(u)$ are obtained from partial derivatives of te constituting
functions with respect to the local unknowns. 


## Automatic differentiation to calculate the Jacobian

Forward mode automatic differentiation allows to evaluate a nonlinear function such that both its value and its derivative are obtained at once. A straightforward implementation is based on dual numbers $\mathbb D$ defined by extending the set of real numbers $\mathbb R$. Similar to introducing the imaginary unit $i$ with $i^2=-1$ to define the complex numbers, one  introduces a special number $\varepsilon$ to define the set of dual numbers as $\mathbb D= \{ a + b\varepsilon \; |\; a,b \in \mathbb R\}$. With assuming $\varepsilon^2=0$, the 
evaluation of a polynomial $p(x)=\sum_{i=0}^n p_i x^i$ on a dual number $a+\varepsilon$ yields
$$
\begin{aligned}
      p(a+\varepsilon) = \sum_{i=0}^n p_i a^i + \sum_{i=1}^n i p_i  a^{i-1} \varepsilon = p(a)+p'(a)\varepsilon.
\end{aligned}
$$
This fact can be generalized to differentiable functions of several variables and to multivariate dual numbers, allowing for the calculation of partial derivatives.  The Julia computer language via the package \cite{ForwardDiff.jl} provides an easily accessible implementation of dual number arithmetic helping to evaluate nonlinear functions along with their derivatives. As a consequence,  it is indeed not necessary to write
extra code for  obtaining partial derivatives, and the API is restricted to the functions shown in section \ref{sec:API}


\subsection{Linear system solution}
The linear systems of equations are described by sparse matrices. From 
the Julia ecosystem we use the package LinearSolve.jl \cite{LinearSolve.jl}  which allows to
dispatch to  various direct and preconditioned iterative solvers.



# Further features and generalization

## Multiphysics, bulk-interface coupling and additional terms
The package always assumes a subdivision of the domain and its boundary into regions: $\Omega = \cup_{r\in \mathcal R_\Omega \Omega_r}$
and $\Gamma = \cup_{r\in \mathcal R_\Gamma \Gamma_r}$. Therefore it is possible to use different constutive functions in
different regions. Moreover, it is possible to define different species sets in different subregions.

This includes the possibilty to define species whose support lies on the boundary of the domain.  Additional terms in the problem
description are able to describbe boundary  storage, flux and sources.

## Method of lines via DifferentialEquations.jl
VoronoiFVM.jl is accompanied by the package VoronoiFVMDiffEq.jl which allows to reformulate the finite volume system
\eqref{eq:impeul} as an \verb|ODEProblem| in the sense of that package. This allows to use any of impressive
list of highly optimized stiff ODE and DAE solvers availablable from that package. 


# References
