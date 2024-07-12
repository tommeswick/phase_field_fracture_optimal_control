---
title: 'Space-time formulation, discretization, and computational performance studies for phase-field fracture optimal control problems - reproduction code'
tags:
  - C++
  - Applied Mathematics
  - Space-time
  - Phase-field fracture
  - Optimal control
  - Reduced optimization approach
  - DOpElib
  - deal.II
authors:
  - name: Denis Khimin
    orcid: 0000-0003-1442-3728
    equal-contrib: true
    affiliation: 1
  - name: Thomas Wick
    corresponding author: true
    orcid: 0000-0002-1102-6332
    equal-contrib: true 
    affiliation: 1
  - name: Marc C. Steinbach
    orcid: 0000-0002-6343-9809
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Leibniz Universität Hannover, Institut für Angewandte Mathematik, Welfengarten 1, 30167 Hannover, Germany
   index: 1
 - name: Institution Name, Country
   index: 2
 - name: Independent Researcher, Country
   index: 3
date: 11 July 2024
bibliography: joss_jcp.bib
---

# Summary
The provided codebase is designed to solve phase-field fracture optimal control 
problems, wherein the goal is to achieve a desired fracture in a brittle material
thorugh the application of external forces. 
This algorithmic framework was developed alongside our recent publication
[@JCP2022] and enables the accurate and efficient simulation of phase-field 
optimal control problems in a space-time fashion. 
In this example, the fracture is controlled through Neumann boundary conditions,
and the cost functioncal tracks the distance between.
Our code is based on the open source libraries DOpElib[@dopelib] 
and deal.II [@deal2020].

# Statement of need
Fracture mechanics and optimal control are critical areas in materials science 
and engineering, addressing challenges in predicting and controlling 
crack propagation. The following numerical experiment was developed as 
an example for the DOpElib library[@dopelib], which has already been referenced 
in several scientific publications[@addref]. While existing examples in DOpElib 
solve PDEs with phase-field fracture, they do not include optimization of fractures. 
Therefore, our novel codebase presents a powerful tool for researchers and practitioners.
This implementation not only advances theoretical research based on [@JCP2022], 
but also supports practical applications in fracture mechanics.

# Mathematics
\def\vcentcolon{\mathrel{\mathop\ordinarycolon}}
\providecommand\ordinarycolon{:}
\providecommand*\coloneqq{\vcentcolon\mathrel{\mkern-1.2mu}=}
Let $\Omega \in \mathbb{R}^2$ be a bounded domain and $I = [0,T]$ a time interval with 
an equidistant discretization 
$I = {0} \cup I_1 \cup I_2 \cup\cdots\cup I_M$,
where $I_m = (t_{m-1},t_m]$.
Further, let $C(I,L^2(\Gamma_N))$ with $\Gamma_N \subset \partial \Omega$
be the function space for the Neumann boundary control $q$ and 
$$
X = \{\boldsymbol{v} = (v_u,v_\varphi): \boldsymbol{v} \in L^2(I,V),\,
\partial_t v_\varphi \in L^2(I,H^1(\Omega)^*)\}
$$
with $V \coloneqq H_D^1(\Omega;\mathbb{R}^2) \times H^1(\Omega)$
be the function space for the combined state variable
$\boldsymbol{u} = (u,\varphi)$ consisting of
a vector valued displacement field $u$ and a scalar valued 
phase-field $\varphi$. The phase-field variable acts as a smooth 
indicator function for the fractured area. The time discrete spaces for the 
control and the state variable are defined as
\begin{align*}
Q^0_k &\coloneqq \{r \in C(I,L^2(\Omega)):\, 
r(0) \in Q \text{ and } r|_{I_m} \in \mathbb{P}_0(I_m,Q) \text{ for }m=1,\dots,M \}\\
X^0_k &\coloneqq \{\boldsymbol{v} \in X:\, 
v(0) \in V \text{ and } v|_{I_m} \in \mathbb{P}_0(I_m,V) \text{ for }m=1,\dots,M \}.
\end{align*}
Finally, let $\varphi_d \in L^2(\Omega)$ be a given desired phase-field representing
a desired crack in the domain $\Omega$, let $q \in L^2(\Gamma_N)$ be a nominal control
and $\alpha > 0$ a Tikhonov regularization parameter, then we are seeking a
control $q \in Q_k^0$ and a state $\boldsymbol{u} \in X_k^0$ that solve
\begin{equation}
\label{NLP}
\begin{aligned}
  \min_{q,\boldsymbol{u}} \quad & \mathcal{J}(q,\boldsymbol{u}) \coloneqq
  \frac{1}{2}\sum_{m=1}^M  \left\lVert\varphi(t_m) 
  - \varphi_d\right\rVert_{L^2(\Omega)}^2
  + \frac{\alpha}{2}\sum_{m=1}^M \left\lVert q(t_m) - q_d
  \right\rVert _{L^2(\Gamma_N)}^2 \\
  \text{s.t.} &\text{ $(q,\boldsymbol{u})$ solve \eqref{constraint} for
 $m=1,\dots,M$.}
\end{aligned}
\end{equation}
Herein the constraint \eqref{constraint} is given by the following PDE, where for a given 
control $q(t_m)$, we seek a 
state $\boldsymbol{u}(t_m)$ such that 
\begin{equation}
\label{constraint}
\begin{aligned}
    0
    &= \gamma (\varphi(t_m), \Phi_{\varphi}(t_m))_{\varphi(t_m) > \varphi(t_{m-1})}
    + \eta (\varphi(t_m), \Phi_{\varphi}(t_m)) \\
    &- \gamma (\varphi(t_{m-1}), \Phi_{\varphi}(t_{m}))_{\varphi(t_m) > \varphi(t_{m-1})}
    - \eta (\varphi(t_{m-1}), \Phi_{\varphi}(t_{m})) \\
    &  +\bigg[(g(\varphi)(\mathbb{C} e(u(t_m)), e(\Phi_u)))_{L^2(\Omega)} \\
    &\, + G_c \varepsilon (\nabla \varphi(t_m), \nabla \Phi_{\varphi}(t_m))_{L^2(\Omega)}
    - \frac{G_c}{\varepsilon} (1 - \varphi(t_m), \Phi_{\varphi}(t_m))_{L^2(\Omega)} \\
    &\, + (1 - \kappa) (\varphi(t_m) \mathbb{C} e(u(t_m)) : e(u(t_m)), \Phi_{\varphi}(t_m))\\
    &- (q(t_m), \Phi_{u}(t_m))_{L^2(\Gamma_N)}\bigg](t_m - t_{m-1})
  \end{aligned}
\end{equation}
In the example which is solved by the provided code, the domain is defied as 
a unit square $\Omega = [0,1]^2$ including an initial crack at $[0.5,1] \times \{0.5\}$.
The desired phase-field $\varphi_d$ extends the initial crack to the left, i.e.,
$\varphi_d = 0$ in $[0.3,0.5] \times \{0.5\}$ and $\varphi_d = 0$ everywhere else.
The control is acting in orthogonal direction on the top boundary $[0,1] \times \{1\}$.

Our code solves \eqref{NLP} with Algorithm 1 from [@JCP2022] repeatedly, where in each 
iteration we extend the desired phase-field to the left. Formally in the $k$-th iteration, we 
set $\varphi_d = 0$ in $[0.3 \cdot 0.98^k,0.5] \times \{0.5\}$ and 1 else. In total we perform 
10 updates.

# Installation instructions
1. Install deal.II 9.5.1 via www.dealii.org \
   Download: https://www.dealii.org/download.html \
   Installation instructions: https://www.dealii.org/current/readme.html 
 
2. Install DOpElib via http://www.dopelib.net \
   Download: https://github.com/winnifried/dopelib \
   Installation instructions: https://winnifried.github.io/dopelib/documentation.html 
        in Chapter 2 of the *.pdf manual. 
 
3. Please put the current code into some new folder on your machine. \
   Follow instructions from DOpElib *.pdf manual in Chapter 4 
   (i.e, Section 4.4 Creating new examples)
   to set up all environment variables (for finding deal.II and DOpElib) correctly. \
   Then: build, compile, run as described in Section 4.4 of the dopelib manual.
 
4. The results of this code (see local folder Results/ ) should then reproduce 
   Example 1 (Section 5.1.1) of Khimin et al., JCP, 2022. 

# References
