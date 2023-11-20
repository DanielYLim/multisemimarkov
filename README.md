---
author: Yongho Lim
date: "`r format(Sys.time(), '%d %B, %Y')`"
header-includes:
- "`\\usepackage{tikz}`{=tex}"
- "`\\usepackage{pgfplots}`{=tex}"
- "`\\usepackage{pgf,tikz,pgfplots}`{=tex}"
- "`\\pgfplotsset{compat=1.15}`{=tex}"
- "`\\usepackage{mathrsfs}`{=tex}"
- "`\\usetikzlibrary{arrows}`{=tex}"
- "`\\usetikzlibrary{patterns}`{=tex}"
- "`\\usetikzlibrary{positioning}`{=tex}"
output:
  pdf_document: default
title: multisemimarkov package
---

`\usepackage{tikz}`{=tex}

`\usepackage{pgfplots}`{=tex}

`\usepackage{pgf,tikz,pgfplots}`{=tex}

`\pgfplotsset{compat=1.15}`{=tex}

`\usepackage{mathrsfs}`{=tex}

`\usetikzlibrary{arrows}`{=tex}

`\usetikzlibrary{patterns}`{=tex}

`\usetikzlibrary{positioning}`{=tex}

```{=html}
<!-- README.md is generated from README.Rmd. Please edit that file -->
```
`{r, include = FALSE} knitr::opts_chunk$set(   collapse = TRUE,   comment = "#>",   fig.path = "man/figures/README-",   out.width = "100%" )`

# multisemimarkov

```{=html}
<!-- badges: start -->
```
```{=html}
<!-- badges: end -->
```
Provides functionality with Multistate semi-Markov model with masked
causes

## Installation

You can install the development version of `multisemimarkov` like so:

``` r
install.packages("multisemimarkov")
```

## Example

This is a basic example which shows you how to solve a common problem:

`{r example} library(multisemimarkov) ## basic example code`

Multi-state model with partially latent cured and not cured states

```{=tex}
\begin{center}
\begin{tikzpicture}
\node[draw, rectangle, rounded corners, text centered] (treatment) {Treatment (1)};
\node[draw, rectangle, above right= 1.5 of treatment, rounded corners, text centered] (notcured) {Not Cured};
\node[draw, rectangle, right=1.5 of notcured, rounded corners, text centered] (recurrence) {Recurrence (2)};
\node[draw, rectangle, below=4 of recurrence, rounded corners, text centered] (cancerdeath) {Cancer Death (3)};
\node[draw, rectangle,  below=4.1 of notcured, rounded corners, text centered] (cured) {Cured};
\draw[->, line width= 1, dashed] (treatment) --  node[above,font=\scriptsize]{}(notcured);
\draw[->, line width= 1, dashed] (treatment) --  node[above,font=\scriptsize]{}(cured);
\draw [->, line width= 1] (notcured) -- node[above,font=\scriptsize]{}(recurrence);
\draw [->, line width= 1] (recurrence) -- node[above,font=\scriptsize]{}(cancerdeath);
\draw [->, line width= 1] (notcured) -- node[above,font=\scriptsize]{}(cancerdeath);
\end{tikzpicture}
\end{center}
```
Multi-State Model Structure with Death Due to Other Causes State

```{=tex}
\begin{center}
\begin{tikzpicture}
% nodes %
\node[draw, rectangle, rounded corners, text centered] (treatment) {Treatment (1)};
\node[draw, rectangle, above right= 1.5 of treatment, rounded corners, text centered] (notcured) {Not Cured};
\node[draw, rectangle, right=1.5 of notcured, rounded corners, text centered] (recurrence) {Recurrence (2)};
\node[draw, rectangle, below=2 of recurrence, rounded corners, text centered] (cancerdeath) {Cancer Death (3)};
\node[draw, rectangle,  below=5 of notcured, rounded corners, text centered] (cured) {Cured};
\node[draw, rectangle, right=0.6 of cured, rounded corners, text centered] (otherdeath) {Death due to Other Causes (4)};
% edges %
\draw[->, line width= 1, dashed] (treatment) --  node[above,font=\scriptsize]{}(notcured);
\draw[->, line width= 1, dashed] (treatment) --  node[above,font=\scriptsize]{}(cured);
\draw [->, line width= 1] (notcured) -- node[above,font=\scriptsize]{}(recurrence);
\draw [->, line width= 1] (recurrence) -- node[above,font=\scriptsize]{}(cancerdeath);
\draw [->, line width= 1] (notcured) -- node[above,font=\scriptsize]{}(cancerdeath);
\draw [->, line width= 1, dashed] (notcured) -- node[above,font=\scriptsize]{}(otherdeath);
\draw [->, line width= 1, dashed] (cured) -- node[above,font=\scriptsize]{}(otherdeath);
\end{tikzpicture}
\end{center}
```
# Data Generation Algorithm for Unmasked Data

Random samples of
$\{(t_{1i}, t_{2i}, \delta_{12i}, \delta_{13i}, \delta_{23i}), i = 1, 2, \dots, n\}$
can be obtained from the following data generation algorithm.

```{=tex}
\begin{enumerate}
    \item Generate cure status for each individual $i$ from Bernoulli distribution with cure probability $1-p$. We set $p = 0.70$.
    \item If individual $i$ is uncured, generate $T_{1i}$ from $F_{10}(t_{1i}) = \Pr(T_{1i} \leq t_{1i} | \text{Not Cured}) = 1-S_{10}(t_{1i})=1-\exp \left[ -\int_{0}^{t_{1i}} \left( \lambda_{12}(u)+\lambda_{13}(u) \right) du \right]$. We assume log-logistic distribution to model the cause specific hazard functions for disease progression events $M = k$, $k = 2, 3$ conditional on not being cured:  
\begin{equation} \label{csh12_13}
\lambda_{1k} (t_{1}) = \frac{\left(\frac{\beta_{1k}}{\alpha_{1k}} \right) \left( \frac{t_{1}}{\alpha_{1k}} \right)^{\beta_{1k}-1}}{1 + \left(  \frac{t_{1}}{\alpha_{1k}} \right)^{\beta_{1k}}}, \qquad k = 2, 3.
\end{equation}
We set $\alpha_{12} = 2.0$, $\beta_{12} = 4.0$ and $\alpha_{13} = 3.5$, $\beta_{13} = 3.0$ to have more observed cancer recurrences than cancer deaths without experiencing recurrence. 
    \item For an uncured individual $i$, set $M_i = 2$ with probability $ \frac{\lambda_{12}(T_{1i}) }{\lambda_{12}(T_{1i}) +\lambda_{13}(T_{1i}) }$, and set $M_i = 3$ otherwise.  
    \item For an uncured individual $i$, if $M_i = 2$, generate $T_{2i}$ from $C(F_{1|2}(T_{1i}), F_{2}(T_{2i}))$ in (\ref{Eq2}), the conditional joint distribution of $T_{1}$ and $T_{2}$  for patients who experience cancer recurrence. We consider the Clayton copula function
\begin{equation}\label{Clayton}
C_{\phi}(u_1, u_2)=(u_1^{-\phi}+u_2^{-\phi}-1)^{-1/\phi}, \qquad  \phi>0,
\end{equation}
where $\phi$ is the dependence parameter set to $\phi = 0.857$. It gives the Kendall's $\tau = \phi/(\phi+2) = 0.3$, a moderate dependence between the sequential gap times. The marginal distribution of $T_2$ for subjects who have experienced recurrence was considered as the Weibull distribution with $F_2(t_2)=1-\exp\left[ -( \frac{t_2}{\alpha_{23}})^{\beta_{23}} \right]$. We set $\alpha_{23} = 2.5$ and $\beta_{23} =1.5$. 
    \item Generate censoring times $\{C_i, i = 1, 2, \dots, n \}$ from Uniform $(0, 10)$. 
    \item If individual $i$ is cured, generate time-to-death due to other causes from $F_{14|\text{Cured}}(t_{1}) = \Pr(T_{1i} \leq t_{1i}, M_{i}=4 | \text{Cured})$ and set the censoring time $C_i$ as the time to death due to other causes if the censoring time obtained in the previous step is greater than the time to death due to other causes. (Time to death due to other causes is considered as a censoring time when the estimation is performed using \eqref{S2Eq12}.) We considered
    \begin{equation}\label{F14Cured}
F_{14 | \text{Cured}}(t_1)=\Pr( T_1 \leq t_1, M = 4 | \text{Cured})=1-\exp\left[ -\left( \frac{t_{1}}{\alpha_{14 | \text{Cured}}}\right)^{\beta_{14 | \text{Cured}}} \right].
\end{equation}
We set $\alpha_{14 | \text{Cured}} = 7.0$ and $\beta_{14 | \text{Cured}} = 4.0$.    
    \item Obtain $t_{1i} = \min(T_{1i}, C_i)$, $\delta_{12i} = I[T_{1i} = t_{1i}, M_i = 2]$, $\delta_{13i} = I[T_{1i} = t_{1i},  M_i = 3]$.
    If $M_i = 2$, obtain $t_{2i} = \min(T_{2i}, C_i-t_{1i})$ and $\delta_{23i} = I[T_{2i} = t_{2i}, M_i = 2]$. 
\end{enumerate}
```
