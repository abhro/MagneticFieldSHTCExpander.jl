The magnetic field is based on the Altschuler & Newkirk (1969) model, where a potential field is defined
$$\Phi(r, \theta, \phi) = R_\odot \sum_{\ell=0}^\infty \sum_{m=0}^\ell P_\ell^m(\cos \theta) \left[g_{\ell m} \cos(m\phi) + h_{\ell m} \sin(m \phi)\right] \frac{\displaystyle\left(\frac{R_\odot}{r}\right)^{\ell+1} - \left(\frac{R_\odot}{R_{ss}}\right)^{\ell+1} \left(\frac{r}{R_{ss}}\right)^\ell}{\displaystyle \ell + 1 + \ell \left(\frac{R_\odot}{R_{ss}}\right)^{2\ell+1}}$$
and the magnetic field,
$$\mathbf{B}(r, \theta, \phi) = B_r(r, \theta, \phi) \hat{\mathbf{r}} + B_\theta(r, \theta, \phi) \hat{\boldsymbol{\theta}} + B_\phi(r, \theta, \phi) \hat{\boldsymbol{\phi}},$$
is the negative gradient of the potential field:
$$\mathbf{B}(r, \theta, \phi) = -\boldsymbol{\nabla}\Phi(r, \theta, \phi)
  = -\frac{\partial \Phi}{\partial r} \hat{\mathbf{r}} - \frac{1}{r} \frac{\partial \Phi}{\partial \theta} \hat{\boldsymbol{\theta}} - \frac{1}{r \sin \theta} \frac{\partial \Phi}{\partial \phi} \hat{\boldsymbol{\phi}}.$$

Note, here $P_\ell^m$ is the associated Legendre Polynomial under quasi-Schmidt normalization.

For notational convenience, we define the potential function as a product solution, where each function is defined as
```math
\begin{align*}
F_\ell^m(r)      &= \frac{\displaystyle \left(\frac{R_\odot}{r}\right)^{\ell+1} - \left(\frac{R_\odot}{R_{ss}}\right)^{\ell+1} \left(\frac{r}{R_{ss}}\right)^\ell}{\displaystyle \ell + 1 + \ell\left(\frac{R_\odot}{R_{ss}}\right)^{2\ell+1}} \\
G_\ell^m(\theta) &= P_\ell^m(\cos \theta) \\[1ex]
H_\ell^m(\phi)   &= g_{\ell m} \cos(m\phi) + h_{\ell m} \sin(m\phi)
\end{align*}
```

So, the potential function is instead
$$\Phi(r, \theta, \phi) = R_\odot \sum_{\ell=0}^\infty \sum_{m=0}^\ell F_\ell^m(r) G_\ell^m(\theta) H_\ell^m(\phi)$$

### Derivatives
We need the first and second derivatives of each single-variable function for use in the Jacobian of the magnetic field.
```math
\begin{align*}
\frac{d}{dr} F_\ell^m(r)
&= \frac{ \displaystyle
  -(\ell+1) \frac{R_\odot^{\ell+1}}{r^{\ell+2}} - \ell \left(\frac{R_\odot}{R_{ss}}\right)^{\ell+1} \frac{r^{\ell-1}}{R_{ss}^\ell}
}{ \displaystyle
  \ell + 1 + \ell \left(\frac{R_\odot}{R_{ss}}\right)^{2\ell+1}
} \\[1ex]
&= \frac{ \displaystyle
  -\frac{\ell+1}{r} \left(\frac{R_\odot}{r}\right)^{\ell+1} - \frac{\ell}{r} \left(\frac{R_\odot}{R_{ss}}\right)^{\ell+1} \left(\frac{r}{R_{ss}}\right)^\ell
}{ \displaystyle
  \ell + 1 + \ell \left(\frac{R_\odot}{R_{ss}}\right)^{2\ell+1}
}
\\[2ex]
\frac{d^2}{dr^2} F_\ell^m(r)
&= \frac{\displaystyle
  (\ell+1)(\ell+2) \frac{R_\odot^{\ell+1}}{r^{\ell+3}} - \ell(\ell-1) \left(\frac{R_\odot}{R_{ss}}\right)^{\ell+1} \frac{r^{\ell-2}}{R_{ss}^\ell}
}{ \displaystyle
  \ell + 1 + \ell \left(\frac{R_\odot}{R_{ss}}\right)^{2\ell+1}
} \\[1ex]
&= \frac{\displaystyle
  \frac{(\ell+1)(\ell+2)}{r^2} \left(\frac{R_\odot}{r}\right)^{\ell+1} - \frac{\ell(\ell-1)}{r^2} \left(\frac{R_\odot}{R_{ss}}\right)^{\ell+1} \left(\frac{r}{R_{ss}}\right)^\ell
}{ \displaystyle
  \ell + 1 + \ell \left(\frac{R_\odot}{R_{ss}}\right)^{2\ell+1}
}
\end{align*}
```
```math
\begin{align*}
\frac{d}{d\theta} G_\ell^m(\theta) &= - \sin(\theta) ~ \frac{dP_\ell^m(\cos\theta)}{d(\cos\theta)} \\[1ex]
\frac{d^2}{d\theta^2} G_\ell^m(\theta) &= \sin^2(\theta) ~ \frac{d^2 P_\ell^m(\cos\theta)}{d(\cos\theta)^2} - \cos(\theta) ~ \frac{dP_\ell^m(\cos\theta)}{d(\cos\theta)}
\end{align*}
```
```math
\begin{align*}
\frac{d}{d\phi} H_\ell^m(\phi) &= m \left[-g_{\ell m} \sin(m\phi) + h_{\ell m} \cos(m\phi)\right] \\[1ex]
\frac{d^2}{d\phi^2} H_\ell^m(\phi) &= -m^2 [g_{\ell m} \cos(m\phi) + h_{\ell m} \sin(m\phi)] = -m^2 H_\ell^m(\phi)
\end{align*}
```

### Jacobian of magnetic field
The Jacobian for an $\mathbb{R}^3$ is a 3×3 matrix.
For $\mathbf{B}$, it is
```math
\begin{align}
J\mathbf{B}(r, \theta, \phi)
&= \begin{bmatrix}
    \displaystyle -\frac{\partial^2 \Phi}{\partial r^2}   &
    \displaystyle -\frac{1}{r}\frac{\partial^2 \Phi}{\partial \theta \, \partial r} &
    \displaystyle -\frac{1}{r \sin \theta} \frac{\partial^2 \Phi}{\partial \phi \, \partial r} \\
    \displaystyle -\frac{1}{r} \frac{\partial^2 \Phi}{\partial r \, \partial \theta}  &
    \displaystyle -\frac{1}{r^2} \frac{\partial^2 \Phi}{\partial \theta^2} &
    \displaystyle -\frac{1}{r^2 \sin\theta} \frac{\partial^2 \Phi}{\partial\phi \, \partial\theta} \\
    \displaystyle -\frac{1}{r \sin \theta} \frac{\partial^2 \Phi}{\partial r \, \partial \phi} &
    \displaystyle -\frac{1}{r^2 \sin \theta} \frac{\partial^2 \Phi}{\partial\theta \, \partial\phi} &
    \displaystyle -\frac{1}{r^2 \sin^2\theta} \frac{\partial^2 \Phi}{\partial \phi^2}
\end{bmatrix}
\\[2ex]
&= \begin{bmatrix}
    \longleftarrow & \boldsymbol{\nabla}^\mathsf{T} B_r & \longrightarrow \\
    \longleftarrow & \boldsymbol{\nabla}^\mathsf{T} B_\theta & \longrightarrow \\
    \longleftarrow & \boldsymbol{\nabla}^\mathsf{T} B_\phi & \longrightarrow
\end{bmatrix}
\end{align}
```

In terms of the SHTC decomposition, the Jacobian is
```math
J\mathbf{B}(r, \theta, \phi)
= R_\odot \sum_{\ell=0}^\infty \sum_{m = 0}^\ell \begin{bmatrix}
    \displaystyle -\frac{d^2F_\ell^m}{dr^2} G_\ell^m H_\ell^m &
    \displaystyle -\frac{1}{r} \frac{dF_\ell^m}{dr} \frac{dG_\ell^m}{d\theta} H_\ell^m &
    \displaystyle -\frac{1}{r\sin\theta} \frac{dF_\ell^m}{dr} G_\ell^m \frac{dH_\ell^m}{d\phi} &
    \\[1ex]
    \displaystyle -\frac{1}{r} \frac{dF_\ell^m}{dr} \frac{dG_\ell^m}{d\theta} H_\ell^m &
    \displaystyle -\frac{1}{r^2} F_\ell^m \frac{d^2G_\ell^m}{d\theta^2} H_\ell^m &
    \displaystyle -\frac{1}{r^2\sin\theta} F_\ell^m \frac{dG_\ell^m}{d\theta} \frac{dH_\ell^m}{d\phi} &
    \\[1ex]
    \displaystyle -\frac{1}{r\sin\theta} \frac{dF_\ell^m}{dr} G_\ell^m \frac{dH_\ell^m}{d\phi} &
    \displaystyle -\frac{1}{r^2\sin\theta} F_\ell^m \frac{dG_\ell^m}{d\theta} \frac{dH_\ell^m}{d\phi} &
    \displaystyle -\frac{1}{r^2\sin^2\theta} F_\ell^m G_\ell^m \frac{d^2H_\ell^m}{d\phi^2}
\end{bmatrix}
```
Factoring out the parts that appear in $\Phi$, we have
```math
J\mathbf{B}(r, \theta, \phi)
= \sum_{\ell=0}^\infty \sum_{m = 0}^\ell - R_\odot F_\ell^m(r) G_\ell^m(\theta) H_\ell^m(\phi) \begin{bmatrix}
    \displaystyle \frac{1}{F_\ell^m} \frac{d^2F_\ell^m}{dr^2} &
    \displaystyle \frac{1}{r} \frac{1}{F_\ell^m G_\ell^m} \frac{dF_\ell^m}{dr} \frac{dG_\ell^m}{d\theta} &
    \displaystyle \frac{1}{r\sin\theta} \frac{1}{F_\ell^m H_\ell^m} \frac{dF_\ell^m}{dr} \frac{dH_\ell^m}{d\phi}
    \\[1ex]
    \displaystyle \frac{1}{r} \frac{1}{F_\ell^m G_\ell^m} \frac{dF_\ell^m}{dr} \frac{dG_\ell^m}{d\theta} &
    \displaystyle \frac{1}{r^2} \frac{1}{G_\ell^m} \frac{d^2G_\ell^m}{d\theta^2} &
    \displaystyle \frac{1}{r^2\sin\theta} \frac{1}{G_\ell^m H_\ell^m} \frac{dG_\ell^m}{d\theta} \frac{dH_\ell^m}{d\phi}
    \\[1ex]
    \displaystyle \frac{1}{r\sin\theta} \frac{1}{F_\ell^m H_\ell^m} \frac{dF_\ell^m}{dr} \frac{dH_\ell^m}{d\phi} &
    \displaystyle \frac{1}{r^2\sin\theta} \frac{1}{G_\ell^m H_\ell^m} \frac{dG_\ell^m}{d\theta} \frac{dH_\ell^m}{d\phi} &
    \displaystyle \frac{1}{r^2\sin^2\theta} \frac{1}{H_\ell^m} \frac{d^2H_\ell^m}{d\phi^2}
\end{bmatrix}
```
Note that since $\mathbf{B}$ is a gradient field, the Jacobian is symmetric (negative Hessian of $\Phi$).

## Desired quantities
The following are the quantities we want to calculate
* $\displaystyle \Phi = \sum_{\ell=0}^\infty \sum_{m=0}^\ell R_\odot F_\ell^m(r) G_\ell^m(\theta) H_\ell^m(\phi) = \sum_{\ell,m} \Phi_\ell^m$
* $\displaystyle \mathbf{B} = - \sum_{\ell,m} R_\odot F_\ell^m G_\ell^m H_\ell^m \left(\frac{\hat{\mathbf{r}}}{F_\ell^m} \frac{dF_\ell^m}{dr} + \frac{\hat{\boldsymbol{\theta}}}{r G_\ell^m} \frac{dG_\ell^m}{d\theta} + \frac{\hat{\boldsymbol{\phi}}}{r \sin(\theta) H_\ell^m} \frac{dH_\ell^m}{d\varphi}\right)$
* $B = |\mathbf{B}|$
* $\boldsymbol{\nabla} B$. Due to idiosyncrasies of Markdown/LaTeX/MathJax, the mathematical formulation for $\boldsymbol{\nabla} B$ is given outside this list as follows:

```math
\begin{align*}
\boldsymbol{\nabla} B
&= \boldsymbol{\nabla} \sqrt{B_r^2 + B_\theta^2 + B_\phi^2} \\[1ex]
&= \frac{1}{B} \begin{bmatrix}
      \longleftarrow & \dfrac{\partial \mathbf{B}^\mathsf{T}}{\partial r} & \longrightarrow \\
      \longleftarrow & \dfrac{1}{r} \dfrac{\partial \mathbf{B}^\mathsf{T}}{\partial \theta} & \longrightarrow \\
      \longleftarrow & \dfrac{1}{r \sin \theta} \dfrac{\partial \mathbf{B}^\mathsf{T}}{\partial \phi} & \longrightarrow
  \end{bmatrix} \mathbf{B}
\\[1ex]
&= \frac{1}{B} \begin{bmatrix}
      \uparrow & \uparrow & \uparrow \\
      \boldsymbol{\nabla} B_r & \boldsymbol{\nabla} B_\theta & \boldsymbol{\nabla} B_\phi \\
      \downarrow & \downarrow & \downarrow
  \end{bmatrix} \mathbf{B} \\[1ex]
&= \frac{1}{B} (J\mathbf{B})^\mathsf{T} \mathbf{B}
\end{align*}
```
