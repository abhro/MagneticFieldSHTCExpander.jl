# Potential field
The magnetic field is based on the Altschuler & Newkirk (1969) and Zhang, Zhao, and Rassoul (2019) model, where a potential field is defined
```math
\Phi(r, θ, φ)
= R_\odot
  \sum_{ℓ=0}^\infty
    \sum_{m=0}^ℓ
      P_ℓ^m(\cos θ)
      \left[g_{ℓ m} \cos(mφ) + h_{ℓ m} \sin(m φ)\right]
      \frac{
        \displaystyle\left(\frac{R_\odot}{r}\right)^{ℓ+1} - \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \left(\frac{r}{R_{ss}}\right)^ℓ
      }{
        \displaystyle ℓ + 1 + ℓ \left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
      }
```
and the magnetic field,
```math
\mathbf{B}(r, θ, φ) = B_r(r, θ, φ) \hat{\mathbf{r}} + B_θ(r, θ, φ) \hat{\boldsymbol{θ}} + B_φ(r, θ, φ) \hat{\boldsymbol{φ}},
```
is the negative gradient of the potential field:
```math
\mathbf{B}(r, θ, φ)
= - \boldsymbol{∇}\Phi(r, θ, φ)
= - \frac{\partial \Phi}{\partial r} \hat{\mathbf{r}}
  - \frac{1}{r} \frac{\partial \Phi}{\partial θ} \hat{\boldsymbol{θ}}
  - \frac{1}{r \sin θ} \frac{\partial \Phi}{\partial φ} \hat{\boldsymbol{φ}}.
```

Note, here $P_ℓ^m$ is the associated Legendre Polynomial under quasi-Schmidt normalization.

For notational convenience, we define the potential function as a product solution, where each function is defined as
```math
\begin{align*}
F(r) &= \frac{
    \displaystyle \left(\frac{R_\odot}{r}\right)^{ℓ+1} - \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \left(\frac{r}{R_{ss}}\right)^ℓ
}{
    \displaystyle ℓ + 1 + ℓ\left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
} \\
G(θ) &= P_ℓ^m(\cos θ) \\[1ex]
H(φ) &= g_{ℓ m} \cos(mφ) + h_{ℓ m} \sin(mφ)
\end{align*}
```
Note that _F_, _G_, and _H_ all implicitly depend on _ℓ_ and _m_.

So, the potential function is instead
```math
\Phi(r, θ, φ) = R_\odot \sum_{ℓ=0}^\infty \sum_{m=0}^ℓ F(r) G(θ) H(φ)
```

## Derivatives of potential field
We need the first and second derivatives of each single-variable function for use in the Jacobian of the magnetic field.
```math
\begin{align*}
\frac{d}{dr} F(r)
&= \frac{ \displaystyle
  -(ℓ+1) \frac{R_\odot^{ℓ+1}}{r^{ℓ+2}} - ℓ \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \frac{r^{ℓ-1}}{R_{ss}^ℓ}
}{ \displaystyle
  ℓ + 1 + ℓ \left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
} \\[1ex]
&= \frac{ \displaystyle
  -\frac{ℓ+1}{r} \left(\frac{R_\odot}{r}\right)^{ℓ+1} - \frac{ℓ}{r} \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \left(\frac{r}{R_{ss}}\right)^ℓ
}{ \displaystyle
  ℓ + 1 + ℓ \left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
}
\\[2ex]
\frac{d^2}{dr^2} F(r)
&= \frac{\displaystyle
  (ℓ+1)(ℓ+2) \frac{R_\odot^{ℓ+1}}{r^{ℓ+3}} - ℓ(ℓ-1) \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \frac{r^{ℓ-2}}{R_{ss}^ℓ}
}{ \displaystyle
  ℓ + 1 + ℓ \left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
} \\[1ex]
&= \frac{\displaystyle
  \frac{(ℓ+1)(ℓ+2)}{r^2} \left(\frac{R_\odot}{r}\right)^{ℓ+1} - \frac{ℓ(ℓ-1)}{r^2} \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \left(\frac{r}{R_{ss}}\right)^ℓ
}{ \displaystyle
  ℓ + 1 + ℓ \left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
}
\end{align*}
```
```math
\begin{align*}
\frac{d}{dθ} G(θ) &= - \sin(θ) ~ \frac{dP_ℓ^m(\cosθ)}{d(\cosθ)} \\[1ex]
\frac{d^2}{dθ^2} G(θ) &= \sin^2(θ) ~ \frac{d^2 P_ℓ^m(\cosθ)}{d(\cosθ)^2} - \cos(θ) ~ \frac{dP_ℓ^m(\cosθ)}{d(\cosθ)}
\end{align*}
```
```math
\begin{align*}
\frac{d}{dφ} H(φ) &= m \left[-g_{ℓ m} \sin(mφ) + h_{ℓ m} \cos(mφ)\right] \\[1ex]
\frac{d^2}{dφ^2} H(φ) &= -m^2 [g_{ℓ m} \cos(mφ) + h_{ℓ m} \sin(mφ)] = -m^2 H_ℓ^m(φ)
\end{align*}
```

### Jacobian of magnetic field
The Jacobian for an $\mathbb{R}^3$ is a 3×3 matrix.
For **B**, it is
```math
\begin{align}
J\mathbf{B}(r, θ, φ)
&= \begin{bmatrix}
    \displaystyle -\frac{\partial^2 \Phi}{\partial r^2}                        &
    \displaystyle -\frac{1}{r}\frac{\partial^2 \Phi}{\partial θ \, \partial r} &
    \displaystyle -\frac{1}{r \sin θ} \frac{\partial^2 \Phi}{\partial φ \, \partial r} \\
    \displaystyle -\frac{1}{r} \frac{\partial^2 \Phi}{\partial r \, \partial θ} &
    \displaystyle -\frac{1}{r^2} \frac{\partial^2 \Phi}{\partial θ^2}           &
    \displaystyle -\frac{1}{r^2 \sinθ} \frac{\partial^2 \Phi}{\partialφ \, \partialθ} \\
    \displaystyle -\frac{1}{r \sin θ} \frac{\partial^2 \Phi}{\partial r \, \partial φ} &
    \displaystyle -\frac{1}{r^2 \sin θ} \frac{\partial^2 \Phi}{\partialθ \, \partialφ} &
    \displaystyle -\frac{1}{r^2 \sin^2θ} \frac{\partial^2 \Phi}{\partial φ^2}
\end{bmatrix}
\\[2ex]
&= \begin{bmatrix}
    \longleftarrow & \boldsymbol{∇}^\mathsf{T} B_r & \longrightarrow \\
    \longleftarrow & \boldsymbol{∇}^\mathsf{T} B_θ & \longrightarrow \\
    \longleftarrow & \boldsymbol{∇}^\mathsf{T} B_φ & \longrightarrow
\end{bmatrix}
\end{align}
```

In terms of the SHTC expansion, the Jacobian is
```math
J\mathbf{B}(r, θ, φ)
= R_\odot \sum_{ℓ=0}^\infty \sum_{m = 0}^ℓ \begin{bmatrix}
    \displaystyle -\frac{d^2F}{dr^2} G H &
    \displaystyle -\frac{1}{r} \frac{dF}{dr} \frac{dG}{dθ} H &
    \displaystyle -\frac{1}{r\sinθ} \frac{dF}{dr} G \frac{dH}{dφ} &
    \\[1ex]
    \displaystyle -\frac{1}{r} \frac{dF}{dr} \frac{dG}{dθ} H &
    \displaystyle -\frac{1}{r^2} F \frac{d^2G}{dθ^2} H &
    \displaystyle -\frac{1}{r^2\sinθ} F \frac{dG}{dθ} \frac{dH}{dφ} &
    \\[1ex]
    \displaystyle -\frac{1}{r\sinθ} \frac{dF}{dr} G \frac{dH}{dφ} &
    \displaystyle -\frac{1}{r^2\sinθ} F \frac{dG}{dθ} \frac{dH}{dφ} &
    \displaystyle -\frac{1}{r^2\sin^2θ} F G \frac{d^2H}{dφ^2}
\end{bmatrix}
```
Factoring out the parts that appear in \Phi, we have
```math
J\mathbf{B}(r, θ, φ)
= \sum_{ℓ=0}^\infty \sum_{m = 0}^ℓ - R_\odot F(r) G(θ) H(φ) \begin{bmatrix}
    \displaystyle \frac{1}{F} \frac{d^2F}{dr^2} &
    \displaystyle \frac{1}{r} \frac{1}{F G} \frac{dF}{dr} \frac{dG}{dθ} &
    \displaystyle \frac{1}{r\sinθ} \frac{1}{F H} \frac{dF}{dr} \frac{dH}{dφ}
    \\[1ex]
    \displaystyle \frac{1}{r} \frac{1}{F G} \frac{dF}{dr} \frac{dG}{dθ} &
    \displaystyle \frac{1}{r^2} \frac{1}{G} \frac{d^2G}{dθ^2} &
    \displaystyle \frac{1}{r^2\sinθ} \frac{1}{G H} \frac{dG}{dθ} \frac{dH}{dφ}
    \\[1ex]
    \displaystyle \frac{1}{r\sinθ} \frac{1}{F H} \frac{dF}{dr} \frac{dH}{dφ} &
    \displaystyle \frac{1}{r^2\sinθ} \frac{1}{G H} \frac{dG}{dθ} \frac{dH}{dφ} &
    \displaystyle \frac{1}{r^2\sin^2θ} \frac{1}{H} \frac{d^2H}{dφ^2}
\end{bmatrix}
```
Note that since **B** is a gradient field, the Jacobian is symmetric (negative Hessian of Φ).

## Desired quantities
The following are the quantities we want to calculate

- $\displaystyle \Phi = \sum_{ℓ=0}^\infty \sum_{m=0}^ℓ R_\odot F(r) G(θ) H(φ) = \sum_{ℓ,m} \Phi_ℓ^m$, where $\Phi_ℓ^m = R_\odot F(r) G(θ) H(φ)$

- $\displaystyle \mathbf{B} = - \sum_{ℓ,m} R_\odot F(r) G(θ) H(φ) \left(\frac{\hat{\mathbf{r}}}{F} \frac{dF}{dr} + \frac{\hat{\boldsymbol{θ}}}{r G} \frac{dG}{dθ} + \frac{\hat{\boldsymbol{φ}}}{r \sin(θ) H} \frac{dH}{d\varphi}\right)$

- _B_ = |**B**|

- **∇**_B_. Due to idiosyncrasies of Markdown/LaTeX/MathJax, the mathematical formulation for **∇**_B_ is given outside this list as follows:

```math
\begin{align*}
\boldsymbol{∇} B
&= \boldsymbol{∇} \sqrt{B_r^2 + B_θ^2 + B_φ^2} \\[1ex]
&= \frac{1}{B} \begin{bmatrix}
      \longleftarrow & \dfrac{\partial \mathbf{B}^\mathsf{T}}{\partial r} & \longrightarrow \\
      \longleftarrow & \dfrac{1}{r} \dfrac{\partial \mathbf{B}^\mathsf{T}}{\partial θ} & \longrightarrow \\
      \longleftarrow & \dfrac{1}{r \sin θ} \dfrac{\partial \mathbf{B}^\mathsf{T}}{\partial φ} & \longrightarrow
  \end{bmatrix} \mathbf{B}
\\[1ex]
&= \frac{1}{B} \begin{bmatrix}
      \uparrow & \uparrow & \uparrow \\
      \boldsymbol{∇} B_r & \boldsymbol{∇} B_θ & \boldsymbol{∇} B_φ \\
      \downarrow & \downarrow & \downarrow
  \end{bmatrix} \mathbf{B} \\[1ex]
&= \frac{1}{B} (J\mathbf{B})^\mathsf{T} \mathbf{B}
\end{align*}
```

# References
* Altschuler, M.D., Newkirk, G. Magnetic fields and the structure of the solar corona. _Sol Phys_ **9**, 131–149 (1969). https://doi.org/10.1007/BF00145734
* Zhang, Ming, Lulu Zhao, and Hamid K. Rassoul. Stochastic Propagation of Solar Energetic Particles in Coronal and Interplanetary Magnetic Fields. Journal of Physics: Conference Series 1225, no. 1 (May 1, 2019): 012010. https://doi.org/10.1088/1742-6596/1225/1/012010.

