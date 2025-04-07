# Potential field
The magnetic field is based on the [AN1969](@citet) and [ZZR2019](@citet) model, where a potential field is defined
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
= - \grad\Phi(r, θ, φ)
= - \pdv{\Phi}{r} \hat{\mathbf{r}}
  - \frac{1}{r} \pdv{\Phi}{θ} \hat{\boldsymbol{θ}}
  - \frac{1}{r \sin θ} \pdv{\Phi}{φ} \hat{\boldsymbol{φ}}.
```

Note, here $P_ℓ^m$ is the associated Legendre Polynomial under quasi-Schmidt normalization.

For notational convenience, we define the potential function as a product solution, where each function is defined as
```math
\begin{align*}
F_ℓ(r) &= \frac{
    \displaystyle \left(\frac{R_\odot}{r}\right)^{ℓ+1} - \left(\frac{R_\odot}{R_{ss}}\right)^{ℓ+1} \left(\frac{r}{R_{ss}}\right)^ℓ
}{
    \displaystyle ℓ + 1 + ℓ\left(\frac{R_\odot}{R_{ss}}\right)^{2ℓ+1}
} \\
G_{ℓm}(θ) &= P_ℓ^m(\cos θ) \\[1ex]
H_{ℓm}(φ) &= g_{ℓm} \cos(mφ) + h_{ℓ m} \sin(mφ)
\end{align*}
```
Note that _F_, _G_, and _H_ all implicitly depend on _ℓ_ and _m_ (only _ℓ_ in _F_'s case).

So, the potential function is instead
```math
\Phi(r, θ, φ) = R_\odot \sum_{ℓ=0}^\infty F_ℓ(r) \sum_{m=0}^ℓ G_{ℓm}(θ) H_{ℓm}(φ)
```

## Derivatives of potential field
We need the first and second derivatives of each single-variable function for use in the Jacobian of the magnetic field.
```math
\begin{align*}
\dv{}{r} F_ℓ(r)
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
\frac{d^2}{dr^2} F_ℓ(r)
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
\frac{d}{dθ} G_{ℓm}(θ) &= - \sin(θ) ~ \frac{dP_ℓ^m(\cosθ)}{d(\cosθ)} \\[1ex]
\frac{d^2}{dθ^2} G_{ℓm}(θ) &= \sin^2(θ) ~ \frac{d^2 P_ℓ^m(\cosθ)}{d(\cosθ)^2} - \cos(θ) ~ \frac{dP_ℓ^m(\cosθ)}{d(\cosθ)}
\end{align*}
```
```math
\begin{align*}
\frac{d}{dφ} H_{ℓm}(φ) &= m \left[-g_{ℓ m} \sin(mφ) + h_{ℓ m} \cos(mφ)\right] \\[1ex]
\frac{d^2}{dφ^2} H_{ℓm}(φ) &= -m^2 [g_{ℓ m} \cos(mφ) + h_{ℓ m} \sin(mφ)] = -m^2 H_ℓ^m(φ)
\end{align*}
```

### Jacobian of magnetic field
The Jacobian for an $\mathbb{R}^3$ is a 3×3 matrix.
For **B**, it is
```math
\begin{align*}
J\mathbf{B}(r, θ, φ)
&= \begin{bmatrix}
    \displaystyle -\frac{∂^2 \Phi}{∂ r^2}                        &
    \displaystyle -\frac{1}{r}\frac{∂^2 \Phi}{∂ θ \, ∂ r} &
    \displaystyle -\frac{1}{r \sin θ} \frac{∂^2 \Phi}{∂ φ \, ∂ r} \\[0.5ex]
    \displaystyle -\frac{1}{r} \frac{∂^2 \Phi}{∂ r \, ∂ θ} &
    \displaystyle -\frac{1}{r^2} \frac{∂^2 \Phi}{∂ θ^2}           &
    \displaystyle -\frac{1}{r^2 \sinθ} \frac{∂^2 \Phi}{∂φ \, ∂θ}  \\[0.5ex]
    \displaystyle -\frac{1}{r \sin θ} \frac{∂^2 \Phi}{∂ r \, ∂ φ} &
    \displaystyle -\frac{1}{r^2 \sin θ} \frac{∂^2 \Phi}{∂θ \, ∂φ} &
    \displaystyle -\frac{1}{r^2 \sin^2θ} \frac{∂^2 \Phi}{∂ φ^2}
\end{bmatrix}
\\[2ex]
&= \begin{bmatrix}
    \longleftarrow & \grad^\TT B_r & \longrightarrow \\[0.5ex]
    \longleftarrow & \grad^\TT B_θ & \longrightarrow \\[0.5ex]
    \longleftarrow & \grad^\TT B_φ & \longrightarrow
\end{bmatrix}
\end{align*}
```

In terms of the SHTC expansion, the Jacobian is
```math
J\mathbf{B}(r, θ, φ)
= R_\odot \sum_{ℓ=0}^\infty \sum_{m = 0}^ℓ \begin{bmatrix}
    \displaystyle -\frac{d^2F_ℓ}{dr^2} G_{ℓm} H_{ℓm} &
    \displaystyle -\frac{1}{r} \dv{F_ℓ}{r} \dv{G_{ℓm}}{θ} H_{ℓm} &
    \displaystyle -\frac{1}{r\sinθ} \dv{F_ℓ}{r} G_{ℓm} \dv{H_{ℓm}}{φ} &
    \\[2ex]
    \displaystyle -\frac{1}{r} \dv{F_ℓ}{r} \dv{G_{ℓm}}{θ} H_{ℓm} &
    \displaystyle -\frac{1}{r^2} F_ℓ \frac{d^2G_{ℓm}}{dθ^2} H_{ℓm} &
    \displaystyle -\frac{1}{r^2\sinθ} F_ℓ \dv{G_{ℓm}}{θ} \dv{H_{ℓm}}{φ} &
    \\[2ex]
    \displaystyle -\frac{1}{r\sinθ} \dv{F_ℓ}{r} G_{ℓm} \dv{H_{ℓm}}{φ} &
    \displaystyle -\frac{1}{r^2\sinθ} F_ℓ \dv{G_{ℓm}}{θ} \dv{H_{ℓm}}{φ} &
    \displaystyle -\frac{1}{r^2\sin^2θ} F_ℓ G_{ℓm} \frac{d^2H_{ℓm}}{dφ^2}
\end{bmatrix}
```
Factoring out the parts that appear in Φ, we have
```math
J\mathbf{B}(r, θ, φ)
= \sum_{ℓ=0}^\infty \sum_{m = 0}^ℓ - R_\odot F_ℓ(r) G_{ℓm}(θ) H_{ℓm}(φ) \begin{bmatrix}
    \displaystyle \frac{1}{F_ℓ} \frac{d^2F_ℓ}{dr^2} &
    \displaystyle \frac{1}{r} \frac{1}{F_ℓ G_{ℓm}} \dv{F_ℓ}{r} \dv{G_{ℓm}}{θ} &
    \displaystyle \frac{1}{r\sinθ} \frac{1}{F_ℓ H_{ℓm}} \dv{F_ℓ}{r} \dv{H_{ℓm}}{φ}
    \\[2ex]
    \displaystyle \frac{1}{r} \frac{1}{F_ℓ G_{ℓm}} \dv{F_ℓ}{r} \dv{G_{ℓm}}{θ} &
    \displaystyle \frac{1}{r^2} \frac{1}{G_{ℓm}} \frac{d^2G_{ℓm}}{dθ^2} &
    \displaystyle \frac{1}{r^2\sinθ} \frac{1}{G_{ℓm} H_{ℓm}} \dv{G_{ℓm}}{θ} \dv{H_{ℓm}}{φ}
    \\[2ex]
    \displaystyle \frac{1}{r  \sinθ} \frac{1}{F_ℓ H_{ℓm}} \dv{F_ℓ}{r} \dv{H_{ℓm}}{φ} &
    \displaystyle \frac{1}{r^2\sinθ} \frac{1}{G_{ℓm} H_{ℓm}} \dv{G_{ℓm}}{θ} \dv{H_{ℓm}}{φ} &
    \displaystyle \frac{1}{r^2\sin^2θ} \frac{1}{H_{ℓm}} \frac{d^2H_{ℓm}}{dφ^2}
\end{bmatrix}
```
Note that since **B** is a gradient field, the Jacobian is symmetric (negative Hessian of Φ).

### Gradient of the field strength

**∇**_B_ = **∇**|**B**| can be determined from the Jacobian of **B**. What follows is a proof using index notation:
```math
\pdv{B}{x_j}
= \pdv{}{x_j} \left(B_i B_i\right)^{1/2}
= \frac{1}{2} \left(B_i B_i\right)^{-1/2} \left(2 B_i \pdv{B_i}{x_j}\right)
= \frac{B_i \pdv{B_i}{x_j}}{\left(B_k B_k\right)^{1/2}}
= \frac{B_i}{B} \pdv{B_i}{x_j}
```
which implies ``∇B = \hat{\mathbf{b}} ⋅ J\mathbf{B}``

## Desired quantities
The following are the quantities we want to calculate

- The potential field

  ```math
  Φ = \sum_{ℓ=0}^\infty \sum_{m=0}^ℓ R_\odot F_ℓ(r) G_{ℓm}(θ) H_{ℓm}(φ) = \sum_{ℓ,m} \Phi_ℓ^m,
  ```

  where $\Phi_ℓ^m = R_\odot F_ℓ(r) G_{ℓm}(θ) H_{ℓm}(φ)$

- The magnetic field

  ```math
  \mathbf{B} = - \sum_{ℓ,m} R_\odot F_ℓ(r) G_{ℓm}(θ) H_{ℓm}(φ) \left(\frac{\hat{\mathbf{r}}}{F_ℓ} \dv{F_ℓ}{r} + \frac{\hat{\boldsymbol{θ}}}{r G_{ℓm}} \dv{G_{ℓm}}{θ} + \frac{\hat{\boldsymbol{φ}}}{r \sin(θ) H_{ℓm}} \frac{dH_{ℓm}}{d\varphi}\right)
  ```

- The magnetic field strength, _B_ = |**B**|

- The gradient of the magnetic field strength, **∇**_B_

  ```math
  \begin{align*}
  \grad B
  &= \grad \sqrt{B_r^2 + B_θ^2 + B_φ^2} \\[1ex]
  &= \frac{1}{B} \begin{bmatrix}
        \longleftarrow & \dpdv{\mathbf{B}^\TT}{r}                     & \longrightarrow \\
        \longleftarrow & \dfrac{1}{r} \dpdv{\mathbf{B}^\TT}{θ}        & \longrightarrow \\
        \longleftarrow & \dfrac{1}{r \sin θ} \dpdv{\mathbf{B}^\TT}{φ} & \longrightarrow
    \end{bmatrix} \mathbf{B}
  \\[1ex]
  &= \frac{1}{B} \begin{bmatrix}
        \uparrow   & \uparrow   & \uparrow  \\
        \grad B_r  & \grad B_θ  & \grad B_φ \\
        \downarrow & \downarrow & \downarrow
    \end{bmatrix} \mathbf{B} \\[1ex]
  &= \frac{1}{B} (J\mathbf{B})^\TT \mathbf{B}
  \end{align*}
  ```

## References
```@bibliography
Pages = ["physical-theory.md"]
```
