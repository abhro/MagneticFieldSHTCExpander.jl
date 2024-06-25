module MagneticFieldSHTCExpander

export magneticfield, collectmagneticfield

using OffsetArrays: OffsetArray, OffsetMatrix
using LinearAlgebra: LowerTriangular, norm
using LegendrePolynomials: Plm
using ForwardDiff: derivative as derivativefd

using Memoization: @memoize
using ThreadSafeDicts: ThreadSafeDict

const R_0 = 1 # sun's radius in solar radius
const R_SS = 2.5 # source surface in solar radius (PFSS)

const QSNORM = Val(:schmidtquasi)

# declare the return type as its own struct thing because it makes
# it easier during development. once everything has been solidified,
# this should be deleted and the returns made concrete.
const BField = @NamedTuple begin
    Φ::Float64                     # magnetic potential at (r, θ, φ)
    B::Vector{Float64}             # mag field at (r, θ, φ)
    #normB::Float64                # norm of mag field at (r, θ, φ)
    jacobianB::Matrix{Float64}     # Jacobian matrix of mag field at (r, θ, φ)
    ∇normB::Vector{Float64}        # gradient of norm of mag field
end

function magneticfield(
        r, θ, φ,            # position in spherical coordinates (r is in solar radius)
        g::AbstractMatrix,  # indices = (0:lmax, 0:lmax)
        h::AbstractMatrix   # indices = (0:lmax, 0:lmax)
    )::BField


    if axes(g) != axes(h)
        throw(DimensionMismatch("g and h must be the same size and index"))
    end

    # ensure radius within solar surface and source surface
    if r < R_0 || r > R_SS
        throw(DomainError(r, "r must be within [R_0, R_ss] = [$R_0, $R_SS]"))
    end

    if θ < 0 || θ > π
        throw(DomainError(θ, "θ must be within [0, π]"))
    end

    if φ < 0 || φ > 2π
        throw(DomainError(φ, "φ must be within [0, 2π]"))
    end

    # set up array of allowed ℓ values
    ℓ_axes = axes(g)[1]

    # PHYSICAL THEORY (see Physical-theory.md)
    # Each ℓ,m term of the potential can be written as a product of three
    # functions depending uniquely on r, θ, φ each.

    cosθ = cos(θ)
    sinθ = sin(θ)
    if sinθ ≈ 0
        @warn("θ is near 0 or π, causing sinθ to be almost 0. This may cause NaN issues")
    end

    # create vectors and matrices up front because array access is faster than
    # computation
    plm, dplm, d2plm = assoc_legendre_func_table(cosθ, ℓ_axes[end])
    sinmφ_vec = OffsetArray(sin.(ℓ_axes * φ), ℓ_axes)
    cosmφ_vec = OffsetArray(cos.(ℓ_axes * φ), ℓ_axes)


    Φ = 0.0
    B = zeros(Float64, 3)
    jacobianB = zeros(Float64, 3, 3)

    jacobianB_ℓm = zeros(Float64, 3, 3)

    for ℓ in ℓ_axes
        # set up F(r). It actually doesn't depend on m, so set it up
        # outside the m loop.
        # first term of numerator
        r_inverse_scaled = (R_0 / r)^(ℓ+1)
        # second term of numerator
        r_positive_scaled = (R_0 / R_SS)^(ℓ+1) * (r / R_SS)^ℓ
        # denominator
        surface_to_surface_scale_denom = ℓ + 1 + ℓ * (R_0/R_SS)^(2ℓ+1)

        F = (r_inverse_scaled - r_positive_scaled) / surface_to_surface_scale_denom
        dF_dr = (
                 -(ℓ+1)/r * r_inverse_scaled - ℓ/r * r_positive_scaled
                ) / surface_to_surface_scale_denom
        d2F_dr2 = ((
                    (ℓ+1) * (ℓ+2) / r^2 * r_inverse_scaled
                    -
                    ℓ * (ℓ-1) / r^2 * r_positive_scaled
                   )
                   /
                   surface_to_surface_scale_denom)

        # XXX Assuming that m starts at 0. how do we rectify?
        for m in 0:ℓ
            # XXX can move this out of the loop?
            G = plm[ℓ,m]
            dG_dθ = -sinθ * dplm[ℓ,m]
            dGm_dθ2 = - sinθ^2 * d2plm[ℓ,m] - cosθ * dplm[ℓ,m]

            H = g[ℓ,m] * cosmφ_vec[m] + h[ℓ,m]  * sinmφ_vec[m]
            dH_dφ = m * (-g[ℓ,m] * sinmφ_vec[m] + h[ℓ,m] * cosmφ_vec[m])
            d2H_dφ = -m^2 * H

            # Φ = ∑_(ℓ, m) R_0 F G H
            Φ_ℓm = R_0 * F * G * H

            Φ += Φ_ℓm

            # B = -∑_(ℓ,m) R_0 F G H (r̂/F * dF/dr
            #                         + θ̂/(r G) * dG/dθ
            #                         + φ̂/(r sinθ H) * dH/dr)
            @. B += - Φ_ℓm .* [dF_dr / F,
                               dG_dθ / (r * G),
                               dH_dφ / (r * sinθ * H)]

            # I am not even going to attempt writing down the math form of the
            # Jacobian here. See Physical-theory.md.
            # we're only setting the lower triangular part then copying it after
            # the loop. (take advantage of symmetry)
            jacobianB_ℓm .= -Φ_ℓm # wasted on the upper triangular part, but whatever

            # first column
            jacobianB_ℓm[1,1] *= d2F_dr2 / F
            jacobianB_ℓm[2,1] *= dF_dr * dG_dθ / (r * F * G)
            jacobianB_ℓm[3,1] *= dF_dr * dH_dφ / (r * sinθ * F * H)

            # second column
            jacobianB_ℓm[2,2] *= dGm_dθ2 / (r^2 * G)
            jacobianB_ℓm[3,2] *= dG_dθ * dH_dφ / (r^2 * sinθ * G * H)

            # third column
            jacobianB_ℓm[3,3] *= d2H_dφ / (r^2 * sinθ^2 * H)

            jacobianB .+= jacobianB_ℓm

        end
    end

    # copy over the upper triangular part of the jacobian
    jacobianB[1,2] = jacobianB[2,1]
    jacobianB[1,3] = jacobianB[3,1]
    jacobianB[2,3] = jacobianB[3,2]

    ∇normB = jacobianB' * B / norm(B)

    return (
            Φ = Φ,
            B = B,
            jacobianB = jacobianB,
            ∇normB = ∇normB,
           )
end

function collectmagneticfield(
        rs::AbstractVector,
        θs::AbstractVector,
        φs::AbstractVector,
        g::AbstractMatrix,
        h::AbstractMatrix
    )

    results = Array{BField}(undef, length(rs), length(θs), length(φs))

    for (iθ, θ) in enumerate(θs)
        for (ir, r) in enumerate(rs), (iφ, φ) in enumerate(φs)
            results[ir,iθ,iφ] = magneticfield(r, θ, φ, g, h)
        end
    end

    return results
end

Plm′(x, l, m; kwargs...) = derivativefd(x -> Plm(x, l, m; kwargs...), x)
Plm′′(x, l, m; kwargs...) = derivativefd(x -> Plm′(x, l, m; kwargs...), x)

@memoize ThreadSafeDict function assoc_legendre_func_table(x, lmax::Integer)
    # declare the three matrices
    # the assoc legendre polynomials are only defined as m ≤ l, so enforce that
    # in a lower triangular matrix.
    plm_mat = OffsetMatrix(LowerTriangular(zeros(lmax+1, lmax+1)), 0:lmax, 0:lmax)
    dplm_mat = OffsetMatrix(LowerTriangular(zeros(lmax+1, lmax+1)), 0:lmax, 0:lmax)
    d2plm_mat = OffsetMatrix(LowerTriangular(zeros(lmax+1, lmax+1)), 0:lmax, 0:lmax)

    for ℓ = 0:lmax
        for m = 0:ℓ
            plm_mat[ℓ,m] = Plm(x, ℓ, m, norm=QSNORM)
            dplm_mat[ℓ,m] = Plm′(x, ℓ, m, norm=QSNORM)
            d2plm_mat[ℓ,m] = Plm′′(x, ℓ, m, norm=QSNORM)
        end
    end

    # return all 3 matrices
    return (plm_mat, dplm_mat, d2plm_mat)
end

end
