module MagneticFieldSHTCExpander

export magneticfield, collectmagneticfield

using OffsetArrays: OffsetArray, OffsetMatrix
using LinearAlgebra: LowerTriangular, norm
using LegendrePolynomials: Plm
using ForwardDiff: derivative as derivativefd

using Memoization: @memoize
using ThreadSafeDicts: ThreadSafeDict

using StaticArrays: SVector, SMatrix

using DocStringExtensions: TYPEDEF, FIELDS

const R₀ = 1 # sun's radius in solar radius
const Rₛₛ = 2.5 # source surface in solar radius (PFSS)

const QSNORM = Val(:schmidtquasi)

# declare the return type as its own struct thing because it makes
# it easier during development. once everything has been solidified,
# this should be deleted and the returns made concrete.
"""
$(TYPEDEF)

Composite data type containing information about the local magnetic field vector.
Note that position information needed for `magneticfield`, ``(r, θ, φ)`` is not stored.

# Fields
$(FIELDS)
"""
Base.@kwdef struct BField
    "magnetic potential at ``(r, θ, φ)``"
    Φ::Float64

    "magnetic field at ``(r, θ, φ)``"
    B::SVector{3,Float64}

    "Jacobian matrix of magnetic field at ``(r, θ, φ)``"
    jacobianB::SMatrix{3,3,Float64,9}
end

function Base.show(io::IO, bfield::BField)
    println(io, "BField(")
    println(io, "  Φ = ", bfield.Φ, ",")
    println(io, "  B = ", bfield.B, ",")
    println(io, "  jacobianB = ", bfield.jacobianB, ",")
    print(io, ")")
end

"""
    magneticfield(r, θ, φ, g, h) -> BField

Return the magnetic field at ``(r, θ, φ)`` as described by ``g_ℓ^m`` and ``h_ℓ^m``.

# Arguments
- `r`: radial distance form the origin (spherical coordinates)
- `θ`: polar angle (ISO/physics spherical coordinates)
- `φ`: azimuthal angle (ISO/physics spherical coordinates)
- `g`, `h`: arrays (matrices) containing the spherical harmonic transform
            coefficients (SHTC). `g[ℓ,m]` should yield ``g_ℓ^m`` and
            `h[ℓ,m]` should yield ``h_ℓ^m``

# See also
[`collectmagneticfield`](@ref)
"""
function magneticfield(
        r, θ, φ,            # position in spherical coordinates (r is in solar radius)
        g::AbstractMatrix,  # indices = (0:ℓmax, 0:ℓmax)
        h::AbstractMatrix;  # indices = (0:ℓmax, 0:ℓmax)
        legendre_cache = assoc_legendre_func_table(cos(θ), axes(g)[1][end]),
    )::BField


    if axes(g) != axes(h)
        throw(DimensionMismatch("g and h must be the same size and index"))
    end

    checkbounds(r, θ, φ)

    # set up array of allowed ℓ values
    ℓ_axes = axes(g, 1)
    m_axes = axes(g, 2)

    # PHYSICAL THEORY (see Physical-theory.md)
    # Each ℓ,m term of the potential can be written as a product of three
    # functions depending uniquely on r, θ, φ each.

    cosθ = cos(θ)
    sinθ = sin(θ)
    if sinθ ≈ 0
        @warn("θ is near 0 or π, causing sin(θ) to be almost 0. This may cause NaN issues")
    end

    # create vectors and matrices up front because array access is faster than
    # computation
    #plm, dplm, d²plm = assoc_legendre_func_table(cosθ, ℓ_axes[end])
    plm, dplm, d²plm = legendre_cache
    sinmφ_vec = OffsetArray(sin.(ℓ_axes * φ), ℓ_axes)
    cosmφ_vec = OffsetArray(cos.(ℓ_axes * φ), ℓ_axes)


    Φ = 0.0
    B = zeros(Float64, 3)
    jacobianB = zeros(Float64, 3, 3)

    jacobianB_ℓm = zeros(Float64, 3, 3)

    for ℓ in ℓ_axes
        # set up F(r). It actually doesn't depend on m, so set it up
        # outside the m loop.
        #
        # first term of numerator
        r_inverse_scaled = (R₀ / r)^(ℓ+1)
        # second term of numerator
        r_positive_scaled = (R₀ / Rₛₛ)^(ℓ+1) * (r / Rₛₛ)^ℓ
        # denominator
        surface_to_surface_scale_denom = ℓ + 1 + ℓ * (R₀/Rₛₛ)^(2ℓ+1)

        # using prime notation for derivatives
        F = (r_inverse_scaled - r_positive_scaled) / surface_to_surface_scale_denom
        F′ = (
              -(ℓ+1)/r * r_inverse_scaled - ℓ/r * r_positive_scaled
             ) / surface_to_surface_scale_denom
        F″ = ((
                (ℓ+1) * (ℓ+2) / r^2 * r_inverse_scaled
                -
                ℓ * (ℓ-1) / r^2 * r_positive_scaled
               )
               /
               surface_to_surface_scale_denom)

        for m in m_axes
            # XXX can move this out of the loop?
            G   = plm[ℓ,m]
            G′  = -sinθ * dplm[ℓ,m]
            G″ = - sinθ^2 * d²plm[ℓ,m] - cosθ * dplm[ℓ,m]

            H   = g[ℓ,m] * cosmφ_vec[m] + h[ℓ,m]  * sinmφ_vec[m]
            H′  = m * (-g[ℓ,m] * sinmφ_vec[m] + h[ℓ,m] * cosmφ_vec[m])
            H″ = -m^2 * H

            # Φ = ∑_(ℓ, m) R₀ F G H
            Φ_ℓm = R₀ * F * G * H

            Φ += Φ_ℓm

            # B = -∑_(ℓ,m) R₀ F G H (r̂/F * F′ + θ̂/(r G) * G′ + φ̂/(r sinθ H) * H′)
            @. B += - Φ_ℓm .* [F′/F, G′/(r * G), H′/(r * sinθ * H)]

            # I am not even going to attempt writing down the math form of the
            # Jacobian here. See Physical-theory.md.
            # we're only setting the lower triangular part then copying it after
            # the loop. (take advantage of symmetry)
            jacobianB_ℓm .= -Φ_ℓm # wasted on the upper triangular part, but whatever

            # first column
            jacobianB_ℓm[1,1] *= F″ / F
            jacobianB_ℓm[2,1] *= F′ * G′ / (r * F * G)
            jacobianB_ℓm[3,1] *= F′ * H′ / (r * sinθ * F * H)

            # second column
            jacobianB_ℓm[2,2] *= G″ / (r^2 * G)
            jacobianB_ℓm[3,2] *= G′ * H′ / (r^2 * sinθ * G * H)

            # third column
            jacobianB_ℓm[3,3] *= H″ / (r^2 * sinθ^2 * H)

            jacobianB .+= jacobianB_ℓm

        end
    end

    # copy over the upper triangular part of the jacobian
    jacobianB[1,2] = jacobianB[2,1]
    jacobianB[1,3] = jacobianB[3,1]
    jacobianB[2,3] = jacobianB[3,2]

    return BField(; Φ, B, jacobianB)
end

"""
    magneticfield(rvec, g, h)

Call `magneticfield` using a `Vector` instead of specifying coordinates.

`rvec` should be a 3-vector containing the ``(r, θ, φ)`` components.
"""
function magneticfield(rvec, g, h)
    return magneticfield(rvec..., g, h)
end

"""
    collectmagneticfield(rs, θs, φs, g, h) -> Array{BField,3}

Evaluate `magneticfield` for multiple position vectors. Returns a 3d array
containing `BField` information, where each index corresponds to the index of
the `rs`, `θs` and `φs` vectors.

# Arguments
- `rs`: vector of radial distances form the origin (spherical coordinates)
- `θs`: vector of polar angles (ISO/physics spherical coordinates)
- `φs`: vector of azimuthal angles (ISO/physics spherical coordinates)
- `g`, `h`: arrays (matrices) containing the spherical harmonic transform
            coefficients (SHTC). `g[ℓ,m]` should yield ``g_ℓ^m`` and
            `h[ℓ,m]` should yield ``h_ℓ^m``

# Examples

```julia
bgrid = collectmagneticfield(rs, θs, φs, g, h)
r = rs[i]
θ = θs[j]
φ = φs[k]
magneticfield(r, θ, φ, g, h) == bgrid[i,j,k]
```

# See also
[`magneticfield`](@ref)
"""
function collectmagneticfield(
        rs::AbstractVector, θs::AbstractVector, φs::AbstractVector,
        g::AbstractMatrix, h::AbstractMatrix,
    )

    results = Array{BField}(undef, length(rs), length(θs), length(φs))

    for (iθ, θ) in enumerate(θs)
        legendre_cache = assoc_legendre_func_table(cos(θ), ℓ_axes[end])
        for (ir, r) in enumerate(rs), (iφ, φ) in enumerate(φs)
            results[ir,iθ,iφ] = magneticfield(r, θ, φ, g, h; legendre_cache)
        end
    end

    return results
end

Plm′(x, l, m; kwargs...) = derivativefd(x -> Plm(x, l, m; kwargs...), x)
Plm″(x, l, m; kwargs...) = derivativefd(x -> Plm′(x, l, m; kwargs...), x)

@memoize ThreadSafeDict function assoc_legendre_func_table(x, ℓmax::Integer)
    # declare the three matrices
    # the assoc legendre polynomials are only defined as m ≤ l, so enforce that
    # in a lower triangular matrix.
    plm_mat   = OffsetMatrix(LowerTriangular(zeros(ℓmax+1, ℓmax+1)), 0:ℓmax, 0:ℓmax)
    dplm_mat  = OffsetMatrix(LowerTriangular(zeros(ℓmax+1, ℓmax+1)), 0:ℓmax, 0:ℓmax)
    d²plm_mat = OffsetMatrix(LowerTriangular(zeros(ℓmax+1, ℓmax+1)), 0:ℓmax, 0:ℓmax)

    for ℓ in 0:ℓmax
        for m in 0:ℓ
            plm_mat[ℓ,m]   = Plm( x, ℓ, m, norm=QSNORM)
            dplm_mat[ℓ,m]  = Plm′(x, ℓ, m, norm=QSNORM)
            d²plm_mat[ℓ,m] = Plm″(x, ℓ, m, norm=QSNORM)
        end
    end

    # return all 3 matrices
    return (plm_mat, dplm_mat, d²plm_mat)
end

function checkbounds(r, θ, φ)
    # ensure radius within solar surface and source surface
    if r < R₀ || r > Rₛₛ
        throw(DomainError(r, "r must be within [R₀, Rₛₛ] = [$R₀, $Rₛₛ]"))
    end

    if θ < 0 || θ > π
        throw(DomainError(θ, "θ must be within [0, π]"))
    end

    if φ < 0 || φ > 2π
        throw(DomainError(φ, "φ must be within [0, 2π]"))
    end
end

end
