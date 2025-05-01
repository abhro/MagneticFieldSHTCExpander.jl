# functions to calculate the magnetic field at multiple points
# i.e., all the `collectmagneticfield` methods

"""
    collectmagneticfield(rs, θs, φs, g, h) -> Array{BField,3}

Evaluate `magneticfield` for multiple position vectors. Returns a 3d array
containing `BField` information, where each index corresponds to the index of
the `rs`, `θs` and `φs` vectors.

### Arguments
- `rs`: vector of radial distances form the origin (spherical coordinates)
- `θs`: vector of polar angles (ISO/physics spherical coordinates)
- `φs`: vector of azimuthal angles (ISO/physics spherical coordinates)
- `g`, `h`: arrays (matrices) containing the spherical harmonic transform
            coefficients (SHTC). `g[ℓ,m]` should yield ``g_ℓ^m`` and
            `h[ℓ,m]` should yield ``h_ℓ^m``

### Examples

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
    ℓmax = last(axes(g, 1))

    for (iθ, θ) in enumerate(θs)
        legendre_cache = assoc_legendre_func_table(cos(θ), ℓmax)
        for (ir, r) in enumerate(rs)
            collectmagneticfield!(view(results, ir, iθ, :), r, θ, φs, g, h; legendre_cache)
        end
    end

    return results
end

# only loop over θ and φ
function collectmagneticfield(
        r, θs::AbstractVector, φs::AbstractVector,
        g::AbstractMatrix, h::AbstractMatrix,
    )

    results = Matrix{BField}(undef, length(θs), length(φs))

    for (iθ, θ) in enumerate(θs)
        collectmagneticfield!(view(results, iθ, :), r, θ, φs, g, h)
    end

    return results
end

# only loop over φ
"""
Essentially `[magneticfield(r, θ, φ) for φ in φs]`
"""
function collectmagneticfield(
        r, θ, φs::AbstractVector,
        g::AbstractMatrix, h::AbstractMatrix;
        legendre_cache = assoc_legendre_func_table.(cos(θ), axes(g)[1][end]),
    )
    return collectmagneticfield!(Vector{BField}(undef, length(φs)), r, θ, φs, g, h; legendre_cache)
end
function collectmagneticfield!(
        results::AbstractVector,
        r, θ, φs::AbstractVector,
        g::AbstractMatrix, h::AbstractMatrix;
        legendre_cache = assoc_legendre_func_table.(cos(θ), axes(g)[1][end]),
    )

    length(results) == length(φs) || throw(DimensionMismatch())
    ℓmax = last(axes(g, 1))

    for (i, φ) in enumerate(φs)
        results[i] = magneticfield(r, θ, φ, g, h; legendre_cache)
    end

    return results
end
