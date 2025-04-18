### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 162d6060-1be4-11f0-23bf-8ffa3c0a6091
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ 411c6323-3d37-42f5-b62d-7d422138d878
using DelimitedFiles

# ╔═╡ 4b8cd9d6-45ad-4d06-8cfe-360f2f8452a2
using OffsetArrays

# ╔═╡ c8c6b71b-9645-426c-afe1-5ee00654b4b2
using Unitful

# ╔═╡ ff935b93-6105-45c4-84ac-adcaacb88eb9
using MagneticFieldSHTCExpander

# ╔═╡ 56fe0d3a-d157-4b08-816d-5ceb08450dcb
using WGLMakie

# ╔═╡ 3991df8f-a7f8-4997-a39c-4ee2b36267b5
ℓ_all, m_all, g_vec, h_vec = let
	big_mat = readdlm("gh20.rad_CT1941_170_SN6122", skipstart=15) |> eachcol
	Vector{Int}(big_mat[1]), Vector{Int}(big_mat[2]), big_mat[3]*u"μT", big_mat[4]*u"μT"
end

# ╔═╡ eb80c1a6-f275-4d5e-81bb-0927f58b6cf1
g, h = let
	g = zeros(eltype(g_vec), 0:20, 0:20)
	h = zero(g)
	for (ℓ, m, g_ℓm, h_ℓm) in zip(ℓ_all, m_all, g_vec, h_vec)
		g[ℓ,m] = g_ℓm
		h[ℓ,m] = h_ℓm
	end
	g,h
end

# ╔═╡ e7f3c1d9-a9e7-48df-aa35-dd7459d9cfb8
radius_range = [1, 3//2, 2];

# ╔═╡ a9787825-eabd-4642-814d-304e5c0fe436
theta_range = 1:179 .|> deg2rad;

# ╔═╡ e5e5254c-82db-4e91-b2c7-d820f922cc27
phi_range = 1:359 .|> deg2rad;

# ╔═╡ 83bed572-4fe2-4f9e-9653-5ba3ef7ff2ac
B = collectmagneticfield(radius_range, theta_range, phi_range, g, h)

# ╔═╡ 4bfc6eb7-210d-4252-bcad-a4b49fa8ade3
magpots = getfield.(B, :Φ)

# ╔═╡ fc7bf94d-2cca-4787-913a-7c9857f22f39
plotting_slice = view(magpots,1,:,:)

# ╔═╡ 03019953-b799-415c-a63b-449638640021


# ╔═╡ 5cafd6b5-4076-4f86-ade7-32d2a04c55fb
let fig = Figure()
	ax = Axis3(fig[1,1])

	x = [sin(θ) * sin(ϕ) for θ in theta_range, ϕ in phi_range]
	y = [sin(θ) * cos(ϕ) for θ in theta_range, ϕ in phi_range]
	z = [cos(θ) for θ in theta_range, ϕ in phi_range]

	surface!(ax, x, y, z, color = ustrip.(u"μT", plotting_slice), shading = NoShading)

	fig
end

# ╔═╡ Cell order:
# ╠═162d6060-1be4-11f0-23bf-8ffa3c0a6091
# ╠═411c6323-3d37-42f5-b62d-7d422138d878
# ╠═4b8cd9d6-45ad-4d06-8cfe-360f2f8452a2
# ╠═c8c6b71b-9645-426c-afe1-5ee00654b4b2
# ╠═3991df8f-a7f8-4997-a39c-4ee2b36267b5
# ╠═eb80c1a6-f275-4d5e-81bb-0927f58b6cf1
# ╠═ff935b93-6105-45c4-84ac-adcaacb88eb9
# ╠═e7f3c1d9-a9e7-48df-aa35-dd7459d9cfb8
# ╠═a9787825-eabd-4642-814d-304e5c0fe436
# ╠═e5e5254c-82db-4e91-b2c7-d820f922cc27
# ╠═83bed572-4fe2-4f9e-9653-5ba3ef7ff2ac
# ╠═4bfc6eb7-210d-4252-bcad-a4b49fa8ade3
# ╠═56fe0d3a-d157-4b08-816d-5ceb08450dcb
# ╠═fc7bf94d-2cca-4787-913a-7c9857f22f39
# ╠═03019953-b799-415c-a63b-449638640021
# ╠═5cafd6b5-4076-4f86-ade7-32d2a04c55fb
