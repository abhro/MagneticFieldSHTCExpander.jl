var documenterSearchIndex = {"docs":
[{"location":"physical-theory/#Potential-field","page":"Physical theory","title":"Potential field","text":"","category":"section"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"The magnetic field is based on the Altschuler and Newkirk (1969) and Zhang et al. (2019) model, where a potential field is defined","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Phi(r θ φ)\n= R_odot\n  sum_ℓ=0^infty\n    sum_m=0^ℓ\n      P_ℓ^m(cos θ)\n      leftg_ℓ m cos(mφ) + h_ℓ m sin(m φ)right\n      frac\n        displaystyleleft(fracR_odotrright)^ℓ+1 - left(fracR_odotR_ssright)^ℓ+1 left(fracrR_ssright)^ℓ\n      \n        displaystyle ℓ + 1 + ℓ left(fracR_odotR_ssright)^2ℓ+1\n      ","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"and the magnetic field,","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"mathbfB(r θ φ) = B_r(r θ φ) hatmathbfr + B_θ(r θ φ) hatboldsymbolθ + B_φ(r θ φ) hatboldsymbolφ","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"is the negative gradient of the potential field:","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"mathbfB(r θ φ)\n= - boldsymbolPhi(r θ φ)\n= - frac Phi r hatmathbfr\n  - frac1r frac Phi θ hatboldsymbolθ\n  - frac1r sin θ frac Phi φ hatboldsymbolφ","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Note, here P_ℓ^m is the associated Legendre Polynomial under quasi-Schmidt normalization.","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"For notational convenience, we define the potential function as a product solution, where each function is defined as","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"beginalign*\nF(r) = frac\n    displaystyle left(fracR_odotrright)^ℓ+1 - left(fracR_odotR_ssright)^ℓ+1 left(fracrR_ssright)^ℓ\n\n    displaystyle ℓ + 1 + ℓleft(fracR_odotR_ssright)^2ℓ+1\n \nG(θ) = P_ℓ^m(cos θ) 1ex\nH(φ) = g_ℓ m cos(mφ) + h_ℓ m sin(mφ)\nendalign*","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Note that F, G, and H all implicitly depend on ℓ and m.","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"So, the potential function is instead","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Phi(r θ φ) = R_odot sum_ℓ=0^infty sum_m=0^ℓ F(r) G(θ) H(φ)","category":"page"},{"location":"physical-theory/#Derivatives-of-potential-field","page":"Physical theory","title":"Derivatives of potential field","text":"","category":"section"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"We need the first and second derivatives of each single-variable function for use in the Jacobian of the magnetic field.","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"beginalign*\nfracddr F(r)\n= frac displaystyle\n  -(ℓ+1) fracR_odot^ℓ+1r^ℓ+2 - ℓ left(fracR_odotR_ssright)^ℓ+1 fracr^ℓ-1R_ss^ℓ\n displaystyle\n  ℓ + 1 + ℓ left(fracR_odotR_ssright)^2ℓ+1\n 1ex\n= frac displaystyle\n  -fracℓ+1r left(fracR_odotrright)^ℓ+1 - fracℓr left(fracR_odotR_ssright)^ℓ+1 left(fracrR_ssright)^ℓ\n displaystyle\n  ℓ + 1 + ℓ left(fracR_odotR_ssright)^2ℓ+1\n\n2ex\nfracd^2dr^2 F(r)\n= fracdisplaystyle\n  (ℓ+1)(ℓ+2) fracR_odot^ℓ+1r^ℓ+3 - ℓ(ℓ-1) left(fracR_odotR_ssright)^ℓ+1 fracr^ℓ-2R_ss^ℓ\n displaystyle\n  ℓ + 1 + ℓ left(fracR_odotR_ssright)^2ℓ+1\n 1ex\n= fracdisplaystyle\n  frac(ℓ+1)(ℓ+2)r^2 left(fracR_odotrright)^ℓ+1 - fracℓ(ℓ-1)r^2 left(fracR_odotR_ssright)^ℓ+1 left(fracrR_ssright)^ℓ\n displaystyle\n  ℓ + 1 + ℓ left(fracR_odotR_ssright)^2ℓ+1\n\nendalign*","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"beginalign*\nfracddθ G(θ) = - sin(θ)  fracdP_ℓ^m(cosθ)d(cosθ) 1ex\nfracd^2dθ^2 G(θ) = sin^2(θ)  fracd^2 P_ℓ^m(cosθ)d(cosθ)^2 - cos(θ)  fracdP_ℓ^m(cosθ)d(cosθ)\nendalign*","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"beginalign*\nfracddφ H(φ) = m left-g_ℓ m sin(mφ) + h_ℓ m cos(mφ)right 1ex\nfracd^2dφ^2 H(φ) = -m^2 g_ℓ m cos(mφ) + h_ℓ m sin(mφ) = -m^2 H_ℓ^m(φ)\nendalign*","category":"page"},{"location":"physical-theory/#Jacobian-of-magnetic-field","page":"Physical theory","title":"Jacobian of magnetic field","text":"","category":"section"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"The Jacobian for an mathbbR^3 is a 3×3 matrix. For B, it is","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"beginalign*\nJmathbfB(r θ φ)\n= beginbmatrix\n    displaystyle -frac^2 Phi r^2                        \n    displaystyle -frac1rfrac^2 Phi θ   r \n    displaystyle -frac1r sin θ frac^2 Phi φ   r 05ex\n    displaystyle -frac1r frac^2 Phi r   θ \n    displaystyle -frac1r^2 frac^2 Phi θ^2           \n    displaystyle -frac1r^2 sinθ frac^2 Phiφ  θ  05ex\n    displaystyle -frac1r sin θ frac^2 Phi r   φ \n    displaystyle -frac1r^2 sin θ frac^2 Phiθ  φ \n    displaystyle -frac1r^2 sin^2θ frac^2 Phi φ^2\nendbmatrix\n2ex\n= beginbmatrix\n    longleftarrow  boldsymbol^mathsfT B_r  longrightarrow 05ex\n    longleftarrow  boldsymbol^mathsfT B_θ  longrightarrow 05ex\n    longleftarrow  boldsymbol^mathsfT B_φ  longrightarrow\nendbmatrix\nendalign*","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"In terms of the SHTC expansion, the Jacobian is","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"JmathbfB(r θ φ)\n= R_odot sum_ℓ=0^infty sum_m = 0^ℓ beginbmatrix\n    displaystyle -fracd^2Fdr^2 G H \n    displaystyle -frac1r fracdFdr fracdGdθ H \n    displaystyle -frac1rsinθ fracdFdr G fracdHdφ \n    15ex\n    displaystyle -frac1r fracdFdr fracdGdθ H \n    displaystyle -frac1r^2 F fracd^2Gdθ^2 H \n    displaystyle -frac1r^2sinθ F fracdGdθ fracdHdφ \n    15ex\n    displaystyle -frac1rsinθ fracdFdr G fracdHdφ \n    displaystyle -frac1r^2sinθ F fracdGdθ fracdHdφ \n    displaystyle -frac1r^2sin^2θ F G fracd^2Hdφ^2\nendbmatrix","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Factoring out the parts that appear in Φ, we have","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"JmathbfB(r θ φ)\n= sum_ℓ=0^infty sum_m = 0^ℓ - R_odot F(r) G(θ) H(φ) beginbmatrix\n    displaystyle frac1F fracd^2Fdr^2 \n    displaystyle frac1r frac1F G fracdFdr fracdGdθ \n    displaystyle frac1rsinθ frac1F H fracdFdr fracdHdφ\n    1ex\n    displaystyle frac1r frac1F G fracdFdr fracdGdθ \n    displaystyle frac1r^2 frac1G fracd^2Gdθ^2 \n    displaystyle frac1r^2sinθ frac1G H fracdGdθ fracdHdφ\n    1ex\n    displaystyle frac1r  sinθ frac1F H fracdFdr fracdHdφ \n    displaystyle frac1r^2sinθ frac1G H fracdGdθ fracdHdφ \n    displaystyle frac1r^2sin^2θ frac1H fracd^2Hdφ^2\nendbmatrix","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Note that since B is a gradient field, the Jacobian is symmetric (negative Hessian of Φ).","category":"page"},{"location":"physical-theory/#Gradient-of-the-field-strength","page":"Physical theory","title":"Gradient of the field strength","text":"","category":"section"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"∇B = ∇|B| can be determined from the Jacobian of B. What follows is a proof using index notation:","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"fracBx_j\n= fracx_j (B_i B_i)^12\n= frac12 (B_i B_i)^-12 left(2 B_i fracB_ix_jright)\n= fracB_i fracB_ix_jleft(B_k B_kright)^12\n= fracB_iB fracB_ix_j","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"which implies B = hatmathbfb  JmathbfB","category":"page"},{"location":"physical-theory/#Desired-quantities","page":"Physical theory","title":"Desired quantities","text":"","category":"section"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"The following are the quantities we want to calculate","category":"page"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"The potential field\nΦ = sum_ℓ=0^infty sum_m=0^ℓ R_odot F(r) G(θ) H(φ) = sum_ℓm Phi_ℓ^m\nwhere Phi_ℓ^m = R_odot F(r) G(θ) H(φ)\nThe magnetic field displaystyle mathbfB = - sum_ℓm R_odot F(r) G(θ) H(φ) left(frachatmathbfrF fracdFdr + frachatboldsymbolθr G fracdGdθ + frachatboldsymbolφr sin(θ) H fracdHdvarphiright)\nThe magnetic field strength, B = |B|\nThe gradient of the magnetic field strength, ∇B\nbeginalign*\nboldsymbol B\n= boldsymbol sqrtB_r^2 + B_θ^2 + B_φ^2 1ex\n= frac1B beginbmatrix\n      longleftarrow  dfrac mathbfB^mathsfT r  longrightarrow \n      longleftarrow  dfrac1r dfrac mathbfB^mathsfT θ  longrightarrow \n      longleftarrow  dfrac1r sin θ dfrac mathbfB^mathsfT φ  longrightarrow\n  endbmatrix mathbfB\n1ex\n= frac1B beginbmatrix\n      uparrow  uparrow  uparrow \n      boldsymbol B_r  boldsymbol B_θ  boldsymbol B_φ \n      downarrow  downarrow  downarrow\n  endbmatrix mathbfB 1ex\n= frac1B (JmathbfB)^mathsfT mathbfB\nendalign*","category":"page"},{"location":"physical-theory/#References","page":"Physical theory","title":"References","text":"","category":"section"},{"location":"physical-theory/","page":"Physical theory","title":"Physical theory","text":"Altschuler, M. D. and Newkirk, G. (1969). Magnetic fields and the structure of the solar corona: I: Methods of calculating coronal fields. Solar Physics 9, 131–149.\n\n\n\nZhang, M.; Zhao, L. and Rassoul, H. K. (2019). Stochastic Propagation of Solar Energetic Particles in Coronal and Interplanetary Magnetic Fields. Journal of Physics: Conference Series 1225, 012010.\n\n\n\n","category":"page"},{"location":"#MagneticFieldSHTCExpander.jl","page":"Introduction","title":"MagneticFieldSHTCExpander.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Recreate magnetic field vectors out of spherical harmonic transform coefficients (SHTC).","category":"page"},{"location":"#API-Reference","page":"Introduction","title":"API Reference","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"magneticfield\ncollectmagneticfield\n\nMagneticFieldSHTCExpander.BField","category":"page"},{"location":"#MagneticFieldSHTCExpander.magneticfield","page":"Introduction","title":"MagneticFieldSHTCExpander.magneticfield","text":"magneticfield(r, θ, φ, g, h) -> BField\n\nReturn the magnetic field at (r θ φ) as described by g_ℓ^m and h_ℓ^m.\n\nArguments\n\nr: radial distance form the origin (spherical coordinates)\nθ: polar angle (ISO/physics spherical coordinates)\nφ: azimuthal angle (ISO/physics spherical coordinates)\ng, h: arrays (matrices) containing the spherical harmonic transform           coefficients (SHTC). g[ℓ,m] should yield g_ℓ^m and           h[ℓ,m] should yield h_ℓ^m\n\nSee also\n\ncollectmagneticfield\n\n\n\n\n\nmagneticfield(rvec, g, h)\n\nCall magneticfield using a Vector instead of specifying coordinates.\n\nrvec should be a 3-vector containing the (r θ φ) components.\n\n\n\n\n\n","category":"function"},{"location":"#MagneticFieldSHTCExpander.collectmagneticfield","page":"Introduction","title":"MagneticFieldSHTCExpander.collectmagneticfield","text":"collectmagneticfield(rs, θs, φs, g, h) -> Array{BField,3}\n\nEvaluate magneticfield for multiple position vectors. Returns a 3d array containing BField information, where each index corresponds to the index of the rs, θs and φs vectors.\n\nArguments\n\nrs: vector of radial distances form the origin (spherical coordinates)\nθs: vector of polar angles (ISO/physics spherical coordinates)\nφs: vector of azimuthal angles (ISO/physics spherical coordinates)\ng, h: arrays (matrices) containing the spherical harmonic transform           coefficients (SHTC). g[ℓ,m] should yield g_ℓ^m and           h[ℓ,m] should yield h_ℓ^m\n\nExamples\n\nbgrid = collectmagneticfield(rs, θs, φs, g, h)\nr = rs[i]\nθ = θs[j]\nφ = φs[k]\nmagneticfield(r, θ, φ, g, h) == bgrid[i,j,k]\n\nSee also\n\nmagneticfield\n\n\n\n\n\n","category":"function"},{"location":"#MagneticFieldSHTCExpander.BField","page":"Introduction","title":"MagneticFieldSHTCExpander.BField","text":"struct BField\n\nComposite data type containing information about the local magnetic field vector. Note that position information needed for magneticfield, (r θ φ) is not stored.\n\nFields\n\nΦ: magnetic potential at (r θ φ)\nB: magnetic field at (r θ φ)\njacobianB: Jacobian matrix of magnetic field at (r θ φ)\n∇normB: gradient of norm of magnetic field (= B =  mathbfB)\n\n\n\n\n\n","category":"type"}]
}
