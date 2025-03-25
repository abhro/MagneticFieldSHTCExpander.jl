using MagneticFieldSHTCExpander
using Test
using OffsetArrays
using Aqua

@testset "MagneticFieldSHTCExpander" begin

    zerog = zeros(0:12, 0:12)
    zeroh = zeros(0:12, 0:12)
    @test magneticfield(1.4, π/3, 0.8*π, zerog, zeroh).Φ == 0

    @testset "Code quality (Aqua)" begin
        Aqua.test_all(MagneticFieldSHTCExpander)
    end
end
