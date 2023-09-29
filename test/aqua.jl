import Aqua
import BioSequencesOneHot

using Test: @testset

@testset verbose = true "aqua" begin
    Aqua.test_all(BioSequencesOneHot; ambiguities = false)
end
