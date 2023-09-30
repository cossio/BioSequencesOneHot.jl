using Test: @test, @testset
using BioSequences: @dna_str, @rna_str, @aa_str, LongRNA
using BioSequencesOneHot: alphabet

@testset "alphabet" begin
    @test alphabet(dna"A"; gap=true) == dna"ACGT-"
    @test alphabet(rna"A"; gap=true) == rna"ACGU-"
    @test alphabet(aa"A"; gap=true) == aa"ACDEFGHIKLMNPQRSTVWY-"

    @test LongRNA{4}(alphabet(dna"ACGT"; gap=true)) == alphabet(LongRNA{4}(dna"ACGT"); gap=true)
end
