using Test: @test, @testset
using BioSequencesOneHot: onehot, potts, alphabet, aaseq, rnaseq, dnaseq
using LinearAlgebra: I
using BioSequences: @dna_str, @rna_str, @aa_str, LongAA, LongRNA, LongDNA

@testset "onehot" begin
    @test onehot(alphabet(LongAA)) == I
    @test onehot(alphabet(LongDNA)) == I
    @test onehot(alphabet(LongRNA)) == I

    @test aaseq(BitMatrix(I(21))) == alphabet(LongAA)
    @test aaseq(BitMatrix(I(20))) == alphabet(LongAA; gap=false)

    @test rnaseq(BitMatrix(I(5))) == alphabet(LongRNA)
    @test rnaseq(BitMatrix(I(4))) == alphabet(LongRNA; gap=false)

    @test dnaseq(BitMatrix(I(5))) == alphabet(LongDNA)
    @test dnaseq(BitMatrix(I(4))) == alphabet(LongDNA; gap=false)

    @test aaseq(BitArray(stack([I(21), I(21)]))) == [alphabet(LongAA), alphabet(LongAA)]
    @test rnaseq(BitArray(stack([I(5), I(5)]))) == [alphabet(LongRNA), alphabet(LongRNA)]
    @test dnaseq(BitArray(stack([I(5), I(5)]))) == [alphabet(LongDNA), alphabet(LongDNA)]
end

@testset "potts" begin
    @test potts(alphabet(LongAA)) == 1:21
    @test potts(alphabet(LongDNA)) == 1:5
    @test potts(alphabet(LongRNA)) == 1:5

    @test potts(alphabet(LongAA; gap=false)) == 1:20
    @test potts(alphabet(LongDNA; gap=false)) == 1:4
    @test potts(alphabet(LongRNA; gap=false)) == 1:4

    @test aaseq(Int8.(1:20)) == alphabet(LongAA; gap=false)
    @test aaseq(Int8.(1:21)) == alphabet(LongAA)

    @test rnaseq(Int8.(1:4)) == alphabet(LongRNA; gap=false)
    @test rnaseq(Int8.(1:5)) == alphabet(LongRNA)

    @test dnaseq(Int8.(1:4)) == alphabet(LongDNA; gap=false)
    @test dnaseq(Int8.(1:5)) == alphabet(LongDNA)
end
