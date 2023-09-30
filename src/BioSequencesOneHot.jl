module BioSequencesOneHot

using BioSequences: LongRNA, LongDNA, LongAA, @rna_str, @aa_str, @dna_str

include("alphabet.jl")
include("potts.jl")
include("onehot.jl")

end
