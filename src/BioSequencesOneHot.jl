module BioSequencesOneHot

using BioSequences: LongRNA, LongDNA, @rna_str

function onehot(seq::LongRNA{4})
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(rna"ACGU-")
end

onehot(seq::LongDNA{4}) = onehot(LongRNA{4}(seq))

function onehot(seqs::Union{AbstractVector{<:LongRNA}, AbstractVector{<:LongDNA}})
    L = only(unique(length.(seqs)))
    return reshape(reduce(hcat, onehot.(seqs)), 5, L, length(seqs))
end

end
