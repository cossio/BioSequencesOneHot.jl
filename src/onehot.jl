function onehot(seq::Union{LongAA, LongRNA{4}, LongDNA{4}}; gap::Bool=true)
    seq_ = collect(seq)
    return reshape(seq_, 1, size(seq_)...) .== collect(alphabet(seq; gap))
end

function onehot(seqs::Union{AbstractVector{<:LongAA}, AbstractVector{<:LongRNA}, AbstractVector{<:LongDNA}}; gap::Bool=true)
    L = only(unique(length.(seqs))) # all sequences must have same length
    return reshape(mapreduce(s -> onehot(s; gap), hcat, seqs), :, L, length(seqs))
end

function aaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 21
        gap = true
    elseif size(X, 1) == 20
        gap = false
    else
        error("Expected 20 or 21 rows; got $(size(X, 1))")
    end

    return aaseq(potts(X); gap)
end

function rnaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 5
        gap = true
    elseif size(X, 1) == 4
        gap = false
    else
        error("Expected 4 or 5 rows; got $(size(X, 1))")
    end

    return rnaseq(potts(X); gap)
end

function dnaseq(X::Union{BitMatrix, BitArray{3}})
    if size(X, 1) == 5
        gap = true
    elseif size(X, 1) == 4
        gap = false
    else
        error("Expected 4 or 5 rows; got $(size(X, 1))")
    end

    return dnaseq(potts(X); gap)
end
