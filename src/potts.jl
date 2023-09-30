#=
Encodes sequences as Potts arrays, where the integer entry indicates the sequence letter.
=#

function potts(X::BitMatrix)
    return vec(Int8.(first.(Tuple.(argmax(X; dims=1)))))
end

function potts(X::BitArray{3})
    return reshape(Int8.(first.(Tuple.(argmax(X; dims=1))), size(X,2), size(X,3)))
end

function potts(s::Union{LongAA, LongAA, LongRNA{4}, LongDNA{4}, AbstractVector{<:LongAA}, AbstractVector{<:LongRNA}, AbstractVector{<:LongDNA}})
    potts(onehot(s))
end

function aaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongAA([alphabet(LongAA; gap)[i] for i in P])
end

function rnaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongRNA{4}([alphabet(LongRNA; gap)[i] for i in P])
end

function dnaseq(P::AbstractVector{Int8}; gap::Bool=true)
    return LongDNA{4}([alphabet(LongDNA; gap)[i] for i in P])
end

function aaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [aaseq(view(P,:,n); gap) for n in axes(P, 2)]
end

function rnaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [rnaseq(view(P,:,n); gap) for n in axes(P, 2)]
end

function dnaseq(P::AbstractMatrix{Int8}; gap::Bool=true)
    return [dnaseq(view(P,:,n); gap) for n in axes(P, 2)]
end
