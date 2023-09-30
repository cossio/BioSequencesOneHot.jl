#=
Define the order in which amino-acids and nucleotides are encoded as integers.
=#

function alphabet(::Type{<:LongAA}; gap::Bool = true)
    AAs = aa"ACDEFGHIKLMNPQRSTVWY"
    if gap
        return AAs * aa"-"
    else
        return AAs
    end
end

function alphabet(::Type{<:LongRNA}; gap::Bool = true)
    NTs = rna"ACGU"
    if gap
        return NTs * rna"-"
    else
        return NTs
    end
end

function alphabet(::Type{<:LongDNA}; gap::Bool = true)
    NTs = dna"ACGT"
    if gap
        return NTs * rna"-"
    else
        return NTs
    end
end

function alphabet(s::Union{<:LongAA, <:LongRNA{4}, <:LongDNA{4}}; gap::Bool = true)
    return alphabet(typeof(s); gap)
end
