module Flatter
using Hecke

function proj(v::Vector{ArbFieldElem}, u::Vector{ArbFieldElem})::Vector{ArbFieldElem}
    retval = dot(v,u) / dot(u,u) .* u
    return retval
end

@doc raw"""
    gso_with_prec(v::Matrix{<:Number}, n::Int)::ArbMatrix

Function returns the result of Gram Schmidt Orthogonalization of `v`, in `n` precision.

# Example
```jlcon
julia> m = [[1 0 331 303]
       [0 1 456 225]
       [0 0 628 0]
       [0 0 0 628]]
4×4 Matrix{Int64}:
 1  0  331  303
 0  1  456  225
 0  0  628    0
 0  0    0  628

julia> u = Flatter.gso_with_prec(v, 64);

julia> u * transpose(matrix(RR, m));

julia> Float64.(ans)
[    201371.0       219111.0       207868.0   190284.0]
[ 1.82424e-14        20148.2        60187.6   -65747.3]
[-2.20325e-14   -4.51653e-14        13.8453   -19.7221]
[ 5.17197e-14     1.8822e-15   -1.57434e-13    2.76887]
```
"""
function gso_with_prec(v::Matrix{<:Number}, n::Int)::ArbMatrix
  RR = ArbField(n)
  return gso_with_prec(matrix(RR,v), RR)
end

function gso_with_prec(v::ArbMatrix, RR::ArbField)::ArbMatrix
    u = Any[]
    for r in axes(v, 1)
        s = zeros(length(v[r,:]))
        for j in 1:r-1
            s = s + proj(v[r,:], u[j])
        end
        push!(u, v[r,:]-s)
    end
    return matrix(RR, mapreduce(permutedims, vcat, u))
end

function proj(v::Vector{QQFieldElem}, u::Vector{QQFieldElem})::Vector{QQFieldElem}
    retval = dot(v,u) / dot(u,u) * u
    return retval
end

function gram_schmidt_orthogonalization(v::QQMatrix)::Any
    u = Any[]
    for r in axes(v, 1)
        s = zeros(length(v[r,:]))
        for j in 1:r-1
            s = s + proj(v[r,:], u[j])
        end
        push!(u, v[r,:]-s)
    end
    return matrix(QQ,mapreduce(permutedims, vcat, u))
end


@doc raw"""
    norm(b::Vector{ArbFieldElem}) -> ArbFieldElem

Function returns the 2norm of a vector, for gso profile reasons
"""
function norm(b::Vector{ArbFieldElem})
    return sqrt(sum(b.*b))::ArbFieldElem
end

@doc raw"""
    profile(gso::ArbMatrix) -> Vector{ArbFieldElem}

Function gives profile of a given prec GSO matrix (row major)
"""
function profile(gso::ArbMatrix)
    l = ArbFieldElem[]
    for i in 1:nrows(gso)
        push!(l, norm(gso[i,:]))
    end
    return log.(2, BigFloat.(l))::Vector{BigFloat}
end

function spread(profile::Vector{ArbFieldElem})
  return maximum(profile) - minimum(profile)
end

function potential(profile::Vector{ArbFieldElem})
  n = length(profile)
  a = 0
  for i in 1:n
    a += (n - i + 1) * profile[i]
  end
  return a
end


@doc raw"""
Profile of a lattice is defined as the volume of a special set D where

D = union of all [l_next, l_current] where l_next < l_current

"""
function drop(v)
    cur = curmin = curmax = v[1]
    curstart = 1
    dropset = []
    for i in range(2,length(v))
        next = v[i]
        if (next < cur) || (next < curmax)
            curmin = min(curmin, next)
        else
            push!(dropset, v[curstart:i-1])
            curmax = max(next, curmax)
            curstart = i
        end
        cur = next

        if i == length(v)
            # end of loop
            push!(dropset, v[curstart:i])
        end
    end
    return dropset
end


end