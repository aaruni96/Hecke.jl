module Flatter
using Hecke
using LinearAlgebra

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


#=
# size reduction algorithm
# algorithm 5 from paper
=#


function rand_upper_triangular(s,n)
    println("entering rand_upper_triangular")
    return matrix(ZZ, triu(rand(1:s,n,n)))
end

function norm_schur(A)
    println("entering norm_schur")
    n = nrows(A)
    println("calculated n")
    A = A.^2
    println("calculated A^2")
    s = BigFloat(sum(A))
    println("calculated s")
    return sqrt(s)
end

function estimate_condition_number(B)
    println("entering estimate_condition_number")
    # this appears to be *MUCH* larger than the condition number
    S = BigFloat.(transpose(B)*B)
    c = norm_schur(S)*norm_schur(inv(S))
    return c
end

function play()
    a, b = rand(5:100,10)
    println(a)
    println(b)
    println("=========")
    B = rand_upper_triangular(a, b)

    # estimate using schur index
    c = estimate_condition_number(B)

    # straight up calculate from eigen values
    S = transpose(B)*B
    ev = BigFloat.(eigenvalues(RealField(), ZZMatrix(S)))
    k = sqrt(maximum(ev)/minimum(ev))
    println(minimum(ev))
    return c,k, c-k
end

function size_reduce_algo_5(B, c)
    n = nrows(B)
    b = zeros(n,n,n)
    b[:,:,1] .= B
    U = matrix(ZZ, Matrix{Int}(I,n,n))
    pprime = Int(ceil(c + log2(n)) + 2)
    for j in 2:n
        println("j = $j")
        for i in j-1:1
            println("i = $i")
            q = round(b[i,j,j-1] / b[i,i,1])
            for k in 1:i
                println("k = $k")
                b[k,j,j-i+1] = b[k,j,j-i] - (q * b[k,i,i-k+1])
                U[k,j] = U[k,j] - ZZ(q * U[k,i])%(ZZ(2)^pprime)
            end
        end
    end
    Bprime = zeros(n,n)
    for j in 1:n
        for i in 1:j
            Bprime[i,j] = b[i,j,j-i+1]
        end
    end
    Bprime = matrix(ZZ, Bprime)
    return Bprime, U
end

function run()
    a, b = rand(5:100,2)
    B = rand_upper_triangular(a, 3)
    c = estimate_condition_number(B)
    c = Int(round(c))+1
    Bprime, U = size_reduce_algo_5(B,c)
    #check
    if Bprime != B*U
        println("Bprime is $Bprime")
        println("BU is $(B*U)")
        error("Bprime should be equal to BU")
    end
    return Bprime, U
end


end