# Author: Jef Pauwels
# Last Modified: 2023-05-01
# Description: This script computes upper bounds on a linear witness for the EACC scenario using the NPA hierarchy.
# The user can specify adaptive vs nonadaptive strategies, the dimension of the communication, and the witness.
# (see https://arxiv.org/abs/2203.05372 for details).

# Requirements: Julia package QuantumNPA.jl,  a Julia implementation of the NPA hierarchy 
# from https://github.com/ewoodhead/QuantumNPA.jl

# Solves the SDP using Mosek.


using Pkg; Pkg.add(["Combinatorics", "JuMP", "SCS", "BlockDiagonals"])

include("QuantumNPA.jl/QuantumNPA.jl")


using .QuantumNPA

using Mosek
using MosekTools

# Basic parameters

d = 3;  # dimension communication
adaptive = false; # adaptive or nonadaptive strategies

# Specify witness, in the format c[b,x,y]

# 2 --> 1 RAC c[b,x,y] = delta_{b x_y}

nX = 4; # number of inputs for Alice
nY = 2; # number of inputs for Bob
nB = 2; # number of outputs for Bob

c = zeros(2, 4, 2)
c[1, :, :] = [1 1; 1 0; 0 1; 0 0] / 8
c[2, :, :] = [0 0; 0 1; 1 0; 1 1] / 8


# Decoding functions 

if adaptive == false
v = Vector{Array{Int,1}}(undef, nB^d)
for i in 1:nB^d
    indices = zeros(Int, d)
    for j in 1:d
        indices[j] = mod(i-1, nB^(d-j+1)) ÷ nB^(d-j) + 1
    end
    v[i] = indices
end
else
count = zeros(Int, nY, d)
for y in 1:nY
    for a in 1:d
        count[y, a] = (y-1)*d + a
    end
end
end

# Create measurement projectors for Alice and Bob

if adaptive == false
PA = projector(1, 1:d, 1:nX, full=true)
PB = projector(2, 1:nB^d, 1:nY, full=true)
else
PA = projector(1, 1:d, 1:nX, full=true)
PB = projector(2, 1:nB,1:nY*d,full = true)
end

# Create the Bell operator

if adaptive == false
W = sum([c[b,x,y] * PA[m,x] * PB[k,y] for x = 1:nX, y = 1:nY, b = 1:nB, m = 1:d, k = 1:nB^d if v[k][m] == b])
else
W = sum([c[b,x,y] * PA[m,x] * PB[b,count[y,m]] for x = 1:nX, y = 1:nY, b = 1:nB, m = 1:d])
end

# Solve the NPA hierarchy

npa_max(W,3,solver=Mosek.Optimizer)



