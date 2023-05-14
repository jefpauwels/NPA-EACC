# This script computes upper bounds on a linear witness for the EACC scenario using the NPA hierarchy.
# It is based on the Julia package QuantumNPA.jl, which is a Julia implementation of the NPA hierarchy 
# from https://github.com/ewoodhead/QuantumNPA.jl
# Solves the SDP using Mosek.

include("QuantumNPA.jl/qnpa.jl")



using .QuantumNPA

using Mosek
using MosekTools

# Basic parameters

d = 2;  # dimension communication

nX = 4; # number of inputs for Alice
nY = 2; # number of inputs for Bob
nB = 2; # number of outputs for Bob

# Decoding functions 

if d == 2
    v = Vector{Array{Int,1}}(undef, nB^2)
    for i in 1:nB
        for j in 1:nB
            v[(i-1)*nB+j] = [i, j]
        end
    end
elseif d == 3
    v = Vector{Array{Int,1}}(undef, nB^3)
    for i in 1:nB
        for j in 1:nB
            for k in 1:nB
                v[(i-1)*nB^2+(j-1)*nB+k] = [i, j, k]
            end
        end
    end
end

# Specify witness, in the format c[b,x,y]

# 2 --> 1 RAC c[b,x,y] = delta_{b x_y}

c = zeros(2, 4, 2)
c[1, :, :] = [1 1; 1 0; 0 1; 0 0] / 8
c[2, :, :] = [0 0; 0 1; 1 0; 1 1] / 8

# Create measurement projectors for Alice and Bob

PA = projector(1, 1:d, 1:nX, full=true)
PB = projector(2, 1:nB^d, 1:nY, full=true)

# Create the Bell operator


W = sum([c[b,x,y] * PA[m,x] * PB[k,y] for x = 1:nX, y = 1:nY, b = 1:nB, m = 1:d, k = 1:nB^d if v[k][m] == b])

# Solve the NPA hierarchy

npa_max(W,2,solver=Mosek.Optimizer)



