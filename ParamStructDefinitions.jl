using StaticArrays
##############################################################
#Define structs to store parameters for the simulation--------
##############################################################

struct CrowdParams{T_L<:Real,T_Float<:Real}
    Rs::Vector{T_Float} #Vector containing the radii of each person
    Lx::T_L #Length in x-direction of dancefloor
    Ly::T_L #Length in y-direction of dancefloor
    N::Int #Number of people in the crowd
    dt::T_Float #Timestep size
    maxSteps::Int #Maximum number of timesteps
    masses::Vector{T_Float} #Vector containing the masses of each person
    cBoxSize::T_Float #Size of a cell list box (take to be square)
    alpha::T_Float #Neighbour list radius parameter
    maxNeighbours::Int #Maximum neighbour list size
    v₀::T_Float #Limiting speed far from the stage
    boundaryScale::T_Float #Scale for the region near the stage where speed drops

    function CrowdParams{T_L,T_Float}(Rs::Real, Lx, Ly, N, dt, maxSteps, masses, cBoxSize, alpha, maxNeighbours, v₀, boundaryScale) where {T_L<:Real,T_Float<:Real}
        new{T_L,T_Float}(Rs*ones(T_Float,N), Lx, Ly, N, dt, maxSteps, masses, cBoxSize, alpha, maxNeighbours, v₀, boundaryScale)
    end #function
    function CrowdParams{T_L,T_Float}(Rs, Lx, Ly, N, dt, maxSteps, masses::Real, cBoxSize, alpha, maxNeighbours, v₀, boundaryScale) where {T_L<:Real,T_Float<:Real}
        new{T_L,T_Float}(Rs, Lx, Ly, N, dt, maxSteps, masses*ones(T_Float,N), cBoxSize, alpha, maxNeighbours, v₀, boundaryScale)
    end #function
    function CrowdParams{T_L,T_Float}(Rs::Real, Lx, Ly, N, dt, maxSteps, masses::Real, cBoxSize, alpha, maxNeighbours, v₀, boundaryScale) where {T_L<:Real,T_Float<:Real}
        new{T_L,T_Float}(Rs*ones(T_Float,N), Lx, Ly, N, dt, maxSteps, masses*ones(T_Float,N), cBoxSize, alpha, maxNeighbours, v₀, boundaryScale)
    end #function
    
end #struct

function CrowdParams(Rs, Lx, Ly, N, dt, maxSteps, masses, cBoxSize, alpha, maxNeighbours, v₀, boundaryScale) 
    #If cBoxSize does not divide into Lx or Ly, give an error
    if (Lx % cBoxSize, Ly % cBoxSize) != (0,0) 
        error("cBoxSize does not divide Lx and Ly")
    else
        return CrowdParams{typeof(Lx),typeof(dt)}(Rs, Lx, Ly, N, dt, maxSteps, masses, cBoxSize, alpha, maxNeighbours, v₀, boundaryScale)
    end #if
end #function

mutable struct MusicParams{T_Float<:Real}
    noise::T_Float
    moshRate::T_Float
    frontRowPull::T_Float
    angleNoise::T_Float

    function MusicParams{T_Float}(noise::Real,moshRate::Real,frontRowPull::Real,angleNoise::Real) where {T_Float}
        new{T_Float}(T_Float(noise),T_Float(moshRate),T_Float(frontRowPull),T_Float(angleNoise))
    end #function
end #struct
MusicParams(noise, moshRate, frontRowPull, angleNoise) = MusicParams{typeof(noise)}(noise, moshRate, frontRowPull, angleNoise)