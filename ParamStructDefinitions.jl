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

    function CrowdParams{T_L,T_Float}(Rs::Real, Lx, Ly, N, dt, maxSteps, masses) where {T_L<:Real,T_Float<:Real}
        new{T_L,T_Float}(Rs*ones(T_Float,N), Lx, Ly, N, dt, maxSteps, masses)
    end #function
    function CrowdParams{T_L,T_Float}(Rs, Lx, Ly, N, dt, maxSteps, masses::Real) where {T_L<:Real,T_Float<:Real}
        new{T_L,T_Float}(Rs, Lx, Ly, N, dt, maxSteps, masses*ones(T_Float,N))
    end #function
    function CrowdParams{T_L,T_Float}(Rs::Real, Lx, Ly, N, dt, maxSteps, masses::Real) where {T_L<:Real,T_Float<:Real}
        new{T_L,T_Float}(Rs*ones(T_Float,N), Lx, Ly, N, dt, maxSteps, masses*ones(T_Float,N))
    end #function
    
end #struct
CrowdParams(Rs, Lx, Ly, N, dt, maxSteps, masses) = CrowdParams{typeof(Lx),typeof(dt)}(Rs, Lx, Ly, N, dt, maxSteps, masses)

struct MusicParams{T_Float<:Real}
    noise::T_Float
    moshRate::T_Float
    frontRowPull::T_Float

    function MusicParams{T_Float}(noise::Real,moshRate::Real,frontRowPull::Real) where {T_Float}
        new{T_Float}(T_Float(noise),T_Float(moshRate),T_Float(frontRowPull))
    end #function
end #struct
MusicParams(noise, moshRate, frontRowPull) = MusicParams{typeof(noise)}(noise, moshRate, frontRowPull)