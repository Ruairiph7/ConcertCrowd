using LinearAlgebra
include("EventTrees.jl")
include("ParamStructDefinitions.jl")

##############################################################
#Event-Driven Molecular dynamics within each dt---------------
##############################################################

function updateParticleToTime!(rs::Vector{SVector{2,Float64}}, vs::Vector{SVector{2,Float64}}, innerδts::Vector{Float64}, desiredδt::Float64, i::Int)
    #NOTE: As we only ever deal with the next collision, we will never be calling
    #this function to update a particle that will collide in the given time
    rs[i] += vs[i]*(desiredδt - innerδts[i])
    innerδts[i] = desiredδt
    return nothing
end #function


function getTimeToWallCollisions(r::SVector{2,Float64}, v::SVector{2,Float64}, R::Float64, params::CrowdParams)
    #Compute time till collision with a wall in the x-direction
    if v[1] > 0
        δtx::Float64 = (params.Lx-R-r[1])/v[1]
    else 
        δtx = (R-r[1])/v[1]
    end #if
    #Compute time till collision with a wall in the y-direction
    if v[2] > 0
        δty::Float64 = (params.Ly-R-r[2])/v[2]
    else 
        δty = (R-r[2])/v[2]
    end #if
    return δtx, δty
end #function

function getTimeToNeighbourUpdate()
    return Inf
end #function

function getTimeToParticleCollision(rᵢ::SVector{2,Float64}, rⱼ::SVector{2,Float64}, vᵢ::SVector{2,Float64}, vⱼ::SVector{2,Float64}, innerδtᵢ::Float64, innerδtⱼ::Float64, Rᵢ::Float64, Rⱼ::Float64)
    #Roll j back/forwards in time to match innerδtᵢ
    r_ij::SVector{2,Float64} = (rⱼ + vⱼ*(innerδtᵢ - innerδtⱼ)) - rᵢ
    v_ij::SVector{2,Float64} = vⱼ - vᵢ
    b::Float64 = r_ij ⋅ v_ij
    sqrtTerm::Float64 = b^2 - (v_ij⋅v_ij)*((r_ij⋅r_ij) - (Rᵢ + Rⱼ)^2) 

    if b < 0 && sqrtTerm >= 0
        δt = ( -b - sqrt(sqrtTerm) )/((v_ij⋅v_ij))
    else 
        δt = Inf
    end #if
    return δt
end #function


function updateTreeForPerson!(eventTree::EventTree, nextEventTimes::Vector{Float64}, i::Int, rs::Vector{SVector{2,Float64}}, vs::Vector{SVector{2,Float64}}, innerδts::Vector{Float64}, n_cols::Vector{Int64}, params::CrowdParams)
    #Store:
        #x wall collisions as node.people = [i,0]
        #y wall collisions as node.people = [i,-1]
        #neighbour updates as node.people = [i,-2]
        #particle collisions as node.people = [i,j,n_col_j]

    (xWallδt, yWallδt) = getTimeToWallCollisions(rs[i],vs[i],params.Rs[i],params)
    neighbourListδt = getTimeToNeighbourUpdate()

    (minδtSoFar,typeIndicator) = findmin((xWallδt,yWallδt,neighbourListδt))
    if typeIndicator == 1
        updateSoFar = [i, 0]
    elseif typeIndicator == 2
        updateSoFar = [i,-1]
    else
        updateSoFar = [i,-2]
    end #if 

    for j = eachindex(rs)
        j == i && continue
        δt_ij = getTimeToParticleCollision(rs[i],rs[j],vs[i],vs[j],innerδts[i],innerδts[j],params.Rs[i],params.Rs[j])
        δt_ij > minδtSoFar && continue
        minδtSoFar = δt_ij
        updateSoFar = [i,j,n_cols[j]]
    end #for j

    #Update nothing if the event occurs after time dt
    if minδtSoFar >= (params.dt - innerδts[i])
        return nothing
    else
        #Add an event to the tree at time (time to event) + (inner time of particle)
        addToTree!(eventTree,updateSoFar,minδtSoFar + innerδts[i])
        nextEventTimes[i] = minδtSoFar + innerδts[i]
        return nothing
    end #if

    return nothing
end #function


function performEvent!(rs::Vector{SVector{2,Float64}}, vs::Vector{SVector{2,Float64}}, innerδts::Vector{Float64}, n_cols::Vector{Int64}, eventData::EventData, eventTree::EventTree, nextEventTimes::Vector{Float64}, params::CrowdParams)

    if length(eventData.details) == 2
        if eventData.details[2] == 0 #It's an x wall collision

            i = eventData.details[1] #Get idx of person involved
            updateParticleToTime!(rs,vs,innerδts,eventData.δt,i)
            vs[i] = SVector(-vs[i][1],vs[i][2])
            n_cols[i] += 1
            updateTreeForPerson!(eventTree, nextEventTimes, i,rs,vs,innerδts,n_cols,params)

        elseif eventData.details[2] == -1 #It's a y wall collision

            i = eventData.details[1] #Get idx of person involved
            updateParticleToTime!(rs,vs,innerδts,eventData.δt,i)
            vs[i] = SVector(vs[i][1],-vs[i][2])
            n_cols[i] += 1
            updateTreeForPerson!(eventTree, nextEventTimes, i,rs,vs,innerδts,n_cols,params)

        elseif eventData.details[2] == -2 #It's a neighbourlist update

        end #if
    elseif length(eventData.details) == 3 #It's a person collision

        i = eventData.details[1]
        j = eventData.details[2]

        if n_cols[j] != eventData.details[3] #j has collided since scheduling the event
            updateTreeForPerson!(eventTree, nextEventTimes, i,rs,vs,innerδts,n_cols,params)
            return nothing
        else
            updateParticleToTime!(rs,vs,innerδts,eventData.δt,i)
            updateParticleToTime!(rs,vs,innerδts,eventData.δt,j)
            
            r_ij = rs[j] - rs[i]
            v_ij = vs[j] - vs[i]
            sigma_ij = params.Rs[i] + params.Rs[j]
            J = 2 * params.masses[i] * params.masses[j] * (v_ij⋅r_ij) * r_ij / ((sigma_ij^2)*(params.masses[i] + params.masses[j]))
            vs[i] += J/params.masses[i]
            vs[j] -= J/params.masses[j]
            n_cols[i] += 1
            n_cols[j] += 1

            updateTreeForPerson!(eventTree, nextEventTimes, i,rs,vs,innerδts,n_cols,params)
            #For j we need to remove their stored event first
            deleteNode!(eventTree, j, nextEventTimes[j])
            updateTreeForPerson!(eventTree, nextEventTimes, j,rs,vs,innerδts,n_cols,params)

        end #if

    else #length(eventData) == 1 ---> its the initial node ([0], dt)
        return nothing
    end #if

    return nothing
end #function


function performEDMD!(rs::Vector{SVector{2,Float64}}, vs::Vector{SVector{2,Float64}}, innerδts::Vector{Float64}, n_cols::Vector{Int64}, trial_Δrs::Vector{SVector{2,Float64}}, eventTree::EventTree, nextEventTimes::Vector{Float64}, params::CrowdParams, musicParams::MusicParams)
    
    #Empty event tree --- MAYBE THIS HAPPENS BY DEFAULT AFTER ALL UPDATES?
    eventTree = EventTree([0],params.dt)

    #Set innerδts, MDSimTime to be zero
    innerδts = zeros(params.N)
    MDSimTime = 0.0

    #Initialise the tree with first updates for each person
    #Also store the fictitious velocities in vs for future ease
    for pIdx = eachindex(rs)
        updateTreeForPerson!(eventTree, nextEventTimes, pIdx, rs, trial_Δrs/params.dt,innerδts,n_cols,params)
        vs[pIdx] = trial_Δrs[pIdx]/params.dt
    end #for pIdx

    while MDSimTime < params.dt
        nextEventData = getTreeMin(eventTree)
        deleteNode!(eventTree, nextEventData)
        performEvent!(rs, vs, innerδts, n_cols, nextEventData, eventTree, nextEventTimes, params)
        MDSimTime = nextEventData.δt
    end #while MDSimTime

    #We have reached the end of the timestep with no further collisions
    #--> update each person to time dt on their current trajectory
    for pIdx = eachindex(rs)
        updateParticleToTime!(rs,vs,innerδts,params.dt,pIdx)
    end #for pIdx

    return nothing
end #function