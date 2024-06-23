include("EDMD.jl")

function openDoors!(rs::Vector{SVector{2,Float64}}, params::CrowdParams)
    packingDiameter = 2*maximum(params.Rs) + 0.2
    maxX = floor(Int,params.Lx / packingDiameter)
    maxY = floor(Int,params.Ly / packingDiameter)
    if length(rs) > maxX*maxY
        error("People do not fit in the box")
    end #if
    pIdx = 0
    donePlacing = false
    for xIdx = 1:maxX
        donePlacing && break
        for yIdx = 1:maxY
            pIdx += 1
            if pIdx > length(rs)
                donePlacing = true
                break
            else
                rs[pIdx] = [(xIdx-0.5)*packingDiameter,(yIdx-0.5)*packingDiameter]
            end #if
        end #for yIdx
    end #for xIdx
    return nothing
end #function

function getTrialDisplacements!(trial_Δrs::Vector{SVector{2,Float64}}, rs::Vector{SVector{2,Float64}}, params::CrowdParams, musicParams::MusicParams)
    for pIdx = eachindex(trial_Δrs)
        #TO-DO:: Have added pull to front, but need to add moshpits
        trial_Δrs[pIdx] = sqrt(2*musicParams.noise*params.dt)*randn(2) - (musicParams.frontRowPull*params.dt*(params.Ly - rs[pIdx][2]))*[0,1]
    end #for pIdx
    return nothing
end #function

function startConcert(params::CrowdParams, musicParams::MusicParams)

    rs = Vector{SVector{2,Float64}}(undef,params.N)
    openDoors!(rs,params)

    vs = similar(rs)
    innerδts = Vector{Float64}(undef,params.N)
    n_cols = zeros(Int64,params.N)
    trial_Δrs = similar(rs) 
    eventTree = EventTree([0],params.dt)
    nextEventTimes = zeros(params.N)
    cellList = CellList(rs, params)
    nghbrLists = NeighbourLists(rs, cellList, params)

    anim = @animate for timestep = 1:params.maxSteps
        scatter([rs[i][1] for i in eachindex(rs)], [rs[i][2] for i in eachindex(rs)],xlimits = [0,params.Lx],ylimits = [0,params.Ly], markersize = 20, aspect_ratio=:equal)
        getTrialDisplacements!(trial_Δrs, rs, params, musicParams)
        #TO-DO:: CURRENTLY ONLY INCLUDING FRONT/MOSHPIT FORCES IN THE BROWNIAN MOTION stepsBeforeFirstOPSave
        #--> COULD ADD THEM TO THE EDMD BIT
        performEDMD!(rs, vs, innerδts, n_cols, trial_Δrs, eventTree, nextEventTimes, cellList, nghbrLists, params, musicParams)


    end #for timestep

    return anim
end #function
