include("EDMD.jl")

function openDoors!(rs::Vector{SVector{2,Float64}}, params::CrowdParams)
    for pIdx = eachindex(rs)
        rs[pIdx] = [params.Rs[pIdx] + (pIdx-1) * params.Rs[pIdx], params.Rs[pIdx]]
        # rs[pIdx] = [params.Lx*rand(),params.Ly*rand()]
    end #for pIdx
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

    anim = @animate for timestep = 1:params.maxSteps
        scatter([rs[i][1] for i in eachindex(rs)], [rs[i][2] for i in eachindex(rs)],xlimits = [0,params.Lx],ylimits = [0,params.Ly], markersize = 20, aspect_ratio=:equal)
        getTrialDisplacements!(trial_Δrs, rs, params, musicParams)
        #TO-DO:: CURRENTLY ONLY INCLUDING FRONT/MOSHPIT FORCES IN THE BROWNIAN MOTION stepsBeforeFirstOPSave
        #--> COULD ADD THEM TO THE EDMD BIT
        performEDMD!(rs, vs, innerδts, n_cols, trial_Δrs, eventTree, nextEventTimes, params, musicParams)


    end #for timestep

    return anim
end #function
