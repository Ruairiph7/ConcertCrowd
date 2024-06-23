using GLMakie
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

#Method when the figure has been initialised
function startConcert(figAndBoxes::Tuple{Figure,Textbox,Textbox,Textbox,Textbox,Axis}, params::CrowdParams, musicParams::MusicParams, maxFps::Int)

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

    #Initialise figure
    fig = figAndBoxes[1]
    tb_noise = figAndBoxes[2]
    tb_moshRate = figAndBoxes[3]
    tb_frontRowPull = figAndBoxes[4]
    tb_markersize = figAndBoxes[5]
    ax = figAndBoxes[6]
    display(fig)
    
    tb_noise.stored_string = string(musicParams.noise)
    tb_moshRate.stored_string = string(musicParams.moshRate)
    tb_frontRowPull.stored_string = string(musicParams.frontRowPull)
    
    running = Observable(true)
    musicParams = Observable(musicParams)
    markerSize = Observable(parse(Float64,figAndBoxes[5].stored_string[])) #Bodge for now

    #Update noise if value is changed
    on(tb_noise.stored_string) do s
        musicParams[].noise = parse(Float64,s)
    end #on

    #Update moshRate if value is changed
    on(tb_moshRate.stored_string) do s
        musicParams[].moshRate = parse(Float64,s)
    end #on

    #Update frontRowPull if value is changed
    on(tb_frontRowPull.stored_string) do s
        musicParams[].frontRowPull = parse(Float64,s)
    end #on

    #Update markersize if value is changed
    on(tb_markersize.stored_string) do s
        markerSize[] = parse(Float64,s)
    end #on

    #Quit if user presses q
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.q
                running[] = false
            end #if
        end #if
    end #on

    while running[] == true
        GLMakie.scatter!(ax,[rs[i][1] for i in eachindex(rs)], [rs[i][2] for i in eachindex(rs)],markersize = markerSize[],color=:blue)
        getTrialDisplacements!(trial_Δrs, rs, params, musicParams[])
        #TO-DO:: CURRENTLY ONLY INCLUDING FRONT/MOSHPIT FORCES IN THE BROWNIAN MOTION stepsBeforeFirstOPSave
        #--> COULD ADD THEM TO THE EDMD BIT
        performEDMD!(rs, vs, innerδts, n_cols, trial_Δrs, eventTree, nextEventTimes, cellList, nghbrLists, params, musicParams[])
        sleep(1/maxFps)
        empty!(ax)
    end #for timestep

    return nothing
end #function

#Method when we need to initialise the figure
function startConcert(params::CrowdParams, musicParams::MusicParams, maxFps::Int)
    fig = Figure()
    #Make text boxes to change musicParams, markersize and N
    #Bodge for now - manually change markersize until it visually matches
    Label(fig[2,2][1,1],"noise:")
    Label(fig[2,2][1,2],"moshRate:")
    Label(fig[2,2][1,3],"frontRowPull:")
    Label(fig[2,2][1,4],"markersize:")
    
    tb_noise = Textbox(fig[2,2][2,1],validator=Float64,placeholder=string(musicParams.noise),tellwidth=false)
    tb_moshRate = Textbox(fig[2,2][2,2],validator=Float64,placeholder=string(musicParams.moshRate),tellwidth=false)
    tb_frontRowPull = Textbox(fig[2,2][2,3],validator=Float64,placeholder=string(musicParams.frontRowPull),tellwidth=false)
    tb_markersize = Textbox(fig[2,2][2,4],validator=Float64,placeholder=string(45),tellwidth=false)
    
    tb_noise.stored_string = string(musicParams.noise)
    tb_moshRate.stored_string = string(musicParams.moshRate)
    tb_frontRowPull.stored_string = string(musicParams.frontRowPull)
    tb_markersize.stored_string = string(45) #Bodge for now

    ax = Axis(fig[1,1:3],limits=(0,params.Lx,0,params.Ly),aspect=AxisAspect(1))

    figAndBoxes = (fig,tb_noise,tb_moshRate,tb_frontRowPull,tb_markersize,ax)

    startConcert(figAndBoxes,params,musicParams,maxFps)

    return figAndBoxes
end #function
