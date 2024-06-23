using GLMakie
include("EDMD.jl")

function openDoors!(rs::Vector{SVector{2,Float64}}, θs::Vector{Float64}, params::CrowdParams)
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
                θs[pIdx] = 2π*rand()
            end #if
        end #for yIdx
    end #for xIdx
    return nothing
end #function

function v(r::SVector{2,Float64},θ::Float64,params::CrowdParams)
    modV = params.v₀ * exp(-params.boundaryScale/r[2])
    return modV * SVector{2,Float64}(cos(θ),sin(θ))
end #function

function noiseScaleFactor(r::SVector{2,Float64},params::CrowdParams)
    return exp(-r[2]/params.boundaryScale_noise)
end #function

function U_r(r::SVector{2,Float64},params::CrowdParams)
    #Get mosh pit potential for r 
    return SVector{2,Float64}(0.0,0.0)
end #function

function U_θ(r::SVector{2,Float64}, θ::Float64, params::CrowdParams)
    #Get mosh pit potential for θ 
    return 0.0
end #function

function performBM!(trial_Δrs::Vector{SVector{2,Float64}}, θs::Vector{Float64}, rs::Vector{SVector{2,Float64}}, params::CrowdParams, musicParams::MusicParams)
    for pIdx = eachindex(trial_Δrs)
        #TO-DO:: Have added pull to front, but need to add moshpits
        #Maybe have v also change based on moshpits?
        trial_Δrs[pIdx] = params.dt*v(rs[pIdx],θs[pIdx],params) + U_r(rs[pIdx],params) + sqrt(2*noiseScaleFactor(rs[pIdx],params)*musicParams.noise*params.dt)*randn(2)
        θs[pIdx] += params.dt*musicParams.frontRowPull*sin((-0.5*π)-θs[pIdx]) + U_θ(rs[pIdx],θs[pIdx],params) + sqrt(2*noiseScaleFactor(rs[pIdx],params)*musicParams.angleNoise*params.dt)*randn()
    end #for pIdx
    return nothing
end #function

#Method when the figure has been initialised
function startConcert(figAndBoxes::Tuple{Figure,Textbox,Textbox,Textbox,Textbox,Textbox,Axis}, params::CrowdParams, musicParams::MusicParams, maxFps::Int)

    rs = Vector{SVector{2,Float64}}(undef,params.N)
    θs = Vector{Float64}(undef,params.N)
    openDoors!(rs, θs, params)

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
    tb_angleNoise = figAndBoxes[6]
    ax = figAndBoxes[7]
    display(fig)

    #Define "person" marker for the plot
    personMarker = BezierPath([
       MoveTo(Point(1,0)),
       EllipticalArc(Point(0,0),1,1,0,0,2pi),
       MoveTo(Point(0.505,0.4)),
       EllipticalArc(Point(0.5,0.4),0.22,0.22,0,0,-2pi),
       MoveTo(Point(0.505,-0.4)),
       EllipticalArc(Point(0.5,-0.4),0.22,0.22,0,0,-2pi),
       ClosePath(),
    ])
    #Set colours for the plot
    peopleColors = rand(params.N)
    peopleColorMap = :darkrainbow
    
    tb_noise.stored_string = string(musicParams.noise)
    tb_moshRate.stored_string = string(musicParams.moshRate)
    tb_frontRowPull.stored_string = string(musicParams.frontRowPull)
    tb_angleNoise.stored_string = string(musicParams.angleNoise)
    
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

    #Update angleNoise if value is changed
    on(tb_angleNoise.stored_string) do s
        musicParams[].angleNoise = parse(Float64,s)
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
        GLMakie.scatter!(ax,[rs[i][1] for i in eachindex(rs)], [rs[i][2] for i in eachindex(rs)],markersize = markerSize[], color=peopleColors, colormap = peopleColorMap, marker = personMarker, rotation = θs)
        performBM!(trial_Δrs, θs, rs, params, musicParams[])
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
    Label(fig[2,2][1,5],"angleNoise:")
    
    tb_noise = Textbox(fig[2,2][2,1],validator=Float64,placeholder=string(musicParams.noise),tellwidth=false)
    tb_moshRate = Textbox(fig[2,2][2,2],validator=Float64,placeholder=string(musicParams.moshRate),tellwidth=false)
    tb_frontRowPull = Textbox(fig[2,2][2,3],validator=Float64,placeholder=string(musicParams.frontRowPull),tellwidth=false)
    tb_markersize = Textbox(fig[2,2][2,4],validator=Float64,placeholder=string(16),tellwidth=false)
    tb_angleNoise = Textbox(fig[2,2][2,5],validator=Float64,placeholder=string(musicParams.angleNoise),tellwidth=false)
    
    tb_noise.stored_string = string(musicParams.noise)
    tb_moshRate.stored_string = string(musicParams.moshRate)
    tb_frontRowPull.stored_string = string(musicParams.frontRowPull)
    tb_markersize.stored_string = string(16) #Bodge for now
    tb_angleNoise.stored_string = string(musicParams.angleNoise)

    ax = Axis(fig[1,1:3],limits=(0,params.Lx,0,params.Ly),aspect=AxisAspect(1))

    figAndBoxes = (fig,tb_noise,tb_moshRate,tb_frontRowPull,tb_markersize,tb_angleNoise,ax)

    startConcert(figAndBoxes,params,musicParams,maxFps)

    return figAndBoxes
end #function
