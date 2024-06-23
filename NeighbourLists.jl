using LinearAlgebra
include("ParamStructDefinitions.jl")

##############################################################
# Cell Lists -------------------------------------------------
##############################################################

#Define CellList struct:
# -- Auto-construct lists of cell neighbours when called
##############################################################
mutable struct CellList{T<:Real}
    list::Vector{Int64} #List of length N giving the cell each person is in
    numBoxesX::Int64 #Number of cells in x direction
    numBoxesY::Int64 #Number of cells in y direction
    cBoxSize::T #Size of a cell box
    cellNeighbours::Vector{SVector} #Vector containing list of neighbour sIDs for each cell

    function CellList{T}(list,numBoxesX,numBoxesY,cBoxSize) where {T<:Real}
        #Construct the cell neighbours:
        cellNeighbours = Any[Int[sID] for sID = 1:numBoxesX*numBoxesY]
        for xIdx = 1:numBoxesX
            for yIdx = 1:numBoxesY
                sID = xIdx + numBoxesX*(yIdx - 1)
                vID = [xIdx, yIdx]
                #Box above
                yIdx != 1 && push!(cellNeighbours[sID],vectorToScalarBoxID(vID + [0,-1],numBoxesX,numBoxesY))
                #Box top right
                yIdx != 1 && xIdx != numBoxesX && push!(cellNeighbours[sID],vectorToScalarBoxID(vID + [1,-1],numBoxesX,numBoxesY))
                #Box right
                xIdx != numBoxesX && push!(cellNeighbours[sID],vectorToScalarBoxID(vID+[1,0],numBoxesX,numBoxesY))
                #Box bottom right
                xIdx != numBoxesX && yIdx != numBoxesY && push!(cellNeighbours[sID],vectorToScalarBoxID(vID+[1,1],numBoxesX,numBoxesY))
                #Box below
                yIdx != numBoxesY && push!(cellNeighbours[sID],vectorToScalarBoxID(vID+[0,1],numBoxesX,numBoxesY))
                #Box bottom left
                yIdx != numBoxesY && xIdx != 1 && push!(cellNeighbours[sID],vectorToScalarBoxID(vID+[-1,1],numBoxesX,numBoxesY))
                #Box left
                xIdx != 1 && push!(cellNeighbours[sID],vectorToScalarBoxID(vID+[-1,0],numBoxesX,numBoxesY))
                #Box top left
                xIdx != 1 && yIdx != 1 && push!(cellNeighbours[sID],vectorToScalarBoxID(vID+[-1,-1],numBoxesX,numBoxesY))
                cellNeighbours[sID] = SVector{length(cellNeighbours[sID]),Int64}(cellNeighbours[sID])
            end #for yIdx
        end #for xIdx
        cellNeighbours = Vector{SVector}(cellNeighbours)
        new{T}(list,numBoxesX,numBoxesY,cBoxSize,cellNeighbours)
    end #function
end #struct
CellList(list,numBoxesX,numBoxesY,cBoxSize) = CellList{typeof(cBoxSize)}(list,numBoxesX,numBoxesY,cBoxSize)
CellList(list,params::CrowdParams) = CellList(list,Int(params.Lx/params.cBoxSize),Int(params.Ly/params.cBoxSize),params.cBoxSize)
function CellList(rs::Vector{<:SVector{2,<:Real}}, params::CrowdParams)
    cellList = CellList(zeros(Int64,params.N),params)
    fillList!(cellList, rs)
    return cellList
end #function

#Functions for converting between scalar and vector cell indices 
########################################################
function scalarToVectorBoxID(sID::Int,numBoxesX::Int,numBoxesY::Int)
    vID = Vector{Int64}(undef,2)
    #x component:
    if sID%numBoxesX != 0
        vID[1] = sID%numBoxesX
    else
        vID[1] = numBoxesX
    end #if
    #y component:
    for j = 1:numBoxesY
        if (j-1)*numBoxesX < sID <= j*numBoxesX
            vID[2] = j
            break
        end #if
    end #for j
    return vID
end #function

function vectorToScalarBoxID(vID::Vector{<:Int},numBoxesX::Int,numBoxesY::Int)
    #Apply PBCs (for neighbours outside the grid)
    if vID[2] == numBoxesY + 1 #If above box, correct down
        true_vID_y = 1
    elseif vID[2] == 0 #If below box, correct up
        true_vID_y = numBoxesY
    else
        true_vID_y = vID[2]
    end #if
    if vID[1] == numBoxesX + 1 #If to right, correct left
        true_vID_x = 1
    elseif vID[1] == 0 #If to left, correct right
        true_vID_x = numBoxesX 
    else
        true_vID_x = vID[1]
    end #if
    return true_vID_x + numBoxesX *(true_vID_y - 1)
end #function

#Functions that operate on objects of type CellList
##############################################################
function findCell(cellList::CellList, rs::Vector{<:SVector{2,<:Real}}, i::Int)
    return ceil(Int,rs[i][1]/cellList.cBoxSize) + cellList.numBoxesX*floor(Int,rs[i][2]/cellList.cBoxSize)
end #function

function findCell!(cellList::CellList, rs::Vector{<:SVector{2,<:Real}}, i::Int)
    cellList.list[i] = findCell(cellList,rs,i)
end #function

function fillList!(cellList::CellList, rs::Vector{<:SVector{2,<:Real}})
    for i = eachindex(cellList.list)
        findCell!(cellList, rs, i)
    end
    return nothing
end #function


##############################################################
# Neighbour Lists --------------------------------------------
##############################################################

#Define NeighbourList and NeighbourLists structs
##############################################################
mutable struct NeighbourList{T<:AbstractVector{<:Int64}}
    list::T #List of neighbours for the person
    i::Int64 #Index of the person
end #struct
NeighbourList(list, i) = NeighbourList{typeof(list)}(list, i)
NeighbourList(maxNeighbours::Int, i) = NeighbourList(zeros(Int64,maxNeighbours), i) 
NeighbourList(maxNeighbours::Int, i, SVector_flag::Any) = NeighbourList(zeros(SVector{maxNeighbours,Int64}), i)

mutable struct NeighbourLists{T<:Real}
    lists::Vector{NeighbourList}
    r_ns::Vector{<:SVector{2,T}}
    alpha::T
    maxNeighbours::Int64
end #struct
NeighbourLists(lists, r_ns, alpha, maxNeighbours) = NeighbourLists{typeof(alpha)}(lists, r_ns, alpha, maxNeighbours)
NeighbourLists(params::CrowdParams) = NeighbourLists(NeighbourList[NeighbourList(params.maxNeighbours,i) for i = 1:params.N],zeros(SVector{2,typeof(params.alpha)},params.N),params.alpha, params.maxNeighbours)
function NeighbourLists(rs::Vector{<:SVector{2,<:Real}}, cellList::CellList, params::CrowdParams)
    neighbourLists = NeighbourLists(params)
    for i = eachindex(neighbourLists.lists)
        fillList!(neighbourLists.lists[i],rs,cellList,params)
        neighbourLists.r_ns[i] = rs[i]
    end #for i
    return neighbourLists
end #function

#Functions that operate on objects of type NeighbourList
##############################################################
function addToList!(neighbourList::NeighbourList, j::Int64)
    idx = findnext(iszero,neighbourList.list,1)
    neighbourList.list[idx] = j
end #function

function removeFromList!(neighbourList::NeighbourList, j::Int64)
    idx = findall(==(j),neighbourList.list)
    length(idx) != 1 && error("Either j was not in list, or was in it more than once")
    neighbourList.list[idx] .= 0
end #function

function fillList!(neighbourList::NeighbourList, r_ns::Vector{<:SVector{2,<:Real}}, cellList::CellList, params::CrowdParams)
    i = neighbourList.i
    sID = cellList.list[i]
    cellNeighbours = cellList.cellNeighbours[sID]
    for cellNeighbour = cellNeighbours
        pIdxs = findall(==(cellNeighbour),cellList.list)
        for j = pIdxs
            j == i && continue
            if norm(r_ns[j]-r_ns[i]) < (1+params.alpha)*(params.Rs[i]+params.Rs[j])
                addToList!(neighbourList,j)
            end #if
        end #for pIdx
    end #for cellNeighbour
end #function

function emptyList!(neighbourList::NeighbourList)
    neighbourList.list = zeros(Int64,length(neighbourList.list))
end #function


