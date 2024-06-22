using AbstractTrees

##############################################################
#Binary Search Tree for the event calendar--------------------
##############################################################

struct EventData{T1<:Int,T2<:Real}
    details::AbstractVector{T1}
    δt::T2
end #struct
EventData(details,δt) = EventData{eltype(details),typeof(δt)}(details,δt)

maybeNothing{T} = Union{T,Nothing}
mutable struct EventTree{T1<:Int,T2<:Real}
    eventData::EventData{T1,T2}
    leftBranch::maybeNothing{EventTree}
    rightBranch::maybeNothing{EventTree}

    function EventTree{T1,T2}(eventData, l=nothing, r=nothing) where {T1<:Int, T2<:Real}
        new{T1,T2}(eventData, l, r)
    end #function
end #struct
EventTree(eventData::EventData) = EventTree{eltype(eventData.details),typeof(eventData.δt)}(eventData)
EventTree(details::AbstractVector{<:Int}, δt::Real) = EventTree(EventData(details,δt))
EventTree() = EventTree([0],Inf)

function isLessThan(event1::EventData, event2::EventData)
    if event1.δt < event2.δt
        return true
    elseif event1.δt > event2.δt
        return false
    elseif event1.δt == event2.δt
        if event1.details[1] < event2.details[1]
            return true
        elseif event1.details[1] > event2.details[1]
            return false
        elseif event1.details[1] == event2.details[1]
            return false
        end #if
    end #if
end #function

function isMoreThan(event1::EventData, event2::EventData)
    if event1.δt > event2.δt
        return true
    elseif event1.δt < event2.δt
        return false
    elseif event1.δt == event2.δt
        if event1.details[1] > event2.details[1]
            return true
        elseif event1.details[1] < event2.details[1]
            return false
        elseif event1.details[1] == event2.details[1]
            return false
        end #if
    end #if
end #function

function addToTree!(node::maybeNothing{EventTree}, eventData::EventData)
    if isLessThan(eventData, node.eventData)
        if isnothing(node.leftBranch)
            node.leftBranch = EventTree(eventData)
        else
            addToTree!(node.leftBranch, eventData)
        end #if
    elseif isMoreThan(eventData, node.eventData)
        if isnothing(node.rightBranch)
            node.rightBranch = EventTree(eventData)
        else 
            addToTree!(node.rightBranch, eventData)
        end #if
    else
        error("Adding an 'identical' event --> should not be possible")
    end #if
    return nothing
end #function
addToTree!(node::maybeNothing{EventTree}, details::AbstractVector{<:Int}, δt::Real) = addToTree!(node, EventData(details,δt))


function getTreeMin(node::EventTree)
    while !isnothing(node.leftBranch)
        node = node.leftBranch
    end #while
    return node.eventData
end #function


function deleteNode(node::maybeNothing{EventTree}, dataToBeDeleted::EventData)
    #Return nothing if node given is nothing 
    #(Needed for the case that dataToBeDeleted is not in the tree)
    if isnothing(node)
        return nothing
    end #if
    #Move down the tree in the right direction to find δtToBeDeleted
    if isLessThan(node.eventData, dataToBeDeleted)
        node.rightBranch = deleteNode(node.rightBranch, dataToBeDeleted)
    elseif isMoreThan(node.eventData, dataToBeDeleted)
        node.leftBranch = deleteNode(node.leftBranch, dataToBeDeleted)
    #We have found δtToBeDeleted
    else
        #Node has no children
        if isnothing(node.leftBranch) && isnothing(node.rightBranch)
            return nothing
        #Node has only a right child
        elseif isnothing(node.leftBranch)
            return node.rightBranch
        #Node has only a left child
        elseif isnothing(node.rightBranch)
            return node.leftBranch
        #Node has both left and right children
        else
            rightMin = getTreeMin(node.rightBranch)
            node = deleteNode(node, rightMin)
            node.eventData = rightMin
            return node
        end #if
    end #if
    return node
end #function

function deleteNode!(node::maybeNothing{EventTree}, dataToBeDeleted::EventData)
    node = deleteNode(node,dataToBeDeleted)
    return nothing
end #function
deleteNode!(node::maybeNothing{EventTree}, details::AbstractVector{<:Int}, δt::Real) = deleteNode!(node,EventData(details,δt))
deleteNode!(node::maybeNothing{EventTree}, idx::Int, δt::Real) = deleteNode!(node,[idx],δt)

##############################################################
#Make use of the AbstractTrees.jl interface 
##############################################################

function AbstractTrees.children(node::EventTree)
    if isnothing(node.leftBranch) && isnothing(node.rightBranch)
        return ()
    elseif isnothing(node.leftBranch) && !isnothing(node.rightBranch)
        return (node.rightBranch,)
    elseif !isnothing(node.leftBranch) && isnothing(node.rightBranch)
        return (node.leftBranch,)
    else
        return (node.leftBranch, node.rightBranch)
    end #if
end #function

AbstractTrees.nodevalue(n::EventTree) = (n.eventData.details,n.eventData.δt)

AbstractTrees.NodeType(::Type{<:EventTree{T1,T2}}) where {T1<:Int,T2<:Real} = HasNodeType()
AbstractTrees.nodetype(::Type{<:EventTree{T1,T2}}) where {T1<:Int,T2<:Real} = EventTree{T1,T2}