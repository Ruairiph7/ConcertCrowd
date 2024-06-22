using AbstractTrees

##############################################################
#Binary Search Tree for the event calendar--------------------
##############################################################

maybeNothing{T} = Union{T,Nothing}
mutable struct EventTree{T1<:Int,T2<:Real}
    people::AbstractVector{T1}
    δt::T2
    leftBranch::maybeNothing{EventTree}
    rightBranch::maybeNothing{EventTree}

    function EventTree{T1,T2}(people, δt, l=nothing, r=nothing) where {T1<:Int, T2<:Real}
        new{T1,T2}(people, δt, l, r)
    end #function
end #struct
EventTree(people, δt) = EventTree{eltype(people),typeof(δt)}(people, δt)
EventTree() = EventTree([0],Inf)

function addToTree!(node::maybeNothing{EventTree}, people::AbstractVector{<:Int}, δt::Real)
    if δt < node.δt
        if isnothing(node.leftBranch)
            node.leftBranch = EventTree(people, δt)
        else
            addToTree!(node.leftBranch, people, δt)
        end #if
    elseif δt >= node.δt
        if isnothing(node.rightBranch)
            node.rightBranch = EventTree(people, δt)
        else 
            addToTree!(node.rightBranch, people, δt)
        end #if
    end #if
    return nothing
end #function


function getTreeMin(node::EventTree)
    while !isnothing(node.leftBranch)
        node = node.leftBranch
    end #while
    return node.people, node.δt
end #function


function deleteNode(node::maybeNothing{EventTree}, δtToBeDeleted::Real)
    #Return nothing if node given is nothing 
    #(Needed for the case that δtToBeDeleted is not in the tree)
    if isnothing(node)
        return nothing
    end #if
    #Move down the tree in the right direction to find δtToBeDeleted
    if node.δt < δtToBeDeleted
        node.rightBranch = deleteNode(node.rightBranch, δtToBeDeleted)
    elseif node.δt > δtToBeDeleted
        node.leftBranch = deleteNode(node.leftBranch, δtToBeDeleted)
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
            node = deleteNode(node, rightMin[2])
            (node.people, node.δt) = rightMin
            return node
        end #if
    end #if
    return node
end #function

function deleteNode!(node::maybeNothing{EventTree}, δtToBeDeleted::Real)
    node = deleteNode(node,δtToBeDeleted)
    return nothing
end #function


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

AbstractTrees.nodevalue(n::EventTree) = (n.δt)

AbstractTrees.NodeType(::Type{<:EventTree{T1,T2}}) where {T1<:Int,T2<:Real} = HasNodeType()
AbstractTrees.nodetype(::Type{<:EventTree{T1,T2}}) where {T1<:Int,T2<:Real} = EventTree{T1,T2}