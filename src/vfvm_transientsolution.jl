"""
$(TYPEDEF)

Abstract type for transient solution
"""
abstract type AbstractTransientSolution{T,N,A,B}  <: AbstractDiffEqArray{T, N, A} end

"""
$(TYPEDEF)

Transient solution structure

## Fields
$(TYPEDFIELDS)
    

## Interface

Object of this type adhere to the `AbstractDiffEqArray`  interface.
For indexing and interpolation, see [https://diffeq.sciml.ai/stable/basics/solution/](https://diffeq.sciml.ai/stable/basics/solution/).

In particular, a TransientSolution `sol` can be accessed as follows:
- `sol[it]` contains the solution for timestep `i`
- `sol[ispec,:,it]` contains the solution for component `ispec` at timestep `i`
- `sol(t)` returns a (linearly) interpolated solution value for `t`.
- `sol.t[it]` is the corresponding time
- `sol[ispec,ix,it]` refers to solution of component `ispec` at node `ix` at moment `it`


"""
mutable struct TransientSolution{T,N,A,B} <: AbstractTransientSolution{T,N,A,B}
    """
    Vector of solutions
    """
    u::A
    """
    Vector of times
    """
    t::B
end

TransientSolution(vec::AbstractVector{T}, ts, ::NTuple{N}) where {T, N} = TransientSolution{eltype(T), N, typeof(vec), typeof(ts)}(vec, ts)

TransientSolution(vec::AbstractVector,ts::AbstractVector) = TransientSolution(vec, ts, (size(vec[1])..., length(vec)))

Base.append!(s::TransientSolution,t::Real,sol::AbstractArray)=push!(s.t,t), push!(s.u,copy(sol))

(sol::TransientSolution)(t) = _interpolate(sol,t)

function _interpolate(sol,t)
    idx=searchsortedfirst(sol.t,t)
    if idx==1 || idx >length(sol)
        return nothing
    end
    if t==sol.t[idx]
        return sol[idx]
    else
        retval=similar(sol[idx])
        dt=sol.t[idx]-sol.t[idx-1]
        a=(t-sol.t[idx-1])/dt
        b=(sol.t[idx]-t)/dt
        retval.=sol[idx-1].*a + sol[idx].*b
    end
end


mutable struct VectorOfDiskArrays{T} <: AbstractVector{T}
    fname::String
    file::Union{JLD2.JLDFile,Nothing}
    n::Int64
end

function Base.push!(v::VectorOfDiskArrays, obj)
    v.n+=1
    if isnothing(v.file)
        jldopen(v.fname,"a+") do file
            file[string(v.n)]=obj
        end
    else
        v.file[string(v.n)]=obj
    end
end

Base.size(v::VectorOfDiskArrays)=(v.n,)
Base.length(v::VectorOfDiskArrays)=v.n
Base.eltype(v::VectorOfDiskArrays{T}) where T = T
function Base.getindex(v::VectorOfDiskArrays,i)
    if isnothing(v.file)
        jldopen(v.fname,"r") do file
            file[string(i)]
        end
    else
        v.file[string(i)]
    end
end

_tempname()=Base.VERSION<v"1.4" ? tempname() : tempname(pwd())*".jld2" 

"""
````
VectorOfDiskArrays(firstobj:AbstractArray;
                   keep_open=true,
                   fname= fname=tempname(pwd())*".jld2")
````
Constructor of vector of arrays stored on disk (via JLD2).

- `keep_open`: if true, disk file is not closed during the existence of the object
- `fname`: file name for the disk file


The disk file is automatically removed if the object is garbage collected.
"""
function VectorOfDiskArrays(obj::AbstractArray{T}; keep_open=true, fname=_tempname()) where T
    file=jldopen(fname,"a+")
    file[string(1)]=obj
    if !keep_open
        close(file)
        file=nothing
    end
    v=VectorOfDiskArrays{T}(fname,file,1)
    finalizer(v-> (isnothing(v.file) ? nothing : close(v.file); rm(v.fname,force=true)),v)
    v
end

Base.append!(s::TransientSolution{T,N,VectorOfDiskArrays{T},B}, t::Real,sol::AbstractArray) where {T,N,B} = push!(s.t,t), push!(s.u,sol)

"""
````
TransientSolution(t0,inival;
                  in_memory=false,
                  keep_open=true,
                  fname=tempname(pwd())*".jld2"
````
Constructor of transient solution with initial value and inital time.

- `in_memory`: if true, data are kept in main memory, otherwise on disk (via JLD2)
- `keep_open`: if true, disk file is not closed during the existence of the object
- `fname`: file name for the disk file
"""
function TransientSolution(t0::Number,
                           inival::AbstractArray{T};
                           in_memory=false,
                           keep_open=true,
                           fname=_tempname()) where T
    if !in_memory && !isa(inival,SparseSolutionArray)
        TransientSolution(VectorOfDiskArrays(inival,keep_open=keep_open,fname=fname),[t0])
    else
        TransientSolution([inival], [t0])
    end
end
