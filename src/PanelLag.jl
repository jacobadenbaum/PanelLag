module PanelLag

using Tables, DataFrames
import Base: diff

export PanelSet, LaggedField, lag, diff, lag!, diff!

struct PanelSet{T, pid, tid, TF}
    df::T                       # Parent Table
    sp::Vector{Int}             # Sort Permutation
    sp_inv::Vector{Int}         # Inverse Permutation
    Δ::TF
end

function PanelSet(df::T, pid::Symbol, tid::Symbol, sp::Vector{Int}, sp_inv::Vector{Int},Δ::TF) where
    {T, N, TF}

    return PanelSet{T, pid, tid, TF}(df, sp, sp_inv, Δ)
end

tuplefy(x) = (x,)
tuplefy(x::AbstractVector) = Tuple(x)
tuplefy(x::Tuple) = x

struct Keys
    pid::Int
    tid::Int
end

Base.Tuple(key::Keys) = (key.pid, key.tid)
Base.isless(key1::Keys, key2::Keys) = Tuple(key1) < Tuple(key2)

function PanelSet(df::AbstractDataFrame, pid, tid; Δ = 1, safe = false)

    # Don't add new columns to parent
    df = copy(df)

    # Calculate a unique panel id
    nm = pid_name(df)
    gd = groupby(df, vcat(tuplefy(pid)...))
    df[!, nm] = groupindices(gd)

    # Check whether keys are unique
    if !safe
        Ns = combine(gd, N1 = tid => length, N2 = tid => length ∘ unique)
        @assert(Ns[!, :N1] == Ns[!, :N2], "The keys must be unique")
    end

    sp      = sortperm(Keys.(df[!, nm], df[!, tid]))
    sp_inv  = invperm(sp)

    return PanelSet(columntable(df), nm, Symbol(tid), sp, sp_inv, Δ)
end

function pid_name(df)::Symbol
    if :_pid in names(df)
        return :_pid
    else
        i = 1
        while true
            nm = Symbol(:_pid, i)
            if ! ( nm in names(df) )
                return nm
            end
            i += 1
        end
    end
end

pid(ps::PanelSet{T, p, t, TF}) where {T, p, t, TF} = p
tid(ps::PanelSet{T, p, t, TF}) where {T, p, t, TF} = t

pid(::Type{PanelSet{T, p, t, TF}}) where {T, p, t, TF} = p
tid(::Type{PanelSet{T, p, t, TF}}) where {T, p, t, TF} = t

struct LaggedField{TI, T, pid, tid, field, TF} <: AbstractArray{TI, 1}
    ps::PanelSet{T, pid, tid, TF}
    l::Int
end

Base.size(lf::LaggedField) = (length(lf.ps.sp), )
Base.IndexStyle(lf::LaggedField) = IndexLinear()
Base.setindex!(lf::LaggedField, i, val)  = throw(MethodError("You cannot mutate a Lagged View"))

function LaggedField(ps::PanelSet{T,pid, tid, TF}, field, l) where {T, pid, tid, TF}
    TI = Union{eltype(getfield(ps.df, field)), Missing}
    return LaggedField{TI, T, pid, tid, field, TF}(ps, l)
end

field(lf::LaggedField{TI, T, p, t, f, TF}) where {TI, T, p, t, f, TF}   = f
pid(lf::LaggedField{TI, T, p, t, f, TF}) where {TI, T, p, t, f, TF}     = p
tid(lf::LaggedField{TI, T, p, t, f, TF}) where {TI, T, p, t, f, TF}     = t

field(::Type{LaggedField{TI, T, p, t, f, TF}}) where {TI, T, p, t, f, TF} = f
pid(lf::Type{LaggedField{TI, T, p, t, f, TF}}) where {TI, T, p, t, f, TF} = p
tid(lf::Type{LaggedField{TI, T, p, t, f, TF}}) where {TI, T, p, t, f, TF} = t

function Base.getindex(ps::PanelSet, i::Int, col)
    return getfield(ps.df, col)[ps.sp[i]]
end

"""
```
same_panel(ps::PanelSet, i1, i2)
```
Checks whether or not rows `i1` and `i2` are part of the same panel
"""
@generated function same_panel(ps::PanelSet, i1, i2)
    # Check for equality field by field
    pf = Expr(:call, :getfield, :(ps.df), QuoteNode(pid(ps)))
    return quote
        @inbounds if coalesce($pf[i1] != $pf[i2],
                    !(ismissing($pf[i1]) & ismissing($pf[i2])))
            return false
        end
        return true
    end
end


function Base.getindex(lf::LaggedField, i)
    i′ = laggedindex(lf, i)
    if ismissing(i′)
        return missing
    else
        return @inbounds getfield(lf.ps.df, field(lf))[i′]
    end
end

function laggedindex(lf::LaggedField, i)::Union{Missing, Int}

    ps = lf.ps
    N  = length(ps.sp)

    # The time variable
    t  = getfield(ps.df, tid(ps))

    # Calculate our target time
    si      = ps.sp_inv[i]
    @inbounds target  = t[i] - lf.l * ps.Δ

    # Loop until we find the right one, or we're convinced it's missing
    si′ = si
    i′  = i
    while true

        # Check whether the new time matches the original time + l * Δ
        if @inbounds t[i′] == target
            return i′

        # With a positive lag, t[si′] should be greater than the target (since
        # we are lagging back from the current observation)
        elseif @inbounds sign(lf.l) * (t[i′] - target) < 0
            return missing
        end

        # Keep lagging it back (checking that we're not out of bounds)
        @inbounds si′ -= sign(lf.l)
        1 <= si′ <= N  || return missing
        @inbounds i′   = ps.sp[si′]

        # Check that we're still in the same panel
        if !same_panel(ps, i, i′)
            return missing
        end
    end
end

function lag(df, pid, tid, field::Symbol, l = 1; Δ = 1)
    ps = PanelSet(df, pid, tid; Δ = Δ)
    return lag(ps, field, l)
end

function lag(ps::PanelSet, field::Symbol, l = 1)
    lf = LaggedField(ps, field, l)
    return collect(lf)
end


function diff(df, pid, tid, field::Symbol, l = 1; Δ = 1)
    ps = PanelSet(df, pid, tid; Δ = Δ)
    return diff(ps, field, l)
end

function diff(ps::PanelSet, field::Symbol, l = 1)
    return getfield(ps.df, field) - lag(ps, field, l)
end


# Inplace Operation
function lag!(df::AbstractDataFrame, pid, tid, fields, l = 1; Δ = 1)
    ps = PanelSet(df, pid, tid; Δ = Δ)
    for field in vcat(fields)
        df[!, lag(field, l)] = lag(ps, field, l)
    end
end

function lag!(df::AbstractDataFrame, ps::PanelSet, fields, l = 1)
    for field in vcat(fields)
        df[!, lag(field, l)] = lag(ps, field, l)
    end
end

function diff!(df::AbstractDataFrame, pid, tid, fields, l = 1; Δ = 1)
    ps = PanelSet(df, pid, tid; Δ = Δ)
    for field in vcat(fields)
        df[!, diff(field, l)] = lag(ps, field, l)
    end
end

function diff!(df::AbstractDataFrame, ps::PanelSet, fields, l = 1)
    for field in vcat(fields)
        df[!, diff(field, l)] = diff(ps, field, l)
    end
end

################################################################################
#################### Manipulate Symbols for Variable Names #####################
################################################################################
lag(s::Symbol, l::Integer) = begin
    if l > 0
        return Symbol(:L, l, :_, s)
    elseif l < 0
        return Symbol(:F, abs(l), :_, s)
    else
        return s
    end
end

diff(s::Symbol, l::Integer) = begin
    if l > 0
        return Symbol(:D, l, :_, s)
    elseif l < 0
        return Symbol(:FD, abs(l), :_, s)
    else
        return s
    end
end

end
