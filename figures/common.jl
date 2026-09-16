# Shared helpers for the figure scripts: reading the historical ACE/gr (xmgr 4)
# and Grace 5 project files that hold the archived curves, and exporting figures
# as both PNG and PDF (CairoMakie, headless).

using CairoMakie

"""
    XmgrSet

One data set of an ACE/gr or Grace project: the `data` matrix has one row per
point (x, y) or (x, y, dy) depending on `type`; `props` keeps the raw style keys
(`symbol`, `line type`, `color`, ...).
"""
struct XmgrSet
    graph::String
    name::String
    type::String
    data::Matrix{Float64}
    legend::String
    props::Dict{String,String}
end

"""
    XmgrGraph

One graph of an ACE/gr or Grace project, with the archived axis annotations and
its data sets.
"""
mutable struct XmgrGraph
    name::String
    title::String
    xlabel::String
    ylabel::String
    logx::Bool
    logy::Bool
    world::Dict{String,Float64}
    sets::Vector{XmgrSet}
end

XmgrGraph(name::String) = XmgrGraph(name, "", "", "", false, false,
                                    Dict{String,Float64}(), XmgrSet[])

_open_maybe_gz(path::AbstractString) =
    endswith(path, ".gz") ? open(`gzip -dc $path`) : open(path)

function read_lines(path::AbstractString)
    isfile(path) || error("Missing figure data: $path")
    io = _open_maybe_gz(path)
    try
        return readlines(io)
    finally
        close(io)
    end
end

"""
    read_xmgr(path)

Parse an ACE/gr (`@TARGET S0`, version 4) or Grace 5 (`@target G0.S0`) project,
gzip-compressed or not. Returns a vector of [`XmgrGraph`](@ref).
"""
function read_xmgr(path::AbstractString)
    graphs = XmgrGraph[]
    index = Dict{String,Int}()
    function graph!(name::AbstractString)
        name = String(name)
        haskey(index, name) && return graphs[index[name]]
        push!(graphs, XmgrGraph(name))
        index[name] = length(graphs)
        return graphs[end]
    end
    gname = "g0"
    graph!(gname)
    cur = nothing
    currows = Vector{Vector{Float64}}()
    function flush_set!()
        cur === nothing && return
        data = isempty(currows) ? Matrix{Float64}(undef, 0, 2) :
               permutedims(reduce(hcat, currows))
        push!(graph!(gname).sets, XmgrSet(gname, cur.name, cur.type, data, cur.legend,
                                          cur.props))
        cur = nothing
        currows = Vector{Vector{Float64}}()
    end
    for raw in read_lines(path)
        s = strip(raw)
        isempty(s) && continue
        if (m = match(r"(?i)^@with\s+(g\d+)", s)) !== nothing
            flush_set!()
            gname = lowercase(m.captures[1])
            graph!(gname)
        elseif (m = match(r"(?i)^@target\s+(?:G(\d+)\.)?S(\d+)", s)) !== nothing
            flush_set!()
            m.captures[1] !== nothing && (gname = "g" * m.captures[1])
            graph!(gname)
            cur = (name = "S" * m.captures[2], type = "xy", legend = "",
                   props = Dict{String,String}())
        elseif (m = match(r"(?i)^@type\s+(\S+)", s)) !== nothing && cur !== nothing
            cur = merge(cur, (type = m.captures[1],))
        elseif (m = match(r"^@\s+s(\d+)\s+legend\s+\"?(.*?)\"?$", s)) !== nothing &&
               cur !== nothing
            cur = merge(cur, (legend = m.captures[2],))
        elseif (m = match(r"^@\s+s\d+\s+(symbol|line type|color|symbol size|symbol fill)\s+(\S+)",
                          s)) !== nothing && cur !== nothing
            cur.props[m.captures[1]] = m.captures[2]
        elseif (m = match(r"^@\s+title\s+\"(.*)\"", s)) !== nothing
            g = graph!(gname)
            isempty(g.title) && (g.title = m.captures[1])
        elseif (m = match(r"^@\s+xaxis\s+label\s+\"(.*)\"", s)) !== nothing
            graph!(gname).xlabel = m.captures[1]
        elseif (m = match(r"^@\s+yaxis\s+label\s+\"(.*)\"", s)) !== nothing
            graph!(gname).ylabel = m.captures[1]
        elseif (m = match(r"^@g\d+\s+type\s+(\S+)", s)) !== nothing
            graph!(gname).logx = occursin("x", m.captures[1])
            graph!(gname).logy = occursin("y", m.captures[1])
        elseif (m = match(r"world\s+(xmin|xmax|ymin|ymax)\s+(\S+)", s)) !== nothing
            graph!(gname).world[m.captures[1]] = parse(Float64, m.captures[2])
        elseif s == "&"
            flush_set!()
        elseif !startswith(s, "@") && cur !== nothing
            vals = tryparse.(Float64, split(s))
            all(v -> v !== nothing, vals) && push!(currows, Float64[v for v in vals])
        end
    end
    flush_set!()
    return graphs
end

"""
    save_fig(fig, name; dir = ...)

Save `fig` as both PNG and PDF under `figures/out/`.
"""
function save_fig(fig, name::AbstractString; dir::AbstractString = joinpath(@__DIR__, "out"))
    mkpath(dir)
    save(joinpath(dir, name * ".png"), fig)
    save(joinpath(dir, name * ".pdf"), fig)
    return joinpath(dir, name)
end

"""
    figure_data(name)

Return the path of a small historical data file distributed with the figure
scripts. See `figures/data/README.md` for its provenance.
"""
figure_data(name::AbstractString) = joinpath(@__DIR__, "data", name)
