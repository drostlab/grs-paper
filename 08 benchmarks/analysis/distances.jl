using Arrow
using CSV
using DataFrames
using Random
using Statistics

trajectory(path, nested) = nested ? occursin('/', path) : true

function tagged(files::AbstractVector)
    out = Dict{Tuple{String, Int32}, Tuple{Vector{Int}, Vector{Int64}}}()
    ids = Dict{String, Int}()
    for (chunk, file) in enumerate(files)
        t = Arrow.Table(file)
        nested = any(p -> occursin('/', String(p)), t.path)
        for i in 1:length(t.name)
            path = String(t.path[i])
            trajectory(path, nested) || continue
            id = get!(ids, string(chunk, ':', path), length(ids) + 1)
            tags, values = get!(() -> (Int[], Int64[]), out, (String(t.name[i]), t.sample[i]))
            push!(tags, id)
            push!(values, Int64(t.value[i]))
        end
    end
    out, length(ids)
end

subset(data, keep) = Dict(k => values[[t in keep for t in tags]] for (k, (tags, values)) in data)

function cells(files::AbstractVector)
    data, n = tagged(files)
    Dict(k => values for (k, (_, values)) in data), n
end

cells(file::AbstractString) = cells([file])

smd(a, b) = (mean(a) - mean(b)) / std(b)

function histogram_distance(a, b)
    lo, hi = min(minimum(a), minimum(b)), max(maximum(a), maximum(b))
    hi > lo || return 0.0
    n = min(hi - lo + 1, 200)
    width = (hi - lo + 1) / n
    pa, pb = zeros(n), zeros(n)
    for v in a
        pa[min(n, Int(fld(v - lo, width)) + 1)] += 1 / length(a)
    end
    for v in b
        pb[min(n, Int(fld(v - lo, width)) + 1)] += 1 / length(b)
    end
    sum(abs, pa .- pb) / 2
end

function wasserstein(x, y)
    x, y = sort(x), sort(y)
    n, m = length(x), length(y)
    i = j = 1
    total = 0.0
    prev = float(min(x[1], y[1]))
    while i <= n || j <= m
        curr = float(i <= n ? (j <= m ? min(x[i], y[j]) : x[i]) : y[j])
        total += abs((i - 1) / n - (j - 1) / m) * (curr - prev)
        while i <= n && x[i] == curr; i += 1 end
        while j <= m && y[j] == curr; j += 1 end
        prev = curr
    end
    total
end

function kolmogorov(x, y)
    x, y = sort(x), sort(y)
    n, m = length(x), length(y)
    i = j = 1
    d = 0.0
    while i <= n && j <= m
        curr = min(x[i], y[j])
        while i <= n && x[i] <= curr; i += 1 end
        while j <= m && y[j] <= curr; j += 1 end
        d = max(d, abs(i / n - j / m))
    end
    d
end

function compare(reference, arm)
    deltas, hist, wass, ks = Float64[], Float64[], Float64[], Float64[]
    for (key, ref) in reference
        haskey(arm, key) || continue
        a = arm[key]
        s = std(ref)
        s > 0 || continue
        push!(deltas, smd(a, ref))
        push!(hist, histogram_distance(a, ref))
        push!(wass, wasserstein(a, ref) / s)
        push!(ks, kolmogorov(a, ref))
    end
    (signed_smd = mean(deltas), abs_smd = mean(abs, deltas), hist = mean(hist),
     wasserstein = mean(wass), ks = mean(ks), cells = length(deltas))
end

function null(data, total, n, keys; seed = 1)
    2n <= total || error(
        "reference has $total trajectories, need $(2n) for a disjoint n=$n null")
    shuffled = shuffle(MersenneTwister(seed), 1:total)
    a = subset(data, BitSet(shuffled[1:n]))
    b = subset(filter(p -> first(p) in keys, data), BitSet(shuffled[n+1:2n]))
    compare(b, a)
end

function sample(data, total, n; seed = 0)
    n <= total || error("reference has $total trajectories, need $n to match the arm")
    subset(data, BitSet(shuffle(MersenneTwister(seed), 1:total)[1:n]))
end

const RUN = r"^(.+)-n(\d+)-s(\d+)$"

function chunks(dirs)
    out = Dict{String, Vector{String}}()
    for dir in sort(collect(dirs))
        m = match(RUN, basename(dir))
        m === nothing && error("unparseable run path: $dir")
        file = joinpath(dir, "counts.arrow")
        push!(get!(Vector{String}, out, m[1]), file)
    end
    out
end

function score(dirs; genes, reference = "grs-ssa", replicates = 0)
    found = chunks(dirs)
    haskey(found, reference) || error("no $reference among the planned runs")
    ref_files = found[reference]
    ref_data, ref_total = tagged(ref_files)
    matched = Dict{Int, Any}()
    rows = DataFrame()
    armkeys = Dict{Int, Set{Tuple{String, Int32}}}()
    for arm in sort!(collect(keys(found)))
        arm == reference && continue
        arm_cells, ntraj = cells(found[arm])
        union!(get!(Set{Tuple{String, Int32}}, armkeys, ntraj), keys(arm_cells))
        push!(rows, merge((genes = genes, arm = arm, trajectories = ntraj,
                           replicate = missing), compare(get!(() -> sample(ref_data, ref_total, ntraj), matched, ntraj), arm_cells));
              cols = :union)
        println(arm, "  n=", ntraj, "  signed ", rows.signed_smd[end])
        flush(stdout)
    end
    for (n, keep) in sort!(collect(armkeys); by = first), replicate in 1:replicates
        2n <= ref_total || continue
        push!(rows, merge((genes = genes, arm = "null", trajectories = n,
                           replicate = replicate), null(ref_data, ref_total, n, keep; seed = replicate));
              cols = :union)
    end
    rows
end

function main(args)
    out, genes = args[1], parse(Int, args[2])
    replicates, reference, dirs = 0, "grs-ssa", String[]
    i = 3
    while i <= length(args)
        if args[i] == "--nulls"
            replicates = parse(Int, args[i+1]); i += 2
        elseif args[i] == "--reference"
            reference = args[i+1]; i += 2
        else
            push!(dirs, args[i]); i += 1
        end
    end
    df = score(dirs; genes, reference, replicates)
    mkpath(dirname(out))
    CSV.write(out, df)
    println("wrote ", out, "  ", nrow(df), " rows")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 3 || error(
        "usage: julia --project=. analysis/distances.jl <out.csv> <genes> " *
        "[--nulls N] [--reference ARM] <run-dir>...")
    main(ARGS)
end
