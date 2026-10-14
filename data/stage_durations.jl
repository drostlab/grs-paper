using XLSX, DataFrames, Statistics, CSV, CairoMakie


# Guan G, Li Z, Ma Y, Ye P, Cao J, Wong MK, Ho VWS, Chan LY, Yan H, Tang C, Zhao Z (2025).
# Cell lineage-resolved embryonic morphological map reveals signaling associated with cell fate and size asymmetry.
# Nat Commun 16, 3700. doi:10.1038/s41467-025-58878-0
# 8 wildtype embryos imaged at 20 °C
# because the Burght2024 omics data annotates cell cycle stage (early/late)
# we split the stages into early and late, assuming early is the first third of each cycle
# the omics data also used embryos grown at 18 °C
# to correct for this, we scale the durations with the Arrhenius law, resulting in x1.24 higher durations than the imaged ones
# Begasse ML, Leaver M, Vazquez F, Grill SW, Hyman AA (2015). Temperature dependence of cell division timing accounts for a shift in the thermal limits of C. elegans and C. briggsae. Cell Reports 10(5):647–653. doi:10.1016/j.celrep.2015.01.006

url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-58878-0/MediaObjects/41467_2025_58878_MOESM3_ESM.zip"
raw = "data/raw/guan2025"
files = ["Supplementary Data 3.xlsx", "Supplementary Data 4.xlsx"]
all(isfile, joinpath.(raw, files)) || (mkpath(raw); run(`unzip -oq $(download(url)) $files -d $raw`))

clean(x) = rstrip(string(x), '\'')

resolution = let
    t = XLSX.readtable(joinpath(raw, "Supplementary Data 3.xlsx"), 1; first_row = 2) |>
        DataFrame
    Dict(clean(r[1]) => Float64(r[8]) for r in eachrow(t) if startswith(clean(r[1]), "WT_"))
end

division = let
    t = XLSX.readtable(joinpath(raw, "Supplementary Data 4.xlsx"), "Cell proceeding cytokinesis";
                       first_row = 2, stop_in_empty_row = false) |> DataFrame |> dropmissing
    d = DataFrame(cell = clean.(t[!, 1]), embryo = clean.(t[!, 2]), t = Int.(t[!, 3]))
    combine(groupby(d, [:cell, :embryo]), :t => maximum => :t)
end

divided(cells, embryo) = let ts = [r.t for r in eachrow(division) if r.cell in cells && r.embryo == embryo]
    isempty(ts) ? NaN : mean(ts) * resolution[embryo]
end

cycles = ["EMS"  => (["P1"], ["EMS"]),
          "E"    => (["EMS"], ["E"]),
          "Ex"   => (["E"], ["Ea", "Ep"]),
          "Exx"  => (["Ea", "Ep"], ["Eal", "Ear", "Epl", "Epr"]),
          "Exxx" => (["Eal", "Ear", "Epl", "Epr"], ["Eala", "Ealp", "Eara", "Earp", "Epla", "Eplp", "Epra", "Eprp"]),
          "MS"   => (["EMS"], ["MS"]),
          "MSx"  => (["MS"], ["MSa", "MSp"]),
          "MSxx" => (["MSa", "MSp"], ["MSaa", "MSap", "MSpa", "MSpp"])]

lengths = DataFrame([(; cycle, embryo, minutes = divided(cells, embryo) - divided(parents, embryo))
                     for (cycle, (parents, cells)) in cycles for embryo in sort(collect(keys(resolution)))])
filter!(r -> isfinite(r.minutes), lengths)

cycle_length = combine(groupby(lengths, :cycle),
    :minutes => mean => :mean,
    :minutes => std => :sd,
    nrow => :n)
L = Dict(r.cycle => r.mean for r in eachrow(cycle_length))

stages = DataFrame(stage = ["earlyEMS", "EMS", "earlyEandMS",
                            "E", "earlyEx", "Ex", "earlyExx", "Exx", "Exxx",
                            "MS", "earlyMSx", "MSx", "earlyMSxx", "MSxx"],
                   minutes = [L["EMS"] / 3, 2L["EMS"] / 3, L["E"] / 3,
                              2L["E"] / 3, L["Ex"] / 3, 2L["Ex"] / 3, L["Exx"] / 3, 2L["Exx"] / 3, L["Exxx"],
                              L["MS"] - L["E"] / 3, L["MSx"] / 3, 2L["MSx"] / 3, L["MSxx"] / 3, 2L["MSxx"] / 3])

# Arrhenius rescaling from 20 °C to 18 °C, Eₐ = 75 kJ/mol for early embryonic cells:
# Rombouts J, Tavella F, Vandervelde A, Phong C, Ferrell JE, Yang Q, Gelens L (2025). Mechanistic origins of temperature scaling in the early embryonic cell cycle. Nature Communications 16:8045. doi:10.1038/s41467-025-62918-0
arrhenius(Ea, T1, T2) = exp(Ea / 8.314 * (1 / (T1 + 273.15) - 1 / (T2 + 273.15)))
stages.minutes .*= arrhenius(75e3, 18, 20)

stages.seconds = round.(Int, 60 .* stages.minutes)
push!(stages, (stage = "Exxxx", minutes = 90.0, seconds = 5400))
CSV.write("data/stage_durations.tsv", stages; delim = '\t')

# plot the stage duration distribution across embryos

stage_colour = Dict(String(r.stage) => String(r.colour) for r in CSV.File("data/stage_colours.tsv"))
stage_lineage = Dict(String(r.stage) => String(r.lineage) for r in CSV.File("data/stage_colours.tsv"))


children(c) = c == "EMS" ? ["E", "MS"] : c in ("Ea", "Ep") ? [c * "l", c * "r"] : [c * "a", c * "p"]
parent_of = Dict(c => (length(ps) == 1 ? only(ps) : only(filter(p -> startswith(c, p), ps)))
                 for (_, (ps, cs)) in cycles for c in cs)
cycle_of = Dict(c => k for (k, (_, cs)) in cycles for c in cs)
origin(e) = divided(["EMS"], e)

function draw!(ax, c::String, e, y, h; depth = 0)
    s, t = divided([parent_of[c]], e) - origin(e), divided([c], e) - origin(e)
    yc = y + h / 2
    isfinite(s) && isfinite(t) && lines!(ax, [s, t], [yc, yc]; color = stage_colour[cycle_of[c]], linewidth = max(5 - depth, 1.5))
    haskey(parent_of, first(children(c))) || return
    isfinite(t) && lines!(ax, [t, t], [y + 3h / 4, y + h / 4]; color = :gray40, linewidth = 0.6)
    for (i, d) in enumerate(children(c))
        draw!(ax, d, e, y + (i - 1) * h / 2, h / 2; depth = depth + 1)
    end
end

function draw!(ax, stages::DataFrame, y, h; shift = 0.0)
    m = Dict(r.stage => r.minutes for r in eachrow(stages))
    yE, yMS = y + 0.75h, y + 0.25h
    ems = m["earlyEMS"] + m["EMS"]
    segments = [("earlyEMS", -ems, y + h / 2), ("EMS", -ems + m["earlyEMS"], y + h / 2)]
    for (lineage, yl) in (("E", yE), ("MS", yMS))
        push!(segments, ("earlyEandMS", 0.0, yl))
        t = m["earlyEandMS"]
        for s in filter(s -> stage_lineage[s] == lineage, stages.stage)
            push!(segments, (s, t, yl))
            t += m[s]
        end
    end
    lines!(ax, [shift, shift], [yMS, yE]; color = :gray40, linewidth = 0.6)
    for (s, t, yl) in segments
        t += shift
        lines!(ax, [t, t + m[s]], [yl, yl]; color = stage_colour[s], linewidth = 6)
        text!(ax, t + m[s] / 2, yl; text = s, align = (:center, :bottom), offset = (0, 4), fontsize = 7, font = :bold)
    end
end

embryos = sort(collect(keys(resolution)))
f = Figure(size = (1400, 900))
ax = Axis(f[1, 1]; xlabel = "minutes after EMS division", yticks = (1:length(embryos), embryos), yreversed = true)
foreach(((i, e),) -> draw!(ax, "EMS", e, i - 0.4, 0.8), enumerate(embryos))
mkpath("figures")
save("figures/stage_durations.svg", f)

f = Figure(size = (900, 220))
ax = Axis(f[1, 1]; xlabel = "minutes after EMS birth", xticks = 0:10:230, yticksvisible = false,
          yticklabelsvisible = false, ygridvisible = false)
draw!(ax, stages, 0, 1; shift = stages.minutes[1] + stages.minutes[2])
ylims!(ax, 0, 1.15)
save("figures/stage_model.svg", f)
