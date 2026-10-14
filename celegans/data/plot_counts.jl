using CSV, DataFrames, CairoMakie, Statistics

counts = CSV.read("data/ems_counts.tsv", DataFrame; stringtype = String)
stages = CSV.read("data/stage_durations.tsv", DataFrame; stringtype = String).stage
stage_table = CSV.read("data/stage_colours.tsv", DataFrame; stringtype = String)
lineage = Dict(zip(stage_table.stage, stage_table.lineage))
gene_table = CSV.read("data/gene_colours.tsv", DataFrame; stringtype = String)
gene_colour = Dict(zip(gene_table.gene, gene_table.colour))
genes = names(counts)[3:end]

gap(a, b) = lineage[a] == lineage[b] ? 1.0 : 1.6
position = Dict(zip(stages, cumsum([1.0; gap.(stages[1:end-1], stages[2:end])])))
cells(s, g) = counts[counts.cell_type .== s, g]
labels = ["$s\n($(count(==(s), counts.cell_type)))" for s in stages]

bins = [0, 1, 2, 3, 5, 9, 17, 33]
binlabels = ["0", "1", "2", "3–4", "5–8", "9–16", "17–32", "33+"]
bin(c) = searchsortedlast(bins, c)
fractions(s, g) = [mean(==(b), bin.(cells(s, g))) for b in eachindex(bins)]

f = Figure(size = (1100, 1500), fontsize = 14)
axes = map(enumerate(genes)) do (i, g)
    ax = Axis(f[i, 1]; ylabel = g, ylabelfont = :bold, yticks = (eachindex(bins), binlabels), yticklabelsize = 9,
              xticks = (getindex.(Ref(position), stages), labels), xticklabelrotation = π / 4,
              xticklabelsvisible = i == length(genes), xgridvisible = false, ygridvisible = false)
    for s in stages
        fr = fractions(s, g)
        b = findall(>(0), fr)
        scatter!(ax, fill(position[s], length(b)), b; markersize = 26 .* sqrt.(fr[b]),
                 color = gene_colour[g], strokecolor = :gray30, strokewidth = 0.5)
    end
    ax
end
linkxaxes!(axes...)
rowgap!(f.layout, 4)
mkpath("figures")
save("figures/ems_counts.svg", f)
