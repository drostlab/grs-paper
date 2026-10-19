### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 8d6d12d2-97d0-42b2-8340-a8acf112e60f
begin
	import Pkg; Pkg.activate(".")

	import GeneRegulatorySystems as GRS

	using CSV, JSON
	
	using Statistics
	import Random

	using CairoMakie
	CairoMakie.activate!()
	
	using PlutoUI
	using Colors
	import CytoscapeJS
	CytoscapeJS.Bonito.Page()

	TableOfContents()
end

# ╔═╡ 09291836-bddd-11f1-aaf5-cd7c7aaa7bc3
md"""
# Modelling the specification of the C. elegans endomesoderm via GeneRegulatorySystems.jl
"""

# ╔═╡ afd5f0c5-eeb2-429b-848c-b055947be9cd
begin
	show_schedule(x::AbstractString; kwargs...) = show_schedule(JSON.parse(x); kwargs...)
	show_schedule(x::AbstractDict; title = "schedule.json") = let
		s = sprint(io -> JSON.print(io, x, 4))
		HTML("""
		<details open>
		  <summary style="cursor:pointer;font-family:monospace;font-weight:300;font-size:0.85em;padding:6px 0">$(title)</summary>
		  <div style="position:relative">
		    <button onclick="navigator.clipboard.writeText(this.nextElementSibling.innerText);this.innerText='✓'"
		            style="position:absolute;top:6px;right:6px;font-size:0.75em;cursor:pointer">copy</button>
		    <pre style="max-height:400px;overflow:auto;font-size:0.8em;margin:0">$(replace(s, "&" => "&amp;", "<" => "&lt;"))</pre>
		  </div>
		</details>
		""")
	end

end

# ╔═╡ 7a93e7dc-4b76-44c6-886a-24d5e8103f15
md"""
## 1. Setup schedule
"""

# ╔═╡ 29d8dde1-8e4b-4628-bdf6-87cfef5e0a4f
md"""
### 1.1. Define the model
"""

# ╔═╡ ed2e5dc7-0489-47d1-adcd-450c239ea374
md"""
#### Base rates

Most of the reaction rates are fixed using information from the literature, and shared across genes.

"""

# ╔═╡ 11fe5a08-04f3-4073-b553-b98a9c980cbd
shared_machinery = Dict(
	:polymerases => 500000,
	:ribosomes => 20000000,
	:proteasomes => 1000000,
)

# ╔═╡ 28e41292-2cbc-4234-a4b3-5d4a2cacb1f4
decay(half_life_minutes) = log(2) / (60half_life_minutes)

# ╔═╡ 187cd4b3-2a64-4ecf-bcc6-24450b5107bb
base_rates = Dict(
	:activation => 0.015,
	:deactivation => 10.0,
	:initiation => 0.2 / shared_machinery[:polymerases],
	:abortion => 4e-4,
	:processing => 1 / 120,
	:premrna_decay => 1e-3,
	:translation => 0.06 / shared_machinery[:ribosomes],
	:mrna_decay => decay(20),
	:protein_decay => decay(20) / shared_machinery[:proteasomes],
)

# ╔═╡ ebfbb616-18d1-4039-a193-b8c0c6383f5d
gene_specific_rates = let
	transient = decay(5)
	degraded, stable = (decay(4.5), decay(120)) ./ shared_machinery[:proteasomes]
	Dict(
		"skn-1" => (activation = 0, mrna_decay = transient, protein_decay = degraded),
		"pop-1" => (activation = 0, protein_decay = 0.0),
		"ref-1" => (activation = 0, protein_decay = stable),
		"med-2" => (mrna_decay = transient, protein_decay = degraded),
		"med-1" => (mrna_decay = transient, protein_decay = degraded),
		"end-3" => (mrna_decay = transient,),
		"end-1" => (protein_decay = stable,),
		"elt-7" => (protein_decay = stable,),
		"elt-2" => (protein_decay = stable,)
	)
end

# ╔═╡ 4a0c8c51-0ac6-4bcb-8681-8b37aec85978
md"""
#### Gene Regulatory Network
"""

# ╔═╡ 693d942c-6270-4ec3-a519-b832600a8f86
gene_colours = Dict(r.gene => parse(Colorant, r.colour) for r in CSV.File("data/gene_colours.tsv"; stringtype = String))

# ╔═╡ 16292fb4-24f9-4041-949b-b4c58c70065b
signalling = [
	Dict(:name => "pop-1.sys-1_binding", :from => ["pop-1", "P2-signal"], :to => ["pop-1.proteins*", "P2-signal"], :rate => 4.9e-6),
	Dict(:name => "pop-1.sys-1_unbinding", :from => ["pop-1.proteins*"], :to => ["pop-1"], :rate => 0.0019),
	Dict(:name => "ref-1.activation", :from => ["ref-1.proteins", "P2-signal"], :to => ["ref-1.proteins*", "P2-signal"], :rate => 4.9e-6),
	Dict(:name => "ref-1.deactivation", :from => ["ref-1.proteins*"], :to => ["ref-1.proteins"], :rate => 0.017),
	Dict(:name => "pop-1.proteins*_decay", :from => ["pop-1.proteins*", "proteasomes"], :to => ["proteasomes"], :rate => gene_specific_rates["pop-1"].protein_decay),
	Dict(:name => "ref-1.proteins*_decay", :from => ["ref-1.proteins*", "proteasomes"], :to => ["proteasomes"], :rate => gene_specific_rates["ref-1"].protein_decay),
]

# ╔═╡ cab2d623-8e0d-415f-8192-7e14b22249c7
Base.show(io::IO, m::MIME"text/html", g::CytoscapeJS.Cytoscape) = show(io, m, CytoscapeJS.Bonito.App(g))

# ╔═╡ e1622993-d86a-4b44-9c25-1dc9a3cf7cc7
md"""
### 1.2. Define the wildtype schedule for the endomesoderm specification
"""

# ╔═╡ 361a46e8-cadc-456c-8aec-67ebd07205f9
stage_colours = Dict(String(r.stage) => parse(Colorant, r.colour) for r in CSV.File("data/stage_colours.tsv"))

# ╔═╡ 9a93134d-84cd-49f9-a337-cf688d491728
PlutoUI.LocalResource("figures/stage_durations.svg")

# ╔═╡ f89e6f07-a453-4d56-9719-8ab4cb89c5f7
stage_duration = Dict(r.stage => r.seconds for r in CSV.File("data/stage_durations.tsv"; stringtype = String))

# ╔═╡ 18571daf-95b1-4293-a9c7-6e77f91b87a0
copies(genes) = Dict(Symbol(g, ".inactive") => 2 for g in genes)

# ╔═╡ 8015c70a-c431-4ae4-afb8-c49f04ae86a4
atlas_spec = let
	num_reps = 20
	sample_rate = 50
	use(name) = Dict(Symbol("\$") => name)
	hexcolour(stage) = "#" * hex(stage_colours[stage])

	capture = Dict(:branch => true, :step => [Dict(:step => [
		Dict(Symbol("{filter}") => "\\.premrnas\$"),
		Dict(Symbol("{resample-each-binomial}") => 0.2),
		Dict(:into => "", :step => Dict(Symbol("{pass}") => true))
	])])
	stage_spec(name) = Dict(:stage => name, :to => stage_duration[name], :step => [use("sample_rate"), capture])
	lineage(name, steps) = Dict(:label => "atlas rep:\${rep} $name lineage", :color => hexcolour(name), :step => steps)
	Dict(
		:seed => "\${rootseed}-sc",
		:into => nothing,
		:sample_rate => sample_rate,
		:do => use("EMS_circuit"),
		:each => Dict(:length => num_reps), :as => "rep", :branch => true,
		:step => [
			Dict(Symbol("{add}") => copies(keys(gene_colours))),
			lineage("EMS", stage_spec.(["earlyEMS", "EMS"])),
			Dict(:branch => true, :step => [
				lineage("MS", stage_spec.(["earlyEandMS", "MS", "earlyMSx", "MSx", "earlyMSxx", "MSxx"])),
				lineage("E", [use("P2_signalling"); stage_spec.(["earlyEandMS", "E", "earlyEx", "Ex", "earlyExx", "Exx", "Exxx"])]),
			]),
		],
	)
end

# ╔═╡ 2f30332d-debf-499a-979e-20a20cc47fae
md"""
### 1.3 Define knockout experiments
"""

# ╔═╡ 0652f10d-d855-4adb-958f-0d8716e2cea6
struct Readout
	gene::String
	stage::String
end

# ╔═╡ 6c2949ef-6309-45b1-a957-df3beae4a22e
Base.show(io::IO, r::Readout) = print(io, uppercase(r.gene), "@", r.stage)

# ╔═╡ dc4d3bba-e07d-48d8-a45f-26e5c4e2c8a2
knockout_data = [(; readout = Readout(r.readout, r.stage), r.mean, r.sd, knockout = sort(String.(split(coalesce(r.knockout, ""), ",", keepempty = false))))
                 for r in CSV.File("data/protein_abundance.tsv"; stringtype = String, missingstring = ["", "NA"])]

# ╔═╡ 87bf13f3-e71b-4a8d-826e-61bc1c635039
knockouts = unique(r.knockout for r in knockout_data)

# ╔═╡ 616b8a62-6edc-4d6a-83c7-b2bffc9a368c
readouts = unique(r.readout for r in knockout_data)

# ╔═╡ 7ce0250c-7a68-44dd-86b8-cfff20c53471
rnai_genes = ["elt-7", "elt-2"]

# ╔═╡ ddbe9327-7b98-404f-8b79-49e1bd074c8b
rnai = [Dict(:name => "$g.rnai", :from => ["$g.mrnas", "$g.sirna"], :to => ["$g.sirna"], :rate => decay(20)) for g in rnai_genes]

# ╔═╡ eabe7e00-bf3f-484f-b64d-980c2d05b4f1
genotype(ko) = Dict(Symbol("{add}") => merge(
	copies(setdiff(keys(gene_colours), setdiff(ko, rnai_genes))),
	Dict(Symbol(g, ".sirna") => 1 for g in ko ∩ rnai_genes)))

# ╔═╡ 91114341-c844-45b9-8679-24079cc65ce3
function plot_knockouts(data)
	f = Figure(size = (1000, 300))
	for (i, readout) in enumerate(unique(r.readout for r in data))
		rows = sort([r for r in data if r.readout == readout]; by = r -> length(r.knockout))
		x = eachindex(rows)
		ax = Axis(f[1, i]; title = string(readout), ylabel = i == 1 ? "protein (a.u.)" : "",
		          xticks = (x, [isempty(r.knockout) ? "WT" : join(r.knockout, "+") for r in rows]), xticklabelrotation = π / 4)
		barplot!(ax, x, [r.mean for r in rows]; color = (gene_colours[readout.gene], 0.7))
		with_sd = [i for (i, r) in enumerate(rows) if !ismissing(r.sd)]
		isempty(with_sd) || errorbars!(ax, with_sd, [rows[i].mean for i in with_sd], [rows[i].sd for i in with_sd]; whiskerwidth = 6, color="gray")
	end
	f
end

# ╔═╡ 80458248-05b1-4380-a1ca-89125da04da0
knockout_spec = let
	use(name) = Dict(Symbol("\$") => name)
	num_reps = 80
	traced = Dict(:branch => true, :step => [Dict(:step => [
		Dict(Symbol("{filter}") => "\\.proteins\$"),
		Dict(Symbol("{resample-each-binomial}") => 0.2),
		Dict(:into => "", :step => Dict(Symbol("{pass}") => true)),
	])])
	stage_spec(name) = Dict(:stage => name, :step => [
		Dict(:to => stage_duration[name], :into => nothing), traced])
	name(ko) = isempty(ko) ? "WT" : join(ko, "+")
	colour(ko) = "#" * hex(isempty(ko) ? stage_colours["E"] : mapreduce(g -> gene_colours[g], (a, b) -> weighted_color_mean(0.5, a, b), ko))
	experiment(ko) = Dict(
		:label => "knockouts ko:$(name(ko))", :color => colour(ko), :knockouts => ko,
		:step => Dict(
			:do => use("EMS_circuit"),
			:each => Dict(:length => num_reps), :as => "rep", :branch => true,
			:step => [genotype(ko); stage_spec.(["earlyEMS", "EMS"]); use("P2_signalling"); stage_spec.(["earlyEandMS", "E", "earlyEx", "Ex", "earlyExx", "Exx", "Exxx", "Exxxx"])],
		),
	)
	Dict(:seed => "\${rootseed}-ko", :branch => true, :into => nothing, :step => experiment.(knockouts))
end

# ╔═╡ a23dfd35-fcc4-4544-973e-8539c5725e33
begin
	function plot_stage_model!(pos; examples = [["elt-2", "elt-7"], ["elt-7", "end-1"]])
		table = CSV.File("data/stage_colours.tsv"; stringtype = String)
		lineage = Dict(r.stage => r.lineage for r in table)
		m = Dict(s => t / 60 for (s, t) in stage_duration)
		ems = m["earlyEMS"] + m["EMS"]
		yMS, yEMS, yE = 0.84, 0.68, 0.52
		yKO = [0.28, 0.07]
		start(y) = [("earlyEMS", 0.0, y), ("EMS", m["earlyEMS"], y)]
		function track(l, y, last)
			segments = [("earlyEandMS", ems, y)]
			t = ems + m["earlyEandMS"]
			for r in table
				lineage[r.stage] == l || continue
				push!(segments, (r.stage, t, y))
				t += m[r.stage]
				r.stage == last && break
			end
			segments
		end
		atlas = [(start(yEMS), true), (track("MS", yMS, "MSxx"), true), (track("E", yE, "Exxx"), false)]
		kos = [[start(y); track("E", y, "Exxxx")] for y in yKO]
		stop(segments, stage) = only(t + m[s] for (s, t, _) in segments if s == stage)
		finish = stop(kos[1], "Exxxx")
		ax = Axis(pos; xlabel = "minutes after EMS cell birth", xticks = 0:20:finish, xgridvisible = false, ygridvisible = false)
		hidespines!(ax)
		hideydecorations!(ax)
		lines!(ax, [0, 1.04finish], [0, 0]; color = :black, linewidth = 1)
		scatter!(ax, [1.04finish], [0]; marker = :rtriangle, color = :black, markersize = 12)
		measurement = (marker = :diamond, color = :white, strokecolor = :gray20)
		readout = (fontsize = 14, font = :italic, color = :gray15)
		input = (fontsize = 14, font = :italic, color = colorant"#3f7fb0")

		event(x, lo, hi) = lines!(ax, [x, x], [lo, hi]; color = (colorant"#6fa8dc", 0.6), linewidth = 1.0)
		top = 1.02
		event(0, yKO[end] - 0.05, top)
		text!(ax, 0, top; text = "+ initial SKN-1, POP-1, REF-1", align = (:left, :top), offset = (4, 0), input...)
		event(ems, yKO[end] - 0.05, yE + 0.3)
		text!(ax, ems, (yE + yKO[1]) / 2 - 0.03; text = "+ P2 signal", align = (:left, :center), offset = (5, 0), input...)

		lines!(ax, [ems, ems], [yE, yMS]; color = :gray25, linewidth = 1.4)
		for (s, t, y) in [reduce(vcat, first.(atlas)); reduce(vcat, kos)]
			lines!(ax, [t, t + m[s]], [y, y]; color = stage_colours[s], linewidth = 9)
		end
		width(s) = 2.2 * length(s)
		for (bar, above_first) in atlas
			a, b = bar[1][2], bar[end][2] + m[bar[end][1]]
			ts = collect(a + 3:6:b)
			scatter!(ax, ts, fill(bar[1][3], length(ts)); measurement..., strokewidth = 0.7, markersize = 8)
			previous = Dict(true => -Inf, false => -Inf)
			for (i, (s, t, y)) in enumerate(bar)
				above = isodd(i) == above_first
				x = t + m[s] / 2
				clash = x - width(s) / 2 < previous[above]
				clash || (previous[above] = x + width(s) / 2)
				shift = clash ? 23 : 7
				text!(ax, x, y; text = s, align = (:center, above ? :bottom : :top), offset = (0, above ? shift : -shift),
					fontsize = 11, font = :bold, color = weighted_color_mean(0.75, stage_colours[s], colorant"black"))
			end
		end
		text!(ax, stop(atlas[2][1], "MSxx"), yEMS; text = "◇ readout · pre-mRNA · 20% capture",
			align = (:left, :center), offset = (12, 0), readout...)

		genes(ko) = join([g in rnai_genes ? "$(g)ⁱ" : g for g in ko], "+")
		for (ko, y, segments) in zip(examples, yKO, kos)
			text!(ax, -6, y; text = rich(rich("KO"; font = :bold), "(", rich(genes(ko); font = :italic), ")"),
				align = (:right, :center), fontsize = 14)
			for r in knockout_data
				r.knockout == ko || continue
				x = stop(segments, r.readout.stage)
				scatter!(ax, [x], [y]; measurement..., strokewidth = 1.4, markersize = 14)
				text!(ax, x, y; text = uppercase(r.readout.gene) * " protein readout", align = (:center, :bottom),
					offset = (0, 11), readout...)
			end
		end
		text!(ax, -22, sum(yKO) / 2; text = "⋮", align = (:center, :center), fontsize = 18)

		for (text, reps, lo, hi) in (("snRNA-seq", atlas_spec[:each][:length], yE, yMS),
				("knockouts", knockout_spec[:step][1][:step][:each][:length], yKO[end], yKO[1]))
			lines!(ax, [-50, -50], [lo - 0.07, hi + 0.07]; color = :gray40, linewidth = 0.8)
			text!(ax, -52, (lo + hi) / 2; text, rotation = π / 2, align = (:center, :bottom), fontsize = 18, font = :bold)
			text!(ax, -48, (lo + hi) / 2; text = "×$reps reps", rotation = π / 2, align = (:center, :top), fontsize = 14)
		end
		xlims!(ax, -66, 1.06finish)
		ylims!(ax, -0.02, 1.12)
	end
	plot_stage_model() = (f = Figure(size = (1200, 420)); plot_stage_model!(f[1, 1]); f)
end

# ╔═╡ c5762b14-f2d4-4acd-aa10-0f8eda3a854b
plot_stage_model()

# ╔═╡ 6cd49d0d-353d-4b5c-a1a2-cb1376ad7d06
md"""
### 1.4 Assemble final schedule
"""

# ╔═╡ 79679362-4810-4b7b-ba82-aee0801a1943
maternal = Dict(Symbol("skn-1.proteins") => 30000, 
                Symbol("pop-1.proteins") => 30000, 
                Symbol("ref-1.proteins") => 30000
               )

# ╔═╡ 3ef37b74-b704-49bb-b3fe-55baf15a6ffb
build(spec) = GRS.Scheduling.prebuild(GRS.Models.parse(JSON.json(spec)), "/.EMS_circuit")

# ╔═╡ 6df7d3b4-adc1-4a38-996d-d460cd5881cc
md"""
### 1.5. Run simulation
"""

# ╔═╡ 89cdf95b-2c77-47c5-8fd9-c5ad2b27f20a
begin
	struct Sample
	    t::Float64
	    tags::Dict{Symbol, Any}
	    counts::Dict{Symbol, Float64}
	end
	is_tag(::Union{Number, AbstractString, Symbol}) = true
	is_tag(v::AbstractVector) = all(is_tag, v)
	is_tag(::Any) = false

	as_tags(bindings) = Dict{Symbol, Any}(k => bindings[k] for k in (:stage, :rep, :knockouts) if haskey(bindings, k))

	matches(s::Sample, selector::AbstractDict) = 
		all(get(s.tags, k, nothing) == v for (k, v) in selector)
	matches(s::Sample, selector) = selector(s.tags)::Bool

	samples_where(samples::Vector{Sample}, selector) = filter(s -> matches(s, selector), samples)
	species(samples::Vector{Sample}, gene, kind) = [get(s.counts, Symbol(gene, ".", kind), 0.0) for s in samples]
end;

# ╔═╡ ca1373c7-44b9-48f0-ac6d-a2db7529970c
begin
	function simulate(schedule!::GRS.Schedule; seed = "", kwargs...)
		schedule! = GRS.Scheduling.reseed(schedule!, seed)
		samples = Sample[]
		lock = ReentrantLock()
		function trace(state = nothing; primitive! = nothing, into = nothing, _...)
			(state === nothing || into === nothing) && return
			tags = primitive! === nothing ? Dict{Symbol, Any}() : as_tags(primitive!.bindings)
			@lock lock push!(samples, Sample(GRS.Models.t(state), tags, Dict{Symbol, Float64}(GRS.Models.counts(state))))
		end
		schedule!(GRS.Models.FlatState(); trace, kwargs...)
		sort!(samples; by = s -> s.t)
	end
	simulate(schedule!::GRS.Schedule, parameters::AbstractDict{Symbol}; kwargs...) =
		simulate(GRS.Models.remake(schedule!, "/.EMS_circuit", parameters); kwargs...)
end

# ╔═╡ f70bac03-3cb0-4cb2-af22-cbd52f9bd147
md"""
## 2. Observations
"""

# ╔═╡ 21334557-28ff-490e-8fa8-f5f4a2abc320
md"""
### 2.1. Transcriptomic data
"""

# ╔═╡ bbdfa49b-6bf5-42e3-87cd-20efe2e6b1a7
ems_counts = CSV.File("data/ems_counts.tsv"; stringtype = String)

# ╔═╡ d78cbb77-8286-47e1-9eeb-d0050ac8d61d
atlas_genes = string.(propertynames(ems_counts)[3:end])

# ╔═╡ 70975d2d-ebd2-470b-885d-6c64c1483f36
atlas_stages = unique(r.cell_type for r in ems_counts)

# ╔═╡ 95509860-902a-404c-ba76-c7c81da94895
PlutoUI.LocalResource("figures/ems_counts.svg")

# ╔═╡ ca8a141a-d0a3-4de0-8465-7a166823c82b
begin
	wildtype(tags) = haskey(tags, :stage) && !haskey(tags, :knockouts)

	atlas_counts(samples::Vector{Sample}, gene, stage) =
		species(samples_where(samples, t -> wildtype(t) && t[:stage] == stage), gene, :premrnas)
	atlas_counts(data::CSV.File, gene, stage) = [r[Symbol(gene)] for r in data if r.cell_type == stage]

	function transcriptome(source, genes, stages)
		f(stat) = [stat(atlas_counts(source, gene, stage)) for gene in genes, stage in stages]
		(; dropout = f(c -> mean(iszero, c)), meanlog = f(c -> mean(log1p, c)), stdlog = f(c -> std(log1p.(c))))
	end
end

# ╔═╡ ef3dc29c-a9d5-4147-9e9d-e6fd09146fc1
observed_transcriptome = transcriptome(ems_counts, atlas_genes, atlas_stages)

# ╔═╡ 0c4ee0a8-276f-4677-bbad-ebbdbd5cd50d
begin
	function plot_transcriptome!(pos, observed, simulated, genes)
		for (i, stat) in enumerate(keys(observed))
			ax = Axis(pos[1, i]; title = string(stat), xlabel = "observed", ylabel = i == 1 ? "simulated" : "", aspect = 1)
			for (g, gene) in enumerate(genes)
				scatter!(ax, observed[stat][g, :], simulated[stat][g, :]; color = gene_colours[gene])
			end
			ablines!(ax, 0, 1; color = :gray, linestyle = :dash)
		end
	end
	plot_transcriptome(args...) = (f = Figure(size = (900, 300)); plot_transcriptome!(f[1, 1], args...); f)
end

# ╔═╡ 88610f5a-b196-459a-b5da-502d0bb37f7a
md"""
### 2.2. Proteomic data under knockouts
"""

# ╔═╡ 8fb2685e-4db1-4d2e-8ff3-c92bdb173354
begin
	reference(readout, readouts) = readout

	function knockout_proteins(samples::Vector{Sample},
			readouts::Vector{Readout}, knockouts::Vector{Vector{String}}; β = 100)
		selects(ko, readout) = tags -> get(tags, :stage, "") == readout.stage && get(tags, :knockouts, nothing) == ko
		level(ko, readout) = mean(species(samples_where(samples, selects(ko, readout)), readout.gene, :proteins))
		wt = String[]
		(; logfoldchange = [log((level(ko, readout) + β) / (level(wt, reference(readout, readouts)) + β)) for readout in readouts, ko in knockouts])
	end

	function knockout_proteins(data::Vector{<:NamedTuple},
			readouts::Vector{Readout}, knockouts::Vector{Vector{String}})
		rows = Dict((d.readout, d.knockout) => d.mean for d in data)
		level(ko, readout) = get(rows, (readout, ko), NaN)
		wt = String[]
		(; logfoldchange = [log(level(ko, readout) / level(wt, reference(readout, readouts))) for readout in readouts, ko in knockouts])
	end
end

# ╔═╡ 5fdc5d3f-1830-4f9a-8073-5abb6fbd8ba8
observed_knockouts = knockout_proteins(knockout_data, readouts, knockouts)

# ╔═╡ 2f952564-8523-44d6-b224-48aac07561f9
function plot_knockouts(observed, simulated, readouts)
	f = Figure(size = (350, 350))
	ax = Axis(f[1, 1]; title = "log fold change", xlabel = "observed", ylabel = "simulated", aspect = 1)
	for (r, readout) in enumerate(readouts)
		scatter!(ax, observed.logfoldchange[r, :], simulated.logfoldchange[r, :]; color = gene_colours[readout.gene])
	end
	ablines!(ax, 0, 1; color = :gray, linestyle = :dash)
	f
end

# ╔═╡ 97d575cb-9efb-4cf6-8459-652ff9f348ee
begin
	function plot_heatmap!(pos, observed, simulated;
		genes = ["med-2", "med-1", "end-3", "end-1", "elt-7", "elt-2", "tbx-35", "ceh-51"],
			stages = ["earlyEMS", "EMS", "earlyEandMS", "E", "earlyEx", "Ex", "earlyExx", "Exx", "Exxx", "MS", "earlyMSx", "MSx", "earlyMSxx", "MSxx"])
		gi, si = indexin(genes, atlas_genes), indexin(stages, atlas_stages)
		for (i, stat) in enumerate((:dropout,))
			value(A) = stat == :dropout ? 1 .- A[gi, si] : A[gi, si]
			D, S = value(observed[stat]), value(simulated[stat])
			lo, hi = extrema(filter(isfinite, [vec(D); vec(S)]))
			shade(gene, x) = isfinite(x) ? weighted_color_mean(clamp((x - lo) / (hi - lo), 0.05, 1), gene_colours[gene], colorant"white") : colorant"grey90"
			ax = Axis(pos[1, i]; title = stat == :dropout ? "fraction expressing (top: data, bottom: simulation)" : "mean log expression",
				xticks = (eachindex(stages), stages), yticks = (eachindex(genes), genes), xticklabelrotation = π / 3,
				yreversed = true, xgridvisible = false, ygridvisible = false, xautolimitmargin = (0, 0), yautolimitmargin = (0, 0))
			cells = [(g, s) for g in eachindex(genes) for s in eachindex(stages)]
			poly!(ax, [Rect2f(s - 0.5, g - 0.42, 1, 0.42) for (g, s) in cells]; color = [shade(genes[g], D[g, s]) for (g, s) in cells])
			poly!(ax, [Rect2f(s - 0.5, g, 1, 0.42) for (g, s) in cells]; color = [shade(genes[g], S[g, s]) for (g, s) in cells])
			vlines!(ax, [findfirst(==("E"), stages), findfirst(==("MS"), stages)] .- 0.5; color = (:grey40, 0.5), linestyle = :dash, linewidth = 1)
		end
	end
	plot_heatmap(args...; kwargs...) = (f = Figure(size = (700, 400)); plot_heatmap!(f[1, 1], args...; kwargs...); f)
end

# ╔═╡ 61267d7f-a38d-4d31-a0af-0c5936c77021
begin
	function plot_knockouts!(pos, data, samples::Vector{Sample})
		level(readout, ko) = species(samples_where(samples, t -> get(t, :stage, "") == readout.stage && get(t, :knockouts, nothing) == ko), readout.gene, :proteins)
		for (i, readout) in enumerate(unique(r.readout for r in data))
			rows = sort([r for r in data if r.readout == readout && !isempty(r.knockout)]; by = r -> length(r.knockout))
			x = eachindex(rows)
			ref = reference(readout, readouts)
			wt_data = first(r.mean for r in data if r.readout == ref && isempty(r.knockout))
			wt_sim = mean(level(ref, String[]))
			sim = [level(readout, r.knockout) for r in rows]
			label(ko) = join([g in rnai_genes ? "$(g)ⁱ" : g for g in ko], "+")
			ax = Axis(pos[1, i]; title = string(readout), ylabel = i == 1 ? "protein relative to WT" : "",
			          xticks = (x, label.(r.knockout for r in rows)), xticklabelrotation = π / 4)

			hlines!(ax, [1]; color = (:gray, 0.6), linestyle = :dash, linewidth = 1)
			barplot!(ax, x .- 0.2, [r.mean / wt_data for r in rows]; width = 0.4, color = (gene_colours[readout.gene], 0.7))
			with_sd = [j for (j, r) in enumerate(rows) if !ismissing(r.sd)]
			isempty(with_sd) || errorbars!(ax, with_sd .- 0.2, [rows[j].mean / wt_data for j in with_sd], [rows[j].sd / wt_data for j in with_sd]; whiskerwidth = 6, color = :gray)
			barplot!(ax, x .+ 0.2, mean.(sim) ./ wt_sim; width = 0.4, color = (:grey50, 0.7))
			isempty(with_sd) || errorbars!(ax, with_sd .+ 0.2, mean.(sim[with_sd]) ./ wt_sim, std.(sim[with_sd]) ./ wt_sim; whiskerwidth = 6, color = :gray)
		end
	end
	plot_knockouts(data, samples::Vector{Sample}) = (f = Figure(size = (1000, 300)); plot_knockouts!(f[1, 1], data, samples); f)
end

# ╔═╡ 9770c432-429e-4321-baa2-06dbd7c8305b
plot_knockouts(knockout_data)

# ╔═╡ 115fad24-4ab7-4fc4-89ab-f6a32d99296a
edges = [
	(from = "skn-1", to = "med-2", kind = :activation, at = 950.0, k = -9.3),
	(from = "med-2", to = "med-1", kind = :activation, at = 4900.0, k = -5.9),
	(from = "skn-1", to = "med-1", kind = :activation, at = 1200.0, k = -9.8),
	(from = "med-1", to = "end-3", kind = :activation, at = 300.0, k = -7.0),
	(from = "med-2", to = "end-3", kind = :activation, at = 430.0, k = -7.3),
	(from = "pop-1.proteins*", to = "end-3", kind = :activation, at = 23000.0, k = -1.5),
	(from = "skn-1", to = "end-3", kind = :activation, at = 2500.0, k = -4.2),
	(from = "pop-1", to = "end-3", kind = :repression, at = 11000.0, k = -5.0),
	(from = "end-3", to = "end-1", kind = :activation, at = 610.0, k = -3.1),
	(from = "med-1", to = "end-1", kind = :activation, at = 1000.0, k = -4.0),
	(from = "pop-1.proteins*", to = "end-1", kind = :activation, at = 8200.0, k = -6.0),
	(from = "skn-1", to = "end-1", kind = :activation, at = 11000.0, k = -4.5),
	(from = "elt-2", to = "end-1", kind = :repression, at = 1000.0, k = -5.0),
	(from = "elt-7", to = "end-1", kind = :repression, at = 30000.0, k = -4.8),
	(from = "pop-1", to = "end-1", kind = :repression, at = 1500.0, k = -2.7),
	(from = "ref-1.proteins", to = "end-1", kind = :repression, at = 20000.0, k = -4.0),
	(from = "elt-7", to = "elt-7", kind = :activation, at = 20000.0, k = -4.0),
	(from = "end-1", to = "elt-7", kind = :activation, at = 4800.0, k = -4.2),
	(from = "end-3", to = "elt-7", kind = :activation, at = 700.0, k = -4.0),
	(from = "elt-2", to = "elt-7", kind = :repression, at = 1500.0, k = -4.0),
	(from = "elt-2", to = "elt-2", kind = :activation, at = 150.0, k = -1.5),
	(from = "elt-7", to = "elt-2", kind = :activation, at = 7000.0, k = -2.0),
	(from = "end-1", to = "elt-2", kind = :activation, at = 5000.0, k = -2.0),
	(from = "end-3", to = "elt-2", kind = :activation, at = 8800.0, k = -1.5),
	(from = "med-1", to = "tbx-35", kind = :activation, at = 660.0, k = -6.0),
	(from = "med-2", to = "tbx-35", kind = :activation, at = 2200.0, k = -4.9),
	(from = "ref-1.proteins*", to = "tbx-35", kind = :repression, at = 2000.0, k = -9.2),
	(from = "ceh-51", to = "ceh-51", kind = :activation, at = 300.0, k = -1.2),
	(from = "tbx-35", to = "ceh-51", kind = :activation, at = 3500.0, k = -9.0),
]

# ╔═╡ 6b607d66-4b66-42b4-b398-c881dfe33681
circuit_spec = let
	length = Dict(r.gene => r.length for r in CSV.File("data/gene_lengths.tsv"; stringtype = String))
	regulation(edges) = Dict(:slots => [Dict(:from => e.from, :at => e.at, :k => e.k, :w => Float64(e.at < 25000)) for e in edges], :aggregate => "generalized_mean", :p => -2)
	inputs(name, kind) = [e for e in edges if e.to == name && e.kind == kind]
	gene(name) = Dict(
		:name => name,
		:color => "#" * hex(gene_colours[name]),
		:unique => false,
		:base_rates => Dict(Symbol("\$") => "base_rates", :transcription => 40 / get(length, name, Inf),
		                    pairs(get(gene_specific_rates, name, (;)))...),
		(kind => regulation(inputs(name, kind)) for kind in (:activation, :repression) if !isempty(inputs(name, kind)))...,
	)
	#method = Dict(:name => "HybridTau", :policy => "CriticalBlend", :nc => 3, :epsilon => nothing, :dt => 10)
	method = "RSSACR"
	Dict(Symbol("{regulation/v1}") => Dict(
		:method => method,
		:genes => gene.(collect(keys(gene_colours))),
		:reactions => signalling,
	))
end

# ╔═╡ 4a343dc4-4aa4-4321-bcba-fe8c14de97f4
circuit = GRS.Models.parse(JSON.json(circuit_spec); base_rates);

# ╔═╡ 4ee67dea-8973-47fd-95d6-a7a4e39ec2d3
network_view = let view = CytoscapeJS.Cytoscape(GRS.Visualisation.Network(circuit); group_colors = gene_colours, height = "600px")
	view.stylesheet[] = [view.stylesheet[]; (; selector = "edge", style = Dict("label" => ""))]
	view
end

# ╔═╡ 98294067-5589-4cf6-9265-4b084a7e3af5
perturbable_spec = let spec = deepcopy(circuit_spec), m = spec[Symbol("{regulation/v1}")]
	m[:reactions] = [m[:reactions]; rnai]
	spec
end

# ╔═╡ 02fa3aff-7a91-4b5c-8390-73985d151495
schedule_spec(experiments...) = Dict(
	:base_rates => base_rates,
	:initial => Dict(Symbol("{add}") => merge(shared_machinery, maternal)),
	:P2_signalling => Dict(Symbol("{add}") => Dict(Symbol("P2-signal") => 3500)),
	:step => Dict(
		:EMS_circuit => perturbable_spec,
		:step => [Dict(Symbol("\$") => "initial"), experiments...],
	),
)

# ╔═╡ fd816447-d0ea-4579-8551-8325af1d90d9
atlas_schedule! = build(schedule_spec(atlas_spec))

# ╔═╡ 0d8be544-e359-4db6-8c06-4da097a23f2f
full_schedule! = build(schedule_spec(atlas_spec, knockout_spec))

# ╔═╡ 9f76477a-80a1-4d2d-8015-24f431d06155
samples = simulate(full_schedule!; seed = "1")

# ╔═╡ 366ad760-58b8-4919-894b-4f5e9b698cb4
simulated_transcriptome = transcriptome(samples, atlas_genes, atlas_stages)

# ╔═╡ 2de5d8fb-90af-412c-a4e6-2ae20705739b
plot_transcriptome(observed_transcriptome, simulated_transcriptome, atlas_genes)

# ╔═╡ c90c57f9-501e-447c-87f4-299733ba3d8a
plot_heatmap(observed_transcriptome, simulated_transcriptome)

# ╔═╡ 90582933-38a7-4bc1-8867-f79015f4035f
simulated_knockouts = knockout_proteins(samples, readouts, knockouts)

# ╔═╡ 732820a8-d0b1-4436-bc17-a97365de9241
plot_knockouts(observed_knockouts, simulated_knockouts, readouts)

# ╔═╡ 0ec32f1f-3614-4946-8672-30e070b26bd1
plot_knockouts(knockout_data, samples)

# ╔═╡ 5a60fd2f-250e-469c-a037-523339d6d4bd
show_schedule(schedule_spec(atlas_spec, knockout_spec); title = "celegans_ems.schedule.json")

# ╔═╡ 6f67a4e3-b5c3-4f1c-9aad-7e73dd862f24
write("schedules/celegans_ems.schedule.json", JSON.json(schedule_spec(atlas_spec, knockout_spec), 4));

# ╔═╡ f65b848e-e4e0-4079-8881-63f69d3d213d
let
	f = Figure(size = (1500, 1400))
	plot_stage_model!(f[1, 1:2])
	ax = Axis(f[2, 1]; aspect = DataAspect(), valign = :top)
	hidedecorations!(ax)
	hidespines!(ax)
	image!(ax, rotr90(Makie.FileIO.load("figures/network_species.png")))
	right = f[2, 2] = GridLayout()
	plot_heatmap!(right[1, 1], observed_transcriptome, simulated_transcriptome)
	plot_transcriptome!(right[2, 1], observed_transcriptome, simulated_transcriptome, atlas_genes)
	plot_knockouts!(right[3, 1], knockout_data, samples)
	for (pos, letter) in zip((f[1, 1:2], f[2, 1], right[1, 1], right[2, 1], right[3, 1]), "ABCDE")
		Label(pos[1, 1, TopLeft()], string(letter); font = :bold, fontsize = 22, padding = (0, 10, 10, 0), halign = :right)
	end
	colsize!(f.layout, 1, Relative(0.5))
	rowsize!(f.layout, 1, Fixed(340))
	rowsize!(right, 1, Relative(0.4))
	rowsize!(right, 2, Relative(0.34))
	save("figures/celegans_figure.pdf", f)
	f
end

# ╔═╡ Cell order:
# ╟─09291836-bddd-11f1-aaf5-cd7c7aaa7bc3
# ╟─afd5f0c5-eeb2-429b-848c-b055947be9cd
# ╠═8d6d12d2-97d0-42b2-8340-a8acf112e60f
# ╟─7a93e7dc-4b76-44c6-886a-24d5e8103f15
# ╟─29d8dde1-8e4b-4628-bdf6-87cfef5e0a4f
# ╟─ed2e5dc7-0489-47d1-adcd-450c239ea374
# ╠═11fe5a08-04f3-4073-b553-b98a9c980cbd
# ╠═28e41292-2cbc-4234-a4b3-5d4a2cacb1f4
# ╠═187cd4b3-2a64-4ecf-bcc6-24450b5107bb
# ╠═ebfbb616-18d1-4039-a193-b8c0c6383f5d
# ╟─4a0c8c51-0ac6-4bcb-8681-8b37aec85978
# ╟─693d942c-6270-4ec3-a519-b832600a8f86
# ╠═16292fb4-24f9-4041-949b-b4c58c70065b
# ╠═6b607d66-4b66-42b4-b398-c881dfe33681
# ╠═4a343dc4-4aa4-4321-bcba-fe8c14de97f4
# ╟─cab2d623-8e0d-415f-8192-7e14b22249c7
# ╠═4ee67dea-8973-47fd-95d6-a7a4e39ec2d3
# ╟─e1622993-d86a-4b44-9c25-1dc9a3cf7cc7
# ╟─361a46e8-cadc-456c-8aec-67ebd07205f9
# ╟─9a93134d-84cd-49f9-a337-cf688d491728
# ╟─f89e6f07-a453-4d56-9719-8ab4cb89c5f7
# ╟─a23dfd35-fcc4-4544-973e-8539c5725e33
# ╠═c5762b14-f2d4-4acd-aa10-0f8eda3a854b
# ╠═18571daf-95b1-4293-a9c7-6e77f91b87a0
# ╠═8015c70a-c431-4ae4-afb8-c49f04ae86a4
# ╟─2f30332d-debf-499a-979e-20a20cc47fae
# ╠═0652f10d-d855-4adb-958f-0d8716e2cea6
# ╠═6c2949ef-6309-45b1-a957-df3beae4a22e
# ╟─dc4d3bba-e07d-48d8-a45f-26e5c4e2c8a2
# ╟─87bf13f3-e71b-4a8d-826e-61bc1c635039
# ╟─616b8a62-6edc-4d6a-83c7-b2bffc9a368c
# ╠═7ce0250c-7a68-44dd-86b8-cfff20c53471
# ╠═ddbe9327-7b98-404f-8b79-49e1bd074c8b
# ╠═98294067-5589-4cf6-9265-4b084a7e3af5
# ╠═eabe7e00-bf3f-484f-b64d-980c2d05b4f1
# ╟─91114341-c844-45b9-8679-24079cc65ce3
# ╟─9770c432-429e-4321-baa2-06dbd7c8305b
# ╠═80458248-05b1-4380-a1ca-89125da04da0
# ╟─6cd49d0d-353d-4b5c-a1a2-cb1376ad7d06
# ╟─79679362-4810-4b7b-ba82-aee0801a1943
# ╠═02fa3aff-7a91-4b5c-8390-73985d151495
# ╠═3ef37b74-b704-49bb-b3fe-55baf15a6ffb
# ╠═fd816447-d0ea-4579-8551-8325af1d90d9
# ╠═0d8be544-e359-4db6-8c06-4da097a23f2f
# ╠═5a60fd2f-250e-469c-a037-523339d6d4bd
# ╠═6f67a4e3-b5c3-4f1c-9aad-7e73dd862f24
# ╟─6df7d3b4-adc1-4a38-996d-d460cd5881cc
# ╠═89cdf95b-2c77-47c5-8fd9-c5ad2b27f20a
# ╠═ca1373c7-44b9-48f0-ac6d-a2db7529970c
# ╠═9f76477a-80a1-4d2d-8015-24f431d06155
# ╟─f70bac03-3cb0-4cb2-af22-cbd52f9bd147
# ╟─21334557-28ff-490e-8fa8-f5f4a2abc320
# ╟─bbdfa49b-6bf5-42e3-87cd-20efe2e6b1a7
# ╟─d78cbb77-8286-47e1-9eeb-d0050ac8d61d
# ╟─70975d2d-ebd2-470b-885d-6c64c1483f36
# ╟─95509860-902a-404c-ba76-c7c81da94895
# ╠═ca8a141a-d0a3-4de0-8465-7a166823c82b
# ╠═ef3dc29c-a9d5-4147-9e9d-e6fd09146fc1
# ╠═366ad760-58b8-4919-894b-4f5e9b698cb4
# ╟─0c4ee0a8-276f-4677-bbad-ebbdbd5cd50d
# ╠═2de5d8fb-90af-412c-a4e6-2ae20705739b
# ╟─88610f5a-b196-459a-b5da-502d0bb37f7a
# ╠═8fb2685e-4db1-4d2e-8ff3-c92bdb173354
# ╠═5fdc5d3f-1830-4f9a-8073-5abb6fbd8ba8
# ╠═90582933-38a7-4bc1-8867-f79015f4035f
# ╠═2f952564-8523-44d6-b224-48aac07561f9
# ╠═732820a8-d0b1-4436-bc17-a97365de9241
# ╟─97d575cb-9efb-4cf6-8459-652ff9f348ee
# ╠═c90c57f9-501e-447c-87f4-299733ba3d8a
# ╟─61267d7f-a38d-4d31-a0af-0c5936c77021
# ╠═0ec32f1f-3614-4946-8672-30e070b26bd1
# ╠═115fad24-4ab7-4fc4-89ab-f6a32d99296a
# ╠═f65b848e-e4e0-4079-8881-63f69d3d213d
