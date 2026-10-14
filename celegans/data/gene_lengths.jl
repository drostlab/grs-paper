using JSON, CSV, DataFrames

# gene spans on C. elegans assembly WBcel235
# Dyer SC, Austine-Orimoloye O, Azov AG, et al. (2024).
# Ensembl 2025.
# Nucleic Acids Res 53(D1), D948–D957. doi:10.1093/nar/gkae1071

genes = ["skn-1", "pop-1", "med-1", "med-2", "end-1", "end-3", "elt-7", "elt-2", "tbx-35", "ceh-51"]

lookup(gene) = JSON.parsefile(download("https://rest.ensembl.org/lookup/symbol/caenorhabditis_elegans/$gene?content-type=application/json"))

lengths = DataFrame(map(genes) do gene
    d = lookup(gene)
    (; gene, id=d["id"], chromosome=d["seq_region_name"], start=d["start"],
    stop=d["end"], length=d["end"]-d["start"]+1, assembly=d["assembly_name"])
end)

CSV.write("data/gene_lengths.tsv", lengths; delim='\t')
