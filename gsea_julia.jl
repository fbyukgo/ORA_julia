using JSON

#reading the pathway json file
function read_json(file)
    open(file,"r") do f
        return JSON.parse(f)
    end
end



all_dict = read_json("kegg_medicus.json")
length(keys(kegg_dict))

myGeneList = ["CYCS", "VDAC3"]
example_list = pathwayGenelist["KEGG_MEDICUS_VARIANT_MUTATION_INACTIVATED_SIGMAR1_TO_CA2_APOPTOTIC_PATHWAY"]

intersect(myGeneList, example_list)

pathwayGenelist = Dict()
for (key, val) in kegg_dict
    pathwayGenelist[key] = val["geneSymbols"]
end

