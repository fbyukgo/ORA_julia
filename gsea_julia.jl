using JSON
include("hypergeometric.jl")

#reading the pathway json file
function read_json(file)
    open(file,"r") do f
        return JSON.parse(f)
    end
end

##
all_dict = read_json("kegg_medicus.json")
nPathways = length(keys(all_dict))
pathwayGenelist = Dict()
for (key, val) in all_dict
    pathwayGenelist[key] = val["geneSymbols"]
end

#=
myGeneList = ["CYCS", "VDAC3"]
example_list = pathwayGenelist["KEGG_MEDICUS_VARIANT_MUTATION_INACTIVATED_SIGMAR1_TO_CA2_APOPTOTIC_PATHWAY"]
intersect(myGeneList, example_list)
=#

function main()
    #run stuff
end
