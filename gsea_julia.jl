using JSON
include("hypergeometric.jl")

#reading the pathway json file
function read_json(file)
    open(file,"r") do f
        return JSON.parse(f)
    end
end

#load the dictionary from a json file 
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
json_dict = read_json("kegg_medicus.json")
nPathways = length(keys(json_dict))
pathwayDict = Dict()
for (key, val) in json_dict
    pathwayDict[key] = val["geneSymbols"]
end



backgroundGenes = []

for (key, val) in pathwayDict
    for gene in val
        if !(gene in backgroundGenes)
            push!(backgroundGenes, val)
        end
    end    
end

example_list = pathwayDict["KEGG_MEDICUS_VARIANT_MUTATION_INACTIVATED_SIGMAR1_TO_CA2_APOPTOTIC_PATHWAY"]


function run_ora(backgroundGenes, geneList, pathwayDict)

    resultDict = Dict()
    alpha = 0.05
    numPathways = length(pathwayDict)
    pvalThr = alpha / numPathways 

    lenList = length(geneList)
    lenBG = length(backgroundGenes)
    for (key, val) in pathwayDict
        pathLength = length(val)
        intersectGenes = intersect(geneList, val)
        interLen = length(intersect(geneList, val))
        pval = hypergeometric(interLen, lenBG, lenList, pathLength)
        if pval <= pvalThr
            resultDict[key] = Dict()
            resultDict[key]["pval_adj"] = pval
            resultDict[key]["no_intersection"] = interLen
            resultDict[key]["intersection"] = intersectGenes
        end
    end 
    return resultDict 
end

results = run_ora(backgroundGenes, example_list, pathwayDict )


results["KEGG_MEDICUS_VARIANT_MUTATION_INACTIVATED_SIGMAR1_TO_CA2_APOPTOTIC_PATHWAY"]
for (key, val) in results
    println(key, val)
end






function main()
    #run stuff
end
