install.packages("networkD3")
install.packages('curl')
library(networkD3)
library(curl)



URL <- paste0("https://cdn.rawgit.com/christophergandrud/networkD3/", "master/JSONdata/energy.json") 
energy <- jsonlite::fromJSON(URL) 
energy
sankeyNetwork(Links = energy$links, Nodes = energy$nodes, Source = "source", Target = "target", Value = "value", NodeID = "name", units = "TWh", fontSize = 12, nodeWidth = 30)
energy$links$energy_type <- sub(" .*", "", energy$nodes[energy$links$source + 1, "name"]) sankeyNetwork(Links = energy$links, Nodes = energy$nodes, Source = "source", Target = "target", Value = "value", NodeID = "name", LinkGroup = "energy_type", NodeGroup = NULL )


