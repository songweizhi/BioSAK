install.packages("riverplot")
library(riverplot)
nodes <- c( LETTERS[1:3] )
edges <- list( A= list( C= 10 ), B= list( C= 10 ) )
r <- makeRiver( nodes, edges, node_xpos= c( 1,1,2 ),
                node_labels= c( A= "Node A", B= "Node B", C= "Node C" ),
                node_styles= list( A= list( col= "yellow" )) )
plot( r )
# equivalent form:
nodes <- data.frame( ID= LETTERS[1:3],
                     x= c( 1, 1, 2 ),
                     col= c( "yellow", NA, NA ),
                     labels= c( "Node A", "Node B", "Node C" ),
                     stringsAsFactors= FALSE )
r <- makeRiver( nodes, edges )
plot( r )
# all nodes but "A" will be red:
r <- makeRiver( nodes, edges, default_style= list( col="red" ) )
plot( r )
# overwrite the node information from "nodes":
r <- makeRiver( nodes, edges, node_styles= list( A=list( col="red" ) ) )
plot( r )



data( minard )
minard
nodes <- minard$nodes
edges <- minard$edges
colnames( nodes ) <- c( "ID", "x", "y" )
colnames( edges ) <- c( "N1", "N2", "Value", "direction" )
# color the edges by troop movement direction
edges$col <- c( "#e5cbaa", "black" )[ factor( edges$direction ) ]
# color edges by their color rather than by gradient between the nodes
edges$edgecol <- "col"
# generate the riverplot object and a style
river <- makeRiver( nodes, edges )
style <- list( edgestyle= "straight", nodestyle= "invisible" )
# plot the generated object
plot( river, lty= 1, default_style= style )
# Add cities
with( minard$cities, points( Longitude, Latitude, pch= 19 ) )
with( minard$cities, text( Longitude, Latitude, Name, adj= c( 0, 0 ) ) )


x <- riverplot.example() 
plot(x) 
plot(x, srt=90, lty=1)




nodes <- c( LETTERS[1:3] ) 
edges <- list( A= list( C= 10 ), B= list( C= 10 ) ) 
r <- makeRiver( nodes, edges, node_xpos= c( 1,1,2 ), node_labels= c( A= "Node A", B= "Node B", C= "Node C" ), node_styles= list( A= list( col= "yellow" )) ) 
plot( r )

# equivalent form: 
nodes <- data.frame( ID= LETTERS[1:3], x= c( 1, 1, 2 ), col= c( "yellow", NA, NA ), labels= c( "Node A", "Node B", "Node C" ), stringsAsFactors= FALSE ) 
r <- makeRiver( nodes, edges ) 
plot( r ) 
# all nodes but "A" will be red: 
r <- makeRiver( nodes, edges, default_style= list( col="red" ) ) 
plot( r ) 
# overwrite the node information from "nodes": 
r <- makeRiver( nodes, edges, node_styles= list( A=list( col="red" ) ) ) 
plot( r )





# To view the default style specification, type
default.style()
ex <- riverplot.example()
ex
ds <- default.style()
ds
plot( ex, default_style= ds )
# nodes with unspecified style will now be semi-transparent red:
ds[["col"]] <- "#FF000099"
plot( ex, default_style= ds )

