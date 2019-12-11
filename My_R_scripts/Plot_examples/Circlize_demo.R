library(circlize)

# mat = read.table('/Users/songweizhi/Dropbox/Research/MetaCHIP_manuscript/real_dataset/26bins_HGT_direction.txt', header = TRUE)

mat = read.table('/Users/songweizhi/Desktop/HGT_candidates_ET_validated_matrix.csv', header = TRUE)


grid.col = c(A = 'brown1', B = 'lawngreen', C = 'mediumorchid', D = 'mediumslateblue', E = 'royalblue', F = 'sandybrown')

chordDiagram(t(mat))
#chordDiagram(t(mat), grid.col = grid.col)

par(mar = rep(0,4), cex = 1.3)    

#for(si in get.all.sector.index()) {
#  circos.axis(h = "top", labels.cex = 0.3, sector.index = si, track.index = 2)
#}

circos.clear()

