library(circlize)
mat = matrix(1:9, 3)
mat

rownames(mat) = letters[1:3]
colnames(mat) = LETTERS[1:3]
df = data.frame(from = letters[1:3], to = LETTERS[1:3], value = 1:3)
df

chordDiagram(mat)
circos.clear()

chordDiagram(df) 
circos.clear()
