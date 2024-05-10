library(bio3d)

pdb <- read.pdb('./1a1u.pdb')

chainA <- atom.select(pdb, 'protein', chain = 'A')
chainB <- atom.select(pdb, 'protein', chain = 'C')

write.pdb(trim.pdb(pdb, chainA), './1a1u_A.pdb')
write.pdb(trim.pdb(pdb, chainB), './1a1u_C.pdb')
