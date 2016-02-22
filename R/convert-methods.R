setOldClass ('phylo')

setAs(from="TreeMan", to="phylo", def=function(from, to) {
  treeman::writeTree(from, file='temp.tre')
  tree <- ape::read.tree(file='temp.tre')
  file.remove('temp.tre')
  return(tree)
})

setAs(from="phylo", to="TreeMan", def=function(from, to) {
  ape::write.tree(from, file='temp.tre')
  tree <- treeman::readTree(file='temp.tre')
  file.remove('temp.tre')
  return(tree)
})