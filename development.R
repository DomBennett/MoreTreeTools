## D.J Bennett
## Testing out new ideas

library (MoreTreeTools)
## Plotting phylogenetic signal
tree <- compute.brlen (rtree (100))
trait <- randCommData (tree, nsites = 1, nspp = 25)[1, ]
blockplot (tree, trait, title = 'Random')
etrait <- evenCommData (tree, nsites = 1, nspp = 25)[1, ]
blockplot (tree, etrait, title = 'Even')

for (psi in c (-100, -10, 1, 10, 100)) {
  ctrait <- genCommData (tree, focal = 25,
                         psi = psi, mean.incid = 25)[1, ]
  blockplot (tree, ctrait, title = paste0 ('Clustered, psi = [',
                                           psi, ']'))
}
