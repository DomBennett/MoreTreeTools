# Dom Bennett
# Simulate three fossil trees for plotting

# LIB
library (MoreTreeTools)

# PROCESS
cat ('\nDE tree ....')
de <- runEDBMM (birth=2, death=1, sig=-1, eps=1, stop.at=500,
                fossils=TRUE)
ed.de <- calcEdgeDiversity (de, n.intervals=4)
ed.de$col <- (log (ed.de$count) - mean (log (ed.de$count))) /
  sd (log (ed.de$count))
# Pan tree
cat ('\nPan tree ....')
pan <- runEDBMM (birth=2, death=1, sig=-1, eps=-1, stop.at=500,
                 fossils=TRUE)
ed.pan <- calcEdgeDiversity (pan, n.intervals=4)
ed.pan$col <- (log (ed.pan$count) - mean (log (ed.pan$count))) /
  sd (log (ed.pan$count))
# hydra tree
cat ('\nHyd tree ....')
hyd <- runEDBMM (birth=2, death=1, sig=-1, eps=0, stop.at=500,
                 fossils=TRUE)
ed.hyd <- calcEdgeDiversity (hyd, n.intervals=4)
ed.hyd$col <- (log (ed.hyd$count) - mean (log (ed.hyd$count))) /
  sd (log (ed.hyd$count))
cat ('\nPacking ....')
edbmm.trees <- list ("DE" = list ('tree' = de, 'edge.diversity' = ed.de),
                     "Pan" = list ('tree' = pan, 'edge.diversity' = ed.pan),
                     "Hyd" = list ('tree' = hyd, 'edge.diversity' = ed.hyd))
cat ('\nOutput ....')
save (edbmm.trees, file="EDBMMtrees.Rd")
cat ('\nDone.')