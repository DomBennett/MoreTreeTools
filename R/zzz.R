# ---------------
# Notes on coding
# ---------------
#
# Vectorized for-loop with plyr:
# -- use mdply if returning a dataframe
# -- use m_ply if returning multiple things (use <<- to write to parent environment)
# -- looping code is contained within .fun()
#
# NodeList and phylo:
# -- Code is developed for both classes
# -- Functions that work on both have visible functions that point to hidden functions specific
#  for each class type
# -- See *-phylo and *-nodelist files for hidden functions
# -- See *-visible for public functions
#
# --------------
# Notes on style
# --------------
#
# Variable names:
#  1. lower case
#  2. separated by '_'
#  e.g. tip_labels
#
# Function/method names:
#  1. camel case
#  2. first letter is always lower case
#  getChildren()
#
# Class names:
#  1. camel case
#  2. first letter is always upper case
#  e.g. NodeList
#
# Hidden functions:
#  1. camel case
#  2. __[asscociation], being a class name or a function
#  e.g. .getChildren__NodeList()
#   .sample__mapNames()