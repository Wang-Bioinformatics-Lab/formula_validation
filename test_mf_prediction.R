install.packages("devtools")
devtools::install_github("mjhelf/MassTools")

library(MassTools)


mz = 200.000659
charge = 1
ppm = 50
min_carbons = 8
max_carbons = 10

min_hydrogens = 4
max_hydrogens = 10
minElements = paste0("C", min_carbons, "H", min_hydrogens, "S0")
maxElements = paste0("C", max_carbons, "H", max_hydrogens, "S0")
top = 3

# FULL FILTER
# Filters = list(DBErange = c(-5,40),
#                                   minElements = "C0",
#                                   maxElements = "C99999",
#                                   parity = "e",
#                                   maxCounts = TRUE,
#                                   SENIOR3 = 0,
#                                   HCratio = TRUE,
#                                   moreRatios = TRUE,
#                                   elementHeuristic = TRUE)

filters = list(minElements = minElements,
                maxElements = maxElements,
                maxCounts = TRUE,
                HCratio = FALSE,
                moreRatios = FALSE)

elements = Rdisop::initializeElements(c("C", "H", "N", "O", "S"))
mfs = calcMF(mz = mz, z = charge, ppm = ppm, top = top, elements = elements, Filters = filters)
print(mfs)

# BIG NOTE TO DO ADD DEPENDENCY R-BASE TO THE CONDA ENVIRONMENT - r-base and rpy in pip dependencies