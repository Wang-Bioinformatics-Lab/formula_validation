#install.packages("devtools")
#devtools::install_github("mjhelf/MassTools")

library(MassTools)

mz = 352.22537
charge = 0
ppm = 10
min_carbons = 1
max_carbons = 24
min_hydrogens = 1
max_hydrogens = 36
min_oxygen = 1
max_oxygen = 6
max_sulfur = 1
max_nitrogen = 1
top = 3

minElements = paste0("C", min_carbons, "H", min_hydrogens, "O", min_oxygen)
maxElements = paste0("C", max_carbons, "H", max_hydrogens, "O", max_oxygen, "N", max_nitrogen, "S", max_sulfur)

elements_list = c("C", "H", "N", "O", "S")


# Create the R code to initialize the 'elements' list
elements = Rdisop::initializeElements(elements_list)


filters = list(
    minElements = minElements,
    maxElements = maxElements,
    maxCounts = TRUE,
    HCratio = FALSE,
    moreRatios = TRUE
)

mfs = calcMF(mz = mz, z = charge, ppm = ppm, top = top, elements = elements, Filters = filters)

print(mfs)



input_data = read.csv("/home/alberto/Escritorio/monoisotopic_mass.csv")


results = list()


for (i in seq_len(nrow(input_data))) {
  monoisotopic_mass = input_data$monoisotopic_mass[i]

  mf_result <- tryCatch({
    calcMF(mz = monoisotopic_mass, z = charge, ppm = ppm, top = top,
           elements = elements, Filters = filters)
  }, error = function(e) {
    message(sprintf("Error at row %d: %s", i, e$message))
    return(NULL)
  })

  if (!is.null(mf_result)) {
    row <- list(mz = monoisotopic_mass)

    for (j in seq_len(top)) {
        if (j <= nrow(mf_result)) {
            row[[paste0("MF_", j)]]     <- mf_result$MF[j]
            row[[paste0("Score_", j)]]  <- mf_result$RdisopScore[j]
            row[[paste0("NRule_", j)]]  <- mf_result$nrule[j]
            row[[paste0("ppm_", j)]]    <- mf_result$ppm[j]
        } else {
            row[[paste0("MF_", j)]]     <- NA
            row[[paste0("Score_", j)]]  <- NA
            row[[paste0("NRule_", j)]]  <- NA
            row[[paste0("ppm_", j)]]    <- NA
        }
    }
  }
  else {
    # Fill in NAs for expected columns if no results
    for (j in 1:top) {
      row[[paste0("MF_", j)]]     <- NA
      row[[paste0("Score_", j)]]  <- NA
      row[[paste0("NRule_", j)]]  <- NA
      row[[paste0("ppm_", j)]]    <- NA
    }
  }
  results[[length(results) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
}

# Combine all into a data frame
output_df <- do.call(rbind, results)

# Write to CSV
write.csv(output_df, "/home/alberto/Escritorio/annotated_formulas.csv", row.names = FALSE)

cat("✅ Results written to 'annotated_formulas.csv'\n")