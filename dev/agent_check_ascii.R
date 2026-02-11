# Script pour identifier les caractères non-ASCII dans les fichiers R
library(tools)

files_to_check <- c(
  "C:/Users/masse/Desktop/KefiR/KefiR/R/cooks.distance_lmer.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/evolreg.r",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/gage_rr.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/m.test.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/pairwise.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/pairwise.boot.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/pde.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sample_bootstrap.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_analyze_ancova_structure.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_ancova_analysis.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_auto_preprocess_g.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_control_independence.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_formulator.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_multi_factor_analysis.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_one_factor_analysis.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_pairwise_medpb2.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_posthoc.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_posthoc_mixed_model.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_select_ss_type.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/sys_variance.R",
  "C:/Users/masse/Desktop/KefiR/KefiR/R/valreg.R"
)

cat("Checking for non-ASCII characters in files...\n\n")

for (file in files_to_check) {
  if (file.exists(file)) {
    result <- showNonASCIIfile(file)
    if (length(result) > 0) {
      cat("File:", basename(file), "\n")
      print(result)
      cat("\n")
    }
  }
}

cat("Check complete.\n")
