# 09_CT_total_opacity.R -- CT total opacity (GMM GGO + Consolidation + Dense Tissue = 100 - aerated lung) per patient at baseline and
#   follow-up and its change; read by 07 (Figure 4j, 5h).
# Produces: opacity/opacity_per_patient.csv.
# Inputs: planA_data_RA24.csv.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
R0 <- file.path(ART_TABLES, "n31"); OD <- file.path(R0, "opacity"); dir.create(OD, recursive = TRUE, showWarnings = FALSE)
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"))
ra$T1_opa <- 100 - ra$T1_aer; ra$T2_opa <- 100 - ra$T2_aer; ra$d_opa <- -ra$d_aer
stopifnot(max(abs(ra$T1_opa - (ra$CT_T1_GGO_pct + ra$CT_T1_Consolidation_pct + ra$CT_T1_Dense_Tissue_pct))) < 0.05,
          max(abs(ra$d_opa - ra$CT_Delta_lung_involvement_pct)) < 0.05)
write.csv(ra[, c("Sample_ID", "T1_opa", "T2_opa", "d_opa")], file.path(OD, "opacity_per_patient.csv"), row.names = FALSE)
