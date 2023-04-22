sfe_uni_global <- data.frame(
    name = c("moran", "geary", "moran.mc", "geary.mc", "sp.mantel.mc",
             "moran.test", "geary.test", "globalG.test", "sp.correlogram",
             "variogram", "variogram_map"),
    description = c("Moran's I", "Geary's C", "Moran's I with permutation testing",
                    "Geary's C with permutation testing",
                    "Mantel-Hubert spatial general cross product statistic",
                    "Moran's I test", "Geary's C test", "Global G test",
                    "Correlogram", "Variogram with model", "Variogram map")
)
sfe_uni_global$variate <- "uni"
sfe_uni_global$scope <- "global"

sfe_uni_local <- data.frame(
    name = c("localmoran", "localmoran_perm", "localC", "localC_perm", "localG",
             "localG_perm", "LOSH", "LOSH.mc", "LOSH.cs", "moran.plot"),
    description = c("Local Moran's I", "Local Moran's I permutation testing",
                    "Local Geary's C", "Local Geary's C permutation testing",
                    "Getis-Ord Gi(*)", "Getis-Ord Gi(*) with permutation testing",
                    "Local spatial heteroscedasticity",
                    "Local spatial heteroscedasticity permutation testing",
                    "Local spatial heteroscedasticity Chi-square test",
                    "Moran scatter plot")
)
sfe_uni_local$variate <- "uni"
sfe_uni_local$scope <- "local"

sfe_bi_global <- data.frame(
    name = c("lee", "lee.mc", "lee.test", "cross_variogram",
             "cross_variogram_map"),
    description = c("Lee's bivariate statistic",
                    "Lee's bivariate static with permutation testing",
                    "Lee's L test", "Cross variogram",
                    "Cross variogram map")
)
sfe_bi_global$variate <- "bi"
sfe_bi_global$scope <- "global"

sfe_bi_local <- data.frame(
    name = c("locallee", "localmoran_bv"),
    description = c("Local Lee's bivariate statistic",
                    "Local bivariate Moran's I")
)
sfe_bi_local$variate <- "bi"
sfe_bi_local$scope <- "local"

sfe_multi <- data.frame(
    name = c("multispati", "localC_multi", "localC_perm_multi"),
    description = c("MULTISPATI PCA", "Multivariate local Geary's C",
                    "Multivariate local Geary's C permutation testing")
)
sfe_multi$variate <- "multi"
sfe_multi$scope <- NA

sfe_methods <- do.call(rbind, list(sfe_uni_global, sfe_uni_local,
                                   sfe_bi_global, sfe_bi_local,
                                   sfe_multi))
# read existing sysdata
load("./R/sysdata.rda")
usethis::use_data(sfe_methods, ditto_colors, internal = TRUE, overwrite = TRUE)
