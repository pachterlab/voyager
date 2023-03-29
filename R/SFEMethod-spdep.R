#' @include res2df.R

spdep_uni <- c(package = "spdep", variate = "uni")
uni_global <- c(spdep_uni, scope = "global", default_attr = NA)
uni_local <- c(spdep_uni, scope = "local")

.to_df_identity <- function(out, nb, p.adjust.method) out

# Construct SFEMethod object for Moran's I
moran <- SFEMethod(
    c(name = "moran", title = "Moran's I", uni_global),
    fun = function(x, listw, zero.policy = NULL)
        spdep::moran(x, listw, n = length(listw$neighbours), S0 = spdep::Szero(listw),
                     zero.policy = zero.policy),
    reorganize_fun = .moran2df
)

geary <- SFEMethod(
    c(name = "geary", title = "Geary's C", uni_global),
    fun = function(x, listw, zero.policy = NULL)
        spdep::geary(x, listw,
                     n = length(listw$neighbours),
                     n1 = length(listw$neighbours) - 1,
                     S0 = spdep::Szero(listw),
                     zero.policy = zero.policy),
    reorganize_fun = .moran2df
)

moran.mc <- SFEMethod(
    c(name = "moran.mc", title = "Moran's I with permutation testing",
      uni_global),
    fun = spdep::moran.mc,
    reorganize_fun = .mcsim2df,
    args_not_check = "nsim"
)

geary.mc <- SFEMethod(
    c(name = "geary.mc", title = "Geary's C with permutation testing",
      uni_global),
    fun = spdep::geary.mc,
    reorganize_fun = .mcsim2df,
    args_not_check = "nsim"
)

sp.mantel.mc <- SFEMethod(
    c(name = "spmantel.mc", title = "Mantel-Hubert spatial general cross product statistic",
      uni_global),
    fun = function(x, listw, ..., zero.policy = NULL)
        spdep::sp.mantel.mc(x, listw, ..., zero.policy = zero.policy),
    reorganize_fun = .mcsim2df,
    args_not_check = "nsim"
)

moran.test <- SFEMethod(
    c(name = "moran.test", title = "Moran's I test", uni_global),
    fun = spdep::moran.test,
    reorganize_fun = .htest2df
)

geary.test <- SFEMethod(
    c(name = "geary.test", title = "Geary's C test", uni_global),
    fun = spdep::geary.test,
    reorganize_fun = .htest2df
)

globalG.test <- SFEMethod(
    c(name = "globalG.test", title = "Global G test", uni_global),
    fun = spdep::globalG.test,
    reorganize_fun = .htest2df
)

.sp.correlogram <- function(x, listw, method = "I", ..., zero.policy = NULL) {
    spdep::sp.correlogram(neighbours = listw$neighbours, var = x,
                          method = method, ..., zero.policy = zero.policy)
}

sp.correlogram <- SFEMethod(
    c(name = "sp.correlogram", title = "Correlogram", uni_global),
    fun = .sp.correlogram,
    reorganize_fun = .correlogram2df,
    args_not_check = c("order", "method")
)

localmoran <- SFEMethod(
    c(name = "localmoran", title = "Local Moran's I", uni_local,
      default_attr = "Ii"),
    fun = spdep::localmoran,
    reorganize_fun = .localmoran2df
)

localmoran_perm <- SFEMethod(
    c(name = "localmoran_perm", title = "Local Moran's I permutation testing",
      uni_local, default_attr = "Ii"),
    fun = spdep::localmoran_perm,
    reorganize_fun = .localmoran2df,
    args_not_check = "nsim"
)

localC <- SFEMethod(
    c(name = "localC", title = "Local Geary's C", uni_local,
      default_attr = "localC"),
    fun = spdep:::localC.default,
    reorganize_fun = .to_df_identity
)

localC_perm <- SFEMethod(
    c(name = "localC_perm", title = "Local Geary's C permutation testing",
      uni_local, default_attr = "localC"),
    fun = spdep:::localC_perm.default,
    reorganize_fun = .localCperm2df,
    args_not_check = "nsim"
)

localG <- SFEMethod(
    c(name = "localG", title = "Getis-Ord Gi(*)", uni_local,
      default_attr = "localG"),
    fun = spdep::localG,
    reorganize_fun = .localG2df
)

localG_perm <- SFEMethod(
    c(name = "localG_perm", title = "Getis-Ord Gi(*) with permutation testing",
      uni_local, default_attr = "localG"),
    fun = spdep::localG_perm,
    reorganize_fun = .localG2df,
    args_not_check = "nsim"
)

LOSH <- SFEMethod(
    c(name = "LOSH", title = "Local spatial heteroscedasticity", uni_local,
      default_attr = "Hi"),
    fun = spdep::LOSH,
    reorganize_fun = .to_df_identity
)

LOSH.mc <- SFEMethod(
    c(name = "LOSH.mc", title = "Local spatial heteroscedasticity permutation testing",
      default_attr = "Hi", uni_local),
    fun = spdep::LOSH.mc,
    reorganize_fun = .LOSHmc2df,
    args_not_check = "nsim"
)

LOSH.cs <- SFEMethod(
    c(name = "LOSH.cs", title = "Local spatial heteroscedasticity Chi-square test",
      uni_local, default_attr = "Hi"),
    fun = spdep::LOSH.cs,
    reorganize_fun = .LOSHmc2df
)

moran.plot <- SFEMethod(
    c(name = "moran.plot", title = "Moran plot", uni_local,
      default_attr = "wx"),
    fun = function(x, listw, ..., zero.policy = NULL)
        spdep::moran.plot(x, listw, zero.policy = zero.policy,
                          plot = FALSE, return_df = TRUE),
    reorganize_fun = .to_df_identity
)

# Multivariate
#localC_multi <- SFEMethod(
#    c(name = "localC_multi", title = "Multivariate local Geary's C",
#      package = "spdep", variate = "multi", scope = "local",
#      default_attr = "localC"),
#    fun = spdep::localC, # need to write my own wrapper
#    reorganize_fun = .to_df_identity
#)
