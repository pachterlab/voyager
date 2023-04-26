#' @include res2df.R

# Construct SFEMethod object for Moran's I
moran <- SFEMethod(
    name = "moran", title = "Moran's I", package = "spdep",
    fun = function(x, listw, zero.policy = NULL)
        spdep::moran(x, listw, n = length(listw$neighbours), S0 = spdep::Szero(listw),
                     zero.policy = zero.policy),
    reorganize_fun = .moran2df
)

geary <- SFEMethod(
    name = "geary", title = "Geary's C", package = "spdep",
    fun = function(x, listw, zero.policy = NULL)
        spdep::geary(x, listw,
                     n = length(listw$neighbours),
                     n1 = length(listw$neighbours) - 1,
                     S0 = spdep::Szero(listw),
                     zero.policy = zero.policy),
    reorganize_fun = .moran2df
)

moran.mc <- SFEMethod(
    name = "moran.mc", title = "Moran's I with permutation testing",
    package = "spdep",
    fun = spdep::moran.mc,
    reorganize_fun = .mcsim2df,
    args_not_check = "nsim"
)

geary.mc <- SFEMethod(
    name = "geary.mc", title = "Geary's C with permutation testing",
    package = "spdep",
    fun = spdep::geary.mc,
    reorganize_fun = .mcsim2df,
    args_not_check = "nsim"
)

sp.mantel.mc <- SFEMethod(
    name = "spmantel.mc", title = "Mantel-Hubert spatial general cross product statistic",
    package = "spdep",
    fun = function(x, listw, ..., zero.policy = NULL)
        spdep::sp.mantel.mc(x, listw, ..., zero.policy = zero.policy),
    reorganize_fun = .mcsim2df,
    args_not_check = "nsim"
)

moran.test <- SFEMethod(
    name = "moran.test", title = "Moran's I test", package = "spdep",
    fun = spdep::moran.test,
    reorganize_fun = .htest2df
)

geary.test <- SFEMethod(
    name = "geary.test", title = "Geary's C test", package = "spdep",
    fun = spdep::geary.test,
    reorganize_fun = .htest2df
)

globalG.test <- SFEMethod(
    name = "globalG.test", title = "Global G test", package = "spdep",
    fun = spdep::globalG.test,
    reorganize_fun = .htest2df
)

.sp.correlogram <- function(x, listw, method = "I", ..., zero.policy = NULL) {
    spdep::sp.correlogram(neighbours = listw$neighbours, var = x,
                          method = method, ..., zero.policy = zero.policy)
}

sp.correlogram <- SFEMethod(
    name = "sp.correlogram", title = "Correlogram", package = "spdep",
    fun = .sp.correlogram,
    reorganize_fun = .correlogram2df,
    args_not_check = c("order", "method")
)

localmoran <- SFEMethod(
    name = "localmoran", title = "Local Moran's I",
    package = "spdep", scope = "local", default_attr = "Ii",
    fun = spdep::localmoran,
    reorganize_fun = .localmoran2df
)

localmoran_perm <- SFEMethod(
    name = "localmoran_perm", title = "Local Moran's I permutation testing",
    package = "spdep", scope = "local", default_attr = "Ii",
    fun = spdep::localmoran_perm,
    reorganize_fun = .localmoran2df,
    args_not_check = "nsim"
)

localC <- SFEMethod(
    name = "localC", title = "Local Geary's C",
    package = "spdep", scope = "local", default_attr = "localC",
    fun = spdep:::localC.default,
    reorganize_fun = .to_df_identity
)

localC_perm <- SFEMethod(
    name = "localC_perm", title = "Local Geary's C permutation testing",
    package = "spdep", scope = "local", default_attr = "localC",
    fun = spdep:::localC_perm.default,
    reorganize_fun = .localCperm2df,
    args_not_check = "nsim"
)

localG <- SFEMethod(
    name = "localG", title = "Getis-Ord Gi(*)",
    package = "spdep", scope = "local", default_attr = "localG",
    fun = spdep::localG,
    reorganize_fun = .localG2df
)

localG_perm <- SFEMethod(
    name = "localG_perm", title = "Getis-Ord Gi(*) with permutation testing",
    package = "spdep", scope = "local", default_attr = "localG",
    fun = spdep::localG_perm,
    reorganize_fun = .localG2df,
    args_not_check = "nsim"
)

LOSH <- SFEMethod(
    name = "LOSH", title = "Local spatial heteroscedasticity",
    package = "spdep", scope = "local", default_attr = "Hi",
    fun = spdep::LOSH,
    reorganize_fun = .to_df_identity
)

LOSH.mc <- SFEMethod(
    name = "LOSH.mc", title = "Local spatial heteroscedasticity permutation testing",
    package = "spdep", scope = "local", default_attr = "Hi",
    fun = spdep::LOSH.mc,
    reorganize_fun = .LOSHmc2df,
    args_not_check = "nsim"
)

LOSH.cs <- SFEMethod(
    name = "LOSH.cs", title = "Local spatial heteroscedasticity Chi-square test",
    package = "spdep", scope = "local", default_attr = "Hi",
    fun = spdep::LOSH.cs,
    reorganize_fun = .LOSHmc2df
)

moran.plot <- SFEMethod(
    name = "moran.plot", title = "Moran plot", package = "spdep", scope = "local",
    default_attr = "wx",
    fun = function(x, listw, ..., zero.policy = NULL)
        spdep::moran.plot(x, listw, zero.policy = zero.policy,
                          plot = FALSE, return_df = TRUE),
    reorganize_fun = .to_df_identity
)
