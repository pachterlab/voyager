#' @include res2df.R

spdep_uni <- c(package = "spdep", variate = "uni")
uni_global <- c(spdep_uni, scope = "global", default_attr = NA)
uni_local <- c(spdep_uni, scope = "local")

# Construct SFEMethod object for Moran's I
moran <- SFEMethod(
    c(name = "moran", title = "Moran's I", uni_global),
    fun = function(x, listw) spdep::moran(x, listw, n = length(listw$neighbours), S0 = Szero(listw)),
    to_df_fun = .moran2df
)

geary <- SFEMethod(
    c(name = "geary", title = "Geary's C", uni_global),
    fun = function(x, listw)
        spdep::geary(x, listw,
                     n = length(listw$neighbours),
                     n1 = length(listw$neighbours) - 1,
                     S0 = Szero(listw)),
    to_df_fun = .moran2df
)

moran.mc <- SFEMethod(
    c(name = "moran.mc", title = "Moran's I with permutation testing",
      uni_global),
    fun = spdep::moran.mc,
    to_df_fun = .mcsim2df
)

geary.mc <- SFEMethod(
    c(name = "geary.mc", title = "Geary's C with permutation testing",
      uni_global),
    fun = spdep::geary.mc,
    to_df_fun = .mcsim2df
)

sp.mantel.mc <- SFEMethod(
    c(name = "spmantel.mc", title = "Mantel-Hubert spatial general cross product statistic",
      uni_global),
    fun = function(x, listw, ...) spdep::sp.mantel.mc(x, listw, ...),
    to_df_fun = .mcsim2df
)

moran.test <- SFEMethod(
    c(name = "moran.test", title = "Moran's I test", uni_global),
    fun = spdep::moran.test,
    to_df_fun = .htest2df
)

geary.test <- SFEMethod(
    c(name = "geary.test", title = "Geary's C test", uni_global),
    fun = spdep::geary.test,
    to_df_fun = .htest2df
)

globalG.test <- SFEMethod(
    c(name = "globalG.test", title = "Global G test", uni_global),
    fun = spdep::globalG.test,
    to_df_fun = .htest2df
)

.sp.correlogram <- function(x, listw, ...) {
    sp.correlogram(neighbours = listw$neighbours, var = x, ...)
}

sp.correlogram <- SFEMethod(
    c(name = "sp.correlogram", title = "Correlogram", uni_global),
    fun = .sp.correlogram,
    to_df_fun = .correlogram2df
)

localmoran <- SFEMethod(
    c(name = "localmoran", title = "Local Moran's I", uni_local,
      default_attr = "Ii"),
    fun = spdep::localmoran,
    to_df_fun = .localmoran2df
)

localmoran_perm <- SFEMethod(
    c(name = "localmoran_perm", title = "Local Moran's I permutation testing",
      uni_local, default_attr = "Ii"),
    fun = spdep::localmoran_perm,
    to_df_fun = .localmoran2df
)
# Deal with localC and localC_perm later since it can be multivariate
# localC <- SFEMethod(
#     c(name = "localC", title = "Local Geary's C", uni_local,
#       default_attr = "localC"),
#     fun = spdep:::localC.default,
#     to_df_fun = .vec2df # Need to implement that
# )

localG <- SFEMethod(
    c(name = "localG", title = "Getis-Ord Gi(*)", uni_local,
      default_attr = "localG"),
    fun = spdep::localG,
    to_df_fun = .localG2df
)

localG_perm <- SFEMethod(
    c(name = "localG_perm", title = "Getis-Ord Gi(*) with permutation testing",
      uni_local, default_attr = "localG"),
    fun = spdep::localG_perm,
    to_df_fun = .localG2df
)

LOSH <- SFEMethod(
    c(name = "LOSH", title = "Local spatial heteroscedasticity", uni_local,
      default_attr = "Hi"),
    fun = spdep::LOSH,
    to_df_fun = function(out, nb, p.adjust.method) out
)

LOSH.mc <- SFEMethod(
    c(name = "LOSH.mc", title = "Local spatial heteroscedasticity permutation testing",
      default_attr = "Hi", uni_local),
    fun = spdep::LOSH.mc,
    to_df_fun = .LOSHmc2df
)

LOSH.cs <- SFEMethod(
    c(name = "LOSH.cs", title = "Local spatial heteroscedasticity Chi-square test",
      uni_local, default_attr = "Hi"),
    fun = spdep::LOSH.cs,
    to_df_fun = .LOSHmc2df
)

moran.plot <- SFEMethod(
    c(name = "moran.plot", title = "Moran plot", uni_local,
      default_attr = "wx"),
    fun = spdep::moran.plot,
    to_df_fun = function(out, nb, p.adjust.method) out
)
