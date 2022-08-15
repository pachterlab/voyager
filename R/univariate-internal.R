# Internal function for univariate metrics
.calc_univar_autocorr <- function(x, listw, fun, BPPARAM, returnDF = FALSE, ...) {
  if (is.null(listw))
    stop("The graph specified is absent from the SFE object.")
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  out_list <- bplapply(seq_len(nrow(x)), function(i) {
    fun(x[i,], listw, ...)
  }, BPPARAM = BPPARAM)
  if (returnDF) {
    out_list <- lapply(out_list, unlist, use.names = TRUE)
    out <- Reduce(rbind, out_list)
    if (!is.matrix(out)) out <- t(as.matrix(out))
    rownames(out) <- rownames(x)
    out <- DataFrame(out)
  } else {
    names(out_list) <- rownames(x)
    out <- out_list
  }
  return(out)
}

#' @importFrom spdep include.self nb2listw
.calc_univar_sfe_fun <- function(fun, returnDF = FALSE, include_self = FALSE) {
  function(x, features = NULL, colGraphName = 1L, sample_id = NULL,
           exprs_values = "logcounts", BPPARAM = SerialParam(),
           zero.policy = NULL, ...) {
    # Am I sure that I want to use logcounts as the default?
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    out <- lapply(sample_id, function(s) {
      features <- .check_features(x, features)[["assay"]]
      listw_use <- colGraph(x, type = colGraphName, sample_id = s)
      if (include_self) {
        nb2 <- include.self(listw_use$neighbours)
        listw_use <- nb2listw(nb2)
      }
      mat <- assay(x, exprs_values)[features, colData(x)$sample_id == s]
      o <- fun(mat, listw_use, BPPARAM = BPPARAM, zero.policy = zero.policy, ...)
      if (returnDF) o$sample_id <- s
      o
    })
    if (returnDF && length(sample_id) > 1L) {
      out <- do.call(rbind, out)
    } else {
      names(out) <- sample_id
    }
    if (length(sample_id) == 1L) out <- out[[1]]
    out
  }
}

.df_univar_autocorr <- function(df, listw, features, fun, BPPARAM,
                                zero.policy, ...) {
  if (is(df, "sf")) df <- st_drop_geometry(df)
  mat <- t(as.matrix(df[, features, drop = FALSE]))
  if (anyNA(mat)) {
    stop("Only numeric columns without NA (within the sample_id) can be used.")
  }
  fun(mat, listw, BPPARAM = BPPARAM, zero.policy = zero.policy, ...)
}

.MoransI2df <- function(out, name) {
  # out should already be a DFrame when this function is called
  names(out)[1] <- name
  out
}

.add_local_res <- function(df, res) {
  names_existing <- intersect(names(df), names(res))
  if (length(names_existing)) df[,names_existing] <- NULL
  cbind(df, res)
}

.local_feature_dfs <- function(res, name) {
  # Output for each sample: list whose names are the features
  res <- lapply(names(res), function(n) {
    r <- as(res[[n]], "DataFrame")
    names(r) <- paste(name, names(res), n, sep = "_")
    r
  })
  do.call(cbind, res)
}

.coldata_univar_fun <- function(fun, to_df_fun = NULL, name, to_df_params = list(),
                                include_self = FALSE, local = FALSE) {
  function(x, features, colGraphName = 1L, sample_id = NULL, BPPARAM = SerialParam(),
           zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    for (s in sample_id) {
      listw_use <- colGraph(x, type = colGraphName, sample_id = s)
      if (include_self) {
        nb2 <- include.self(listw_use$neighbours)
        listw_use <- nb2listw(nb2)
      }
      res <- .df_univar_autocorr(colData(x)[colData(x)$sample_id == s,],
                                 listw_use, features, fun,
                                 BPPARAM, zero.policy, ...)
      if (local) {
        res <- .local_feature_dfs(res, name)
        colData(x)[colData(x)$sample_id == s, names(res)] <- res
      } else {
        x <- .add_fd_dimData(x, MARGIN = 2, res, features, s, to_df_fun, name,
                             to_df_params)
      }
    }
    x
  }
}

.colgeom_univar_fun <- function(fun, to_df_fun = NULL, name, to_df_params = list(),
                                include_self = FALSE, local = FALSE) {
  function(x, features, colGeometryName = 1L, colGraphName = 1L, sample_id = NULL,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    for (s in sample_id) {
      listw_use <- colGraph(x, type = colGraphName, sample_id = s)
      if (include_self) {
        nb2 <- include.self(listw_use$neighbours)
        listw_use <- nb2listw(nb2)
      }
      cg <- colGeometry(x, type = colGeometryName, sample_id = s)
      res <- .df_univar_autocorr(cg, listw_use, features, fun, BPPARAM,
                                 zero.policy, ...)
      if (local) {
        res <- as.data.frame(res)
        names(res) <- paste(name, names(res), sep = "_")
        cg <- cbind(cg, res)
        colGeometry(x, colGeometryName, sample_id = s) <- cg
      } else {
        colGeometry(x, colGeometryName, sample_id = "all") <-
          .add_fd(x, colGeometry(x, colGeometryName, sample_id = "all"),
                  res, features, s, to_df_fun, name, to_df_params)
      }
    }
    x
  }
}

.annotgeom_univar_fun <- function(fun, to_df_fun = NULL, name, to_df_params = list(),
                                  include_self = FALSE, local = FALSE) {
  function(x, features, annotGeometryName = 1L, annotGraphName = 1L, sample_id = NULL,
           BPPARAM = SerialParam(), zero.policy = NULL, ...) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    for (s in sample_id) {
      listw_use <- annotGraph(x, type = annotGraphName, sample_id = s)
      if (include_self) {
        nb2 <- include.self(listw_use$neighbours)
        listw_use <- nb2listw(nb2)
      }
      ag <- annotGeometry(x, type = annotGeometryName, sample_id = s)
      ag <- .rm_empty_geometries(ag, MARGIN = 3)
      res <- .df_univar_autocorr(ag, listw_use, features, fun, BPPARAM,
                                 zero.policy, ...)
      if (local) {
        res <- as.data.frame(res)
        names(res) <- paste(name, names(res), sep = "_")
        # buggy with empty geometries
        ag <- cbind(ag, res)
        annotGeometry(x, annotGeometryName, sample_id = "all") <- ag
      } else {
        annotGeometry(x, annotGeometryName, sample_id = "all") <-
          .add_fd(x, annotGeometry(x, annotGeometryName, sample_id = "all"),
                  res, features, s, to_df_fun, name, to_df_params)
      }
    }
    x
  }
}

.sfe_univar_autocorr <- function(x, features, colGraphName, sample_id,
                                 exprs_values, fun, BPPARAM, zero.policy,
                                 name) {
  sample_id <- .check_sample_id(x, sample_id, one = FALSE)
  for (s in sample_id) {
    out <- fun(x, features, colGraphName, s, exprs_values, BPPARAM,
               zero.policy)
    out$sample_id <- NULL # Not necessary here, but necessary when just calling calculateMoransI
    out <- .MoransI2df(out, name)
    out <- .add_name_sample_id(out, s)
    features <- .symbol2id(x, features)
    rowData(x)[features, names(out)] <- out
  }
  x
}

# For the `nsim` and `alternative` arguments
.calc_univar_sfe_fun_mc <- function(fun) {
  function(x, features, colGraphName = 1L, sample_id = NULL, nsim,
           exprs_values = "logcounts",
           BPPARAM = SerialParam(), zero.policy = NULL,
           alternative = "greater", ...) {
    .calc_univar_sfe_fun(fun)(x, features, colGraphName, sample_id,
                              exprs_values, BPPARAM, zero.policy, nsim = nsim,
                              alternative = alternative, ...)
  }
}

.coldata_univar_fun_mc <- function(fun, to_df_fun, name) {
  function(x, features, colGraphName = 1L, sample_id = NULL, nsim,
           BPPARAM = SerialParam(), zero.policy = NULL, alternative = "greater",
           ...) {
    .coldata_univar_fun(fun, to_df_fun, name)(x, features, colGraphName,
                                              sample_id, BPPARAM,
                                              zero.policy, nsim = nsim,
                                              alternative = alternative, ...)
  }
}

.MoranMC2df <- function(out, name) {
  # Convert results to DFrame. I'll write my own plotting function based on ggplot2.
  out <- lapply(out, function(o) {
    o$res <- I(list(o$res))
    DataFrame(unclass(o))
  })
  rns <- names(out)
  out <- Reduce(rbind, out)
  rownames(out) <- rns
  colnames(out) <- paste(name, colnames(out), sep = "_")
  out
}

.colgeom_univar_fun_mc <- function(fun, to_df_fun, name) {
  function(x, features, colGeometryName = 1L, colGraphName = 1L,
           sample_id = NULL, nsim, BPPARAM = SerialParam(), zero.policy = NULL,
           alternative = "greater",
           ...) {
    .colgeom_univar_fun(fun, to_df_fun, name)(x, features, colGeometryName,
                                              colGraphName, sample_id, BPPARAM,
                                              zero.policy, nsim = nsim,
                                              alternative = alternative, ...)
  }
}

.annotgeom_univar_fun_mc <- function(fun, to_df_fun, name) {
  function(x, features, annotGeometryName = 1L, annotGraphName = 1L,
           sample_id = NULL, nsim, BPPARAM = SerialParam(), zero.policy = NULL,
           alternative = "greater",
           ...) {
    .annotgeom_univar_fun(fun, to_df_fun, name)(
      x, features, annotGeometryName, annotGraphName, sample_id, BPPARAM,
      zero.policy, nsim = nsim, alternative = alternative, ...)
  }
}

.sfe_univar_mc <- function(x, features, colGraphName, sample_id, nsim,
                           exprs_values, fun, BPPARAM, zero.policy, alternative,
                           name, ...) {
  sample_id <- .check_sample_id(x, sample_id, one = FALSE)
  for (s in sample_id) {
    out <- fun(x, features, colGraphName, s, nsim, exprs_values, BPPARAM,
               zero.policy, alternative, ...)
    out <- .MoranMC2df(out, name)
    features <- .symbol2id(x, features)
    out <- .add_name_sample_id(out, s)
    rowData(x)[features, names(out)] <- out
  }
  x
}

.correlogram2df <- function(out, name, method) {
  if (method %in% c("I", "C")) {
    out <- lapply(out, function(o) {
      colnames(o$res) <- c(method, "expectation", "variance")
      o
    })
  }
  out <- lapply(out, function(o) o$res)
  out_df <- DataFrame(res = I(out))
  names(out_df) <- name
  out_df
}

.sfe_univar_local_fun <- function(fun) {
  function(x, features, colGraphName, sample_id, exprs_values, fun, BPPARAM,
           zero.policy, name, ...) {
    sample_id <- .check_sample_id(x, sample_id, one = FALSE)
    for (s in sample_id) {
      out <- fun(x, features, colGraphName, s, exprs_values, BPPARAM,
                 zero.policy)

    }
    x
  }
}
