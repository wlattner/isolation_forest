#' Detect anomalies using an ensemble of decision trees.
#'
#' @param df A \code{data.frame} with observations.
#' @param n_trees The number of trees to fit.
#' @param subsample_size The number of examples to use in fitting each tree.
#' @return ret An object with class \code{isolation_forest}.
#' @export
isolation_forest <- function(df, n_trees = 10, subsample_size = 256) {
  df[] <- sapply(df, as.numeric)
  height_limit <- ceiling(log2(subsample_size))
  trees <- vector("list", n_trees)
  for (i in seq_along(trees)) {
    ind <- sample(1:nrow(df), subsample_size, replace = FALSE)
    trees[[i]] <- .fit_tree(df[ind, ], height_limit)
  }

  # the "expected" depth of an arbitrary example, this will be used to scale
  # the path depth scores at prediciton time
  c_n <- .adjust_path_length(subsample_size)

  ret <- list(trees = trees, c_n = c_n)
  structure(ret, class = "isolation_forest")
}

.fit_tree <- function(df, height_limit) {
  n <- list(right = NULL, left = NULL, is_leaf = FALSE, n_examples = nrow(df))

  unique_vals <- lapply(df, unique)
  n_unique <- sapply(unique_vals, length)

  if (nrow(df) <= 1 || height_limit < 1 || !any(n_unique > 1)) {
    n$is_leaf <- TRUE
    return(n)
  }

  n$split_var <- sample(which(n_unique > 1), size = 1)
  candidate_splits <- unique_vals[[n$split_var]]
  # remove the min value from the candidate splits as this will createn a
  # partition with zero examples in the right child
  candidate_splits <- candidate_splits[-which.min(candidate_splits)]
  if (length(candidate_splits) == 1) {
    # avoid bad (unexpected) behavior with sample and sample.int
    n$split_val <- candidate_splits[[1]]
  } else {
    n$split_val <- sample(candidate_splits, size = 1)
  }
  is_gte <- df[, n$split_var] >= n$split_val
  right_df <- df[is_gte, ]
  left_df <- df[!is_gte, ]

  n$right <- .fit_tree(right_df, height_limit - 1)
  n$left <- .fit_tree(left_df, height_limit - 1)
  return(n)
}

# do we know if this is right? it should be the expected path length for an
# unsuccessful search in a binary search tree with n items
.adjust_path_length <- function(n) {
  2 * (log(n) + 0.5772156649) - 2 * (n - 1) / n
}

#' @export
predict.isolation_forest <- function(object, new_data) {
  path_lengths <- sapply(object$trees, .path_length, df = new_data)
  # we now have a matrix of trees x rows
  avg_lengths <- rowMeans(path_lengths)
  2 ** (-avg_lengths / object$c_n)
}

.path_length <- function(tree, df) {
  l <- vector("integer", nrow(df))
  for (i in 1:nrow(df)) {
    node <- tree
    l[i] <- 0
    while (node$is_leaf) {
      if (df[i, node$split_var] < node$split_val) {
        node <- node$left
      } else {
        node <- node$right
      }
      l[i] <- l[i] + 1
    }
    # adjust the height for the unbuilt subtree
    l[i] <- l[i] + .adjust_path_length(node$n_examples)
  }
  l
}
