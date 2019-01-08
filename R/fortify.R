#' Convert tbl_tree to plottable data.frame
#'
#' @details
#' `fortify_cellpop` prepares plottable data.frame.
#' @param model output of `scale_branches()`
#' @param data tbl_tree
#' @inheritParams filter_origins
#' @rdname fortify
#' @export
fortify_cellpop = function(model, data, method = "fdr", q = 0.05) {
  if (missing(data)) data = model
  mutant = filter_origins(model, method = method, q = q)$node
  meta_info = group_clade(model, mutant) %>%
    dplyr::select(.data$node, .data$mutations, .data$p_driver, .data$group)
  ggtree_fortify(data) %>%
    dplyr::left_join(meta_info, by = "node")
}

#' @details
#' `ggtree_fortify` prepares plottable data.frame.
#' @rdname fortify
#' @export
ggtree_fortify = function(data) {
  data = if (is.data.frame(data)) {
    as_phylo(data)
  } else if (is.list(data)) {
    as_multiphylo(data)
  } else {
    stop("Unknown class(data): ", class(data))
  }
  ggtree::fortify(data)
}

add_xy_coord = function(data, ladderize = TRUE) {
  data = add_x_coord(data)
  if (ladderize) {
    data = .ladderize(data)
  }
  add_y_coord(data)
}

# More efficient than ggtree::getYcoord for >200 tips
get_x_coord = function(data) {
  len = dplyr::coalesce(data$branch.length, 0)
  x = numeric(length(len))
  nextnode = data$node
  root = nextnode[nextnode == data$parent]
  while (any(nextnode != root)) {
    x = x + len[nextnode]
    nextnode = data$parent[nextnode]
  }
  x
}

add_x_coord = function(data) {
  data$x = get_x_coord(data)
  data
}

# More efficient than ggtree::getYcoord for >2000 tips
add_y_coord = function(data, step = 1) {
  n = dplyr::n  # for speed and warning
  num_tips = NROW(data) - sum(data$node %in% data$parent)
  current = seq_len(num_tips)
  data$y = NA_real_
  data$y[data$node <= num_tips] = current * step
  data = dplyr::arrange(data, .data$node)
  while (anyNA(data$y)) {
    df_pairs = data %>%
      dplyr::filter(.data$node %in% current) %>%
      dplyr::group_by(.data$parent) %>%
      dplyr::filter(n() == 2L)
    df_new = dplyr::summarise(df_pairs, y = mean(!!as.name("y"))) %>% dplyr::ungroup()
    data$y[df_new$parent] = df_new$y
    current = c(setdiff(current, df_pairs$node), df_new$parent)
  }
  data
}

.ladderize = function(data, right = FALSE) {
  phylo = ape::ladderize(as_phylo(data), right = right)
  tibble::as_tibble(phylo$edge) %>%
    dplyr::full_join(data, by = c("parent", "node"))
}

# shortcut for fortify.phylo
get_y_coord = function(data) {
  add_y_coord(data)$y
}
