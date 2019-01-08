#' Convert tbl_tree to plottable data.frame
#'
#' @details
#' `fortify_cellpop` prepares plottable data.frame.
#' @param ... passed to `layout_rect` or `ggtree::fortify`.
#' @param model output of `scale_branches()`
#' @param data tbl_tree
#' @inheritParams filter_origins
#' @rdname fortify
#' @export
fortify_cellpop = function(model, data, ..., method = "fdr", q = 0.05) {
  if (missing(data)) data = model
  mutant = filter_origins(model, method = method, q = q)$node
  meta_info = group_clade(model, mutant) %>%
    dplyr::select(.data$node, .data$mutations, .data$p_driver, .data$group)
  layout_rect(data) %>%
    dplyr::left_join(meta_info, by = "node", suffix = c(".x", "")) %>%
    dplyr::select(-dplyr::ends_with(".x"))
}

#' @details
#' `ggtree_fortify` is a wrapper of `ggtree::fortify`.
#' @rdname fortify
#' @export
ggtree_fortify = function(data, ...) {
  data = if (is.data.frame(data)) {
    as_phylo(data)
  } else if (is.list(data)) {
    as_multiphylo(data)
  } else {
    stop("Unknown class(data): ", class(data))
  }
  ggtree::fortify(data, ...)
}

#' @details
#' `layout_rect` is a simplified version of `ggtree::fortify`,
#' whose output lacks "branch" and "angle" columns.
#' @param ladderize Use `ape::ladderize` or not.
#' @rdname fortify
#' @export
layout_rect = function(data, ladderize = TRUE, ...) UseMethod("layout_rect")

#' @rdname fortify
#' @export
layout_rect.list = function(data, ladderize = TRUE, ...) {
  if (is.null(names(data))) names(data) = paste0("Tree #", seq_along(data))
  purrr::map_dfr(data, layout_rect, ladderize = ladderize, .id = ".id")
}

#' @rdname fortify
#' @export
layout_rect.data.frame = function(data, ladderize = TRUE, ...) {
  data$isTip = !(data$node %in% data$parent)
  data$x = get_x_coord(data)
  data$y = get_y_coord(data, ladderize = ladderize)
  data
}

# More efficient than ggtree::getYcoord for >200 tips
get_x_coord = function(data) {
  # assuming data$node == seq_len(NROW(data))
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

# More efficient than ggtree::getYcoord for >2000 tips
get_y_coord = function(data, ladderize = TRUE, step = 1) {
  n = dplyr::n  # for speed and warning
  if (ladderize) data = .ladderize(data)
  data = data[, c("parent", "node")]
  num_tips = NROW(data) - sum(data$node %in% data$parent)
  current = seq_len(num_tips)
  data$y = NA_real_
  data$y[data$node <= num_tips] = current * step
  data = dplyr::arrange(data, .data$node)
  # assuming data$node == seq_len(NROW(data))
  while (anyNA(data$y)) {
    df_pairs = data %>%
      dplyr::filter(.data$node %in% current) %>%
      dplyr::group_by(.data$parent) %>%
      dplyr::filter(n() == 2L)
    df_new = dplyr::summarise(df_pairs, y = mean(!!as.name("y"))) %>% dplyr::ungroup()
    data$y[df_new$parent] = df_new$y
    current = c(setdiff(current, df_pairs$node), df_new$parent)
  }
  data$y
}

.ladderize = function(data, right = FALSE) {
  phylo = ape::ladderize(as_phylo(data), right = right)
  tibble::as_tibble(phylo$edge) %>%
    dplyr::full_join(data, by = c("parent", "node"))
}
