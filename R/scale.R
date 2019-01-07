#' Scale branches
#'
#' @details
#' `scale_branches` is a shortcut.
#' @param x data.frame
#' @param detector function
#' @rdname scale
#' @export
scale_branches = function(x, detector = detect_driver_pois) {
  x = add_extra_columns(x)
  num_edges = nrow(x)
  while (nrow(x) > 1L) {
    x = nest_tippairs(x, detector = detector)
  }
  while (nrow(x) < num_edges) {
    x = unnest_children(x)
  }
  x %>%
    dplyr::select(-.data$children) %>%
    dplyr::arrange(.data$node) %>%
    as_tbl_tree()
}

#' @details
#' `scale_branches_record` records scaling process for explanation.
#' @rdname scale
#' @export
scale_branches_record = function(x) {
  x = add_extra_columns(x)
  num_edges = nrow(x)
  l = list(x)
  while (nrow(x) > 1L) {
    x = nest_tippairs(x)
    l = c(l, list(x))
  }
  purrr::map(l, ~ {
    while (nrow(.x) < num_edges) {
      .x = unnest_children(.x)
    }
    .x %>%
      dplyr::select(-.data$children) %>%
      dplyr::arrange(.data$node) %>%
      as_tbl_tree()
  })
}

# Prepare necessary columns for scaling
add_extra_columns = function(x) {
  x %>%
    dplyr::arrange(.data$parent) %>%
    dplyr::mutate(
      is_tip = !(.data$node %in% .data$parent),
      mutations = .data$branch.length,
      branch.length = pmax(.data$branch.length, 0.01),
      term_length = 0,
      children = list(NULL)
    ) %>%
    as_tbl_tree()
}

#' @rdname scale
#' @export
filter_scale_tips = function(x, detector) {
  n = dplyr::n
  x %>%
    dplyr::filter(.data$is_tip) %>%
    dplyr::mutate(total_length = .data$branch.length + .data$term_length) %>%
    dplyr::group_by(.data$parent) %>%
    dplyr::filter(n() > 1L) %>%
    dplyr::mutate(
      p_driver = detector(!!as.name("total_length")),
      term_length = min(!!as.name("total_length")),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      scale = .data$term_length / .data$total_length,
      branch.length = .data$branch.length * .data$scale,
      children = purrr::map2(.data$children, .data$scale, rescale_descendants),
      scale = NULL,
      total_length = NULL
    )
}

#' @rdname scale
#' @export
nest_tippairs = function(x, detector) {
  nested = filter_scale_tips(x, detector = detector) %>%
    dplyr::group_by(.data$parent, .data$term_length) %>%
    tidyr::nest(.key = "children")
  x %>%
    dplyr::filter(is.na(.data$branch.length) | !.data$parent %in% nested$parent) %>%
    dplyr::left_join(nested, by = c(node = "parent"), suffix = c("", ".y")) %>%
    dplyr::mutate(
      is_tip = .data$is_tip | .data$node %in% nested$parent,
      term_length = pmax(.data$term_length, .data$term_length.y, na.rm = TRUE),
      children = ifelse(purrr::map_lgl(.data$children, is.null), .data$children.y, .data$children),
      term_length.y = NULL, children.y = NULL
    )
}

#' @rdname scale
#' @export
unnest_children = function(x) {
  idx = !purrr::map_lgl(x$children, is.null)
  .outer = tibble::tibble(
    parent = x$node[idx],
    data = x$children[idx]
  )
  .inner = dplyr::mutate(x, children = list(NULL))
  dplyr::bind_rows(tidyr::unnest(.outer), .inner) %>%
    dplyr::mutate(is_tip = !is.na(.data$label))
}

# Rescale descendant branches recursively
rescale_descendants = function(x, scale) {
  if (is.null(x) || scale == 1) {
    x
  } else {
    dplyr::mutate(
      x,
      branch.length = .data$branch.length * scale,
      children = purrr::map(.data$children, rescale_descendants, scale)
    )
  }
}
