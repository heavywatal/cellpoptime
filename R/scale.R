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
  while (nrow(x) > 1L) {
    x = nest_tippairs(x, detector = detector)
  }
  flatten_tbl_tree(x)
}

#' @details
#' `scale_branches_record` records scaling process for explanation.
#' @rdname scale
#' @export
scale_branches_record = function(x, detector = detect_driver_pois) {
  x = add_extra_columns(x)
  l = list(x)
  while (nrow(x) > 1L) {
    x = nest_tippairs(x, detector = detector)
    l[[length(l) + 1L]] = x
  }
  purrr::map(l, flatten_tbl_tree)
}

# Prepare necessary columns for scaling
add_extra_columns = function(x) {
  x %>%
    dplyr::arrange(.data$parent) %>%
    dplyr::mutate(
      isTip = !(.data$node %in% .data$parent),
      mutations = .data$branch.length,
      branch.length = pmax(.data$branch.length, 0.01),
      weight = as.integer(.data$isTip),
      term_length = 0L,
      children = list(NULL)
    ) %>%
    as_tbl_tree()
}

filter_scale_tips = function(x, detector, threshold = 0.01) {
  n = dplyr::n
  x %>%
    dplyr::filter(.data$isTip) %>%
    dplyr::mutate(total_length = .data$branch.length + .data$term_length) %>%
    dplyr::group_by(.data$parent) %>%
    dplyr::filter(n() > 1L) %>%
    dplyr::mutate(
      p_driver = detector(!!as.name("total_length")),
      weight = ifelse(!!as.name("p_driver") > threshold, !!as.name("weight"), NA_real_),
      term_length = !!as.name("weight") * !!as.name("total_length"),
      weight = sum(!!as.name("weight"), na.rm = TRUE),
      term_length = sum(!!as.name("term_length"), na.rm = TRUE) / !!as.name("weight")
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

nest_tippairs = function(x, detector) {
  nested = filter_scale_tips(x, detector = detector) %>%
    dplyr::group_by(.data$parent, .data$weight, .data$term_length) %>%
    tidyr::nest(.key = "children")
  x %>%
    dplyr::filter(is.na(.data$branch.length) | !.data$parent %in% nested$parent) %>%
    dplyr::left_join(nested, by = c(node = "parent"), suffix = c("", ".y")) %>%
    dplyr::mutate(
      isTip = .data$isTip | .data$node %in% nested$parent,
      weight = pmax(.data$weight, .data$weight.y, na.rm = TRUE),
      term_length = pmax(.data$term_length, .data$term_length.y, na.rm = TRUE),
      children = ifelse(.is_null(.data$children), .data$children.y, .data$children),
      weight.y = NULL, term_length.y = NULL, children.y = NULL
    )
}

# Rescale descendant branches recursively
rescale_descendants = function(x, scale) {
  if (!is.null(x) && scale != 1) {
    x$branch.length = x$branch.length * scale
    x$children = purrr::map(x$children, rescale_descendants, scale = scale)
  }
  x
}

flatten_tbl_tree = function(x) {
  x$weight = NULL
  x$term_length = NULL
  x %>%
    unnest_children() %>%
    dplyr::arrange(.data$node) %>%
    as_tbl_tree()
}

unnest_children = function(x) {
  has_children = !.is_null(x$children)
  if (any(has_children)) {
    unnested = tidyr::unnest(tibble::tibble(
      parent = x$node[has_children],
      data = x$children[has_children]
    ))
    x$children = NULL
    dplyr::bind_rows(unnest_children(unnested), x)
  } else {
    x$children = NULL
    x
  }
}

.is_null = function(x) {
  vapply(x, is.null, FALSE, USE.NAMES = FALSE)
}
