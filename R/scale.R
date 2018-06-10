#' Scale branches
#'
#' @description
#' `add_extra_columns` prepares necessary columns for scaling
#' @param x data.frame
#' @rdname scale
#' @export
add_extra_columns = function(x) {
  x %>%
    dplyr::arrange(.data$parent) %>%
    dplyr::mutate(
      is_tip = !is.na(.data$label),
      mutations = .data$branch.length,
      branch.length = pmax(.data$branch.length, 0.01),
      term_length = 0,
      exp_desc = ifelse(.data$is_tip, 1L, NA_integer_), # expected number of descendant cells
      children = list(NULL)
    )
}

#' @rdname scale
#' @export
filter_scale_tips = function(x) {
  x %>%
    dplyr::filter(.data$is_tip) %>%
    dplyr::group_by(.data$parent) %>%
    dplyr::filter(dplyr::n() > 1L) %>%
    dplyr::mutate(total_length = min(.data$branch.length + .data$term_length)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      scale = .data$total_length / (.data$branch.length + .data$term_length),
      branch.length = .data$branch.length * .data$scale,
      term_length = .data$term_length * .data$scale,
      children = purrr::map2(.data$children, .data$scale, rescale_children)
    ) %>%
    dplyr::mutate(term_length = .data$total_length) %>%
    dplyr::select(-.data$scale, -.data$total_length)
}

#' @rdname scale
#' @export
nest_tippairs = function(x) {
  nested = filter_scale_tips(x) %>%
    dplyr::group_by(.data$parent, .data$term_length) %>%
    dplyr::mutate(p_driver = detect_driver(.data$exp_desc)) %>%
    tidyr::nest(.key="children") %>%
    dplyr::mutate(exp_desc = purrr::map_int(.data$children, ~{
      #TODO weight by mutatoins
      sum(.x$exp_desc)
    }))
  x %>%
    dplyr::filter(is.na(.data$branch.length) | !.data$parent %in% nested$parent) %>%
    dplyr::left_join(nested, by=c(node="parent"), suffix=c("", ".y")) %>%
    dplyr::mutate(
      is_tip = .data$is_tip | .data$node %in% nested$parent,
      term_length = pmax(.data$term_length, .data$term_length.y, na.rm=TRUE),
      exp_desc = dplyr::coalesce(.data$exp_desc, .data$exp_desc.y),
      children = ifelse(purrr::map_lgl(.data$children, is.null), .data$children.y, .data$children),
      term_length.y = NULL, exp_desc.y = NULL, children.y = NULL)
}

#' @rdname scale
#' @export
unnest_children = function(x) {
  .outer = x %>%
    dplyr::transmute(parent=.data$node, .data$children) %>%
    dplyr::filter(!purrr::map_lgl(.data$children, is.null)) %>%
    tidyr::unnest()
  .inner = x %>% dplyr::mutate(children=list(NULL))
  dplyr::bind_rows(.outer, .inner) %>%
    dplyr::mutate(is_tip = !is.na(.data$label))
}
