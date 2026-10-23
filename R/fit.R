aim_fit_tree = function(las, inventory)
{
  if (nrow(inventory) > 1L) stop("Cannot process multiple trees")

  res <- arbor:::aim_fit_tree_cpp(las@data, aim_inventory_df(inventory))
  lapply(res, function(x) { class(x) <- c("aim_fit", class(x)) ; x })
}

#' @method print aim_fit
#' @export
print.aim_fit <- function(x, ...) {
  # Helper to safely retrieve an element or return default
  get_val <- function(key, default = NULL) {
    if (!is.null(x[[key]])) x[[key]] else default
  }

  # Status
  success <- get_val("success")
  if (!is.null(success)) {
    cat(sprintf("Status              : %s\n", if (success) "SUCCESS" else "FAILED"))
  }

  # Geometry Parameters
  shape <- get_val("shape_type")
  if (!is.null(shape)) cat(sprintf("Shape Type          : %s\n", shape))

  radius <- get_val("radius")
  if (!is.null(radius)) cat(sprintf("Radius              : %.4f\n", radius))

  # Center coordinates
  cx <- get_val("center_x")
  cy <- get_val("center_y")
  cz <- get_val("center_z")
  if (!all(sapply(list(cx, cy, cz), is.null))) {
    cat(sprintf("Center (X, Y, Z)    : [%s, %s, %s]\n",
                if (!is.null(cx)) sprintf("%.3f", cx) else "NA",
                if (!is.null(cy)) sprintf("%.3f", cy) else "NA",
                if (!is.null(cz)) sprintf("%.3f", cz) else "NA"))
  }

  # Fit Quality & Coverage
  arc <- get_val("covered_arc_degree")
  if (!is.null(arc)) cat(sprintf("Covered Arc         : %g\u00b0\n", arc))

  inlier_pct <- get_val("percentage_inlier")
  if (!is.null(inlier_pct)) cat(sprintf("Inliers             : %.2f%%\n", inlier_pct))

  merge_id <- get_val("to_be_merge_id")
  if (!is.null(merge_id)) cat(sprintf("Merge ID            : %s\n", paste(merge_id, collapse = ", ")))

  nodes <- get_val("nodes")
  if (!is.null(nodes)) {
    if (is.matrix(nodes) || is.data.frame(nodes)) {
      cat(sprintf("Nodes               : matrix [%d x %d]\n", nrow(nodes), ncol(nodes)))
    } else {
      cat(sprintf("Nodes               : length %d\n", length(nodes)))
    }
  } else {
    cat("Nodes               : <NULL>\n")
  }

  inliers <- get_val("inliers")
  if (!is.null(inliers)) {
    cat(sprintf("Inlier Indices      : vector [length %d]\n", length(inliers)))
  } else {
    cat("Inlier Indices      : <NULL>\n")
  }

  # Return object invisibly per standard R print conventions
  invisible(x)
}
