aim_fit_slice = function(slice, xc, yc, zc, radius, zoffset, plot = FALSE)
{
  X <- Y <- Z <- . <- hag <- treeID <- NULL

  if (methods::is(slice, "LAS"))
    slice <- slice@data

  if (zc + zoffset < 1) return(NULL)

  slice = slice[Z > zc + zoffset - 0.15 & Z < zc + zoffset + 0.15, .(X,Y,Z,hag,treeID)]

  if (nrow(slice) < 3)  return(NULL)
  if (max(slice$hag) < 0.5) return(NULL)

  xyz <- as.matrix(slice[, .(X,Y,Z)])
  res <- arbor:::fit_circloid_cpp(xyz, tolerance = 0.03, complexity = 3)

  if (plot)
  {
    show_fitting(xyz, res)
    graphics::points(xc, yc, pch = 19, cex = 4, col = "red")
  }

  if (res$covered_arc_degree > 180)
  {
    inliner = slice[res$inliers]
    seeds = slice[res$inliers, .(X,Y,Z)]
    res$to_be_merge_id = unique(inliner$treeID)
    res$seeds = seeds
    class(res) <- c("aim_fit", class(res))
  }
  else
  {
    res = NULL
  }

  res
}

aim_fit_tree = function(las, inventory)
{
  X <- Y <- Z <- . <- hag <- treeID <- NULL

  if (nrow(inventory) > 1L) stop("Cannot process multiple trees")

  xyz <- sf::st_coordinates(inventory)
  xc  <- xyz[,1]
  yc  <- xyz[,2]
  zc  <- xyz[,3]
  radius <- inventory$DHP/2
  search_radius_factor <- 1.5
  offsets <- seq(-2,2, by = 0.5)
  res <- vector("list", length(offsets))

  pc <- las@data[abs(X - xc) < radius*5 & abs(Y - yc) < radius*5, .(X,Y,Z,hag,treeID)]

  it <- 0
  while (it <= 3)
  {
    it = it + 1
    slice = pc[abs(X - xc) < radius*search_radius_factor & abs(Y - yc) < radius*search_radius_factor]

    for (i in seq_along(offsets))
    {
      res[[i]] <- aim_fit_slice(slice, xc, yc, zc, radius, offsets[i])
    }

    res <- Filter(Negate(is.null), res)

    if (length(res) == 0)
    {
      search_radius_factor = search_radius_factor + 1.5
    }
    else
    {
      it = .Machine$integer.max
    }
  }

  # No fit, no tree detected
  if (length(res) == 0)
  {
    slice <- slice[Z > zc - 0.5 & Z < zc + 0.5 & (X - xc)^2 < 2*radius &(Y - yc)^2 < 2*radius]
    candidates <- unique(slice$treeID)
    if (length(candidates) == 1)
    {
      msg = paste0("No tree detected for tree " , inventory$ID_Arbre, " but a single candidate ID was found and used.")
      res <- list(list(to_be_merge_id = candidates))
      class(res[[1]]) <- c("aim_fit", class(res))
      warning(msg, call. = FALSE)
    }
    else
    {
      msg <- paste0("No tree detected for tree " , inventory$ID_Arbre)
      res <- list(list(to_be_merge_id = NA_integer_))
      class(res[[1]]) <- c("aim_fit", class(res))
      warning(msg, call. = FALSE)
    }
  }

  res
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
