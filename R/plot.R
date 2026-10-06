#' View inventory in 3D
#'
#' @param inventory created by \link{aim_inventory}
#' @param x A numeric vector for translation offsets. Like in the lidR package.
#'
#' @export
aim_add_inventory3d = function(x, inventory)
{
  xyz = sf::st_coordinates(inventory)
  r  = inventory$DHP/2
  aim_circle3d(x, xyz[,1], xyz[,2], r, xyz[,3])
  rgl::texts3d(xyz[,1]-x[1], xyz[,2]-x[2]+0.1, xyz[,3]+0.1, texts = inventory$ID_Arbre, col = "white")
}

aim_circle3d <- function(x, center_x, center_y, radius, height, col = "green")
{
  n <- length(center_x)

  stopifnot(
    length(center_y) == n,
    length(radius) == n,
    length(height) == n
  )

  theta <- seq(0, 2 * pi, length.out = 50)

  # Coordinates: n circles x 50 points
  xx <- outer(center_x - x[1], rep(1, length(theta))) +
    outer(radius, cos(theta))

  yy <- outer(center_y - x[2], rep(1, length(theta))) +
    outer(radius, sin(theta))

  zz <- outer(height, rep(1, length(theta)))

  # Connect each point to the next point within each circle
  # (including the last point back to the first)
  i <- cbind(
    rep(seq_len(n), each = length(theta)),
    rep(seq_len(length(theta)), n)
  )

  j <- cbind(
    rep(seq_len(n), each = length(theta)),
    rep(c(2:length(theta), 1), n)
  )

  rgl::segments3d(
    x = cbind(xx[i], xx[j]),
    y = cbind(yy[i], yy[j]),
    z = cbind(zz[i], zz[j]),
    lwd = 5,
    col = col
  )
}

aim_view = function(las, inventory, res)
{
  X <- Y <- Z <- . <- hag <- treeID <- NULL

  if (nrow(inventory) > 1L)
    stop("Cannot process multiple trees")

  xyz = sf::st_coordinates(inventory)
  xc = xyz[,1]
  yc = xyz[,2]
  zc = xyz[,3]
  radius = inventory$DHP/2
  h2 = lidR::filter_poi(las,  Z < zc+2.2, hag > 0.25, abs(X - xc) < radius*5, abs(Y - yc) < radius*5)
  x <- arbor::plot_instance(h2, size = 2)
  aim_circle3d(x, xc, yc, radius,  zc)
  lapply(res, function(y)
  {
    nodes = y$nodes
    nodes[,1] = nodes[,1] - x[1]
    nodes[,2] = nodes[,2] - x[2]
    rgl::lines3d(nodes, col = "red", lwd = 4)
  })

  return(invisible())
}

show_fitting = function(pt, res, tol = 0.03)
{
  inliers = pt[res$inliers,]

  graphics::plot(pt, asp = 1, main = res$shape_type)
  graphics::points(res$center_x, res$center_y, pch = 3, cex = 2)
  graphics::points(inliers, col = "blue", pch = 18)
  graphics::lines(res$nodes, lwd = 2, col = "red")

  graphics::symbols(res$center_x, res$center_y,
          circles = res$radius,
          inches = FALSE, add = TRUE, fg = "purple")

  graphics::symbols(res$center_x, res$center_y,
          circles = res$radius + tol,
          inches = FALSE, add = TRUE, fg = "purple")

  if (res$radius - tol > 0)
  {
    graphics::symbols(res$center_x, res$center_y,
            circles = res$radius - tol,
            inches = FALSE, add = TRUE, fg = "purple")
  }
  # Draw 36 sectors (10-degree spacing)
  theta <- seq(0, 350, by = 10) * pi / 180

  graphics::segments(
    x0 = res$center_x,
    y0 = res$center_y,
    x1 = res$center_x + (res$radius + tol) * cos(theta),
    y1 = res$center_y + (res$radius + tol) * sin(theta),
    col = "grey80"
  )
}
