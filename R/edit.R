circloid_snaper <- function(xyz)
{
  res = arbor:::fit_circloid_cpp(xyz)
  
  plot(xyz, asp = 1, main = res$shape_type)
  points(res$center_x, res$center_y, pch = 3, cex = 2)
  lines(res$nodes, lwd = 2, col = "red")
  symbols(res$center_x, res$center_y,  circles = res$radius,  inches = FALSE, add = TRUE, fg = "purple")
  c(res$center_x, res$center_y)
}

#' Interactively edit the position of inventory trees on a point cloud
#'
#' Opens a 3D \code{rgl} window showing the point cloud (restricted to the
#' stem base) and the inventory trees (point, DBH circle and label). Each tree
#' can be selected, then moved by clicking on a cloud point. The clicked
#' position is be refined by fitting a circle to snap to the tree.
#' The function blocks until the 3D window is closed, then returns the updated
#' inventory.
#'
#' @section Interaction:
#' \itemize{
#'   \item Left drag: rotate.
#'   \item Right drag: pan.
#'   \item Wheel: zoom.
#'   \item Right click: select a tree, then snap it in \strong{XY} mode.
#'   \item Middle click: select a tree, then snap it in \strong{XYZ} mode.
#'   \item Click on the selected tree again: cancel the selection.
#'   \item Close the 3D window (or interrupt the console): save and exit.
#' }
#'
#' @section Snapping modes:
#' In both modes the tree is centred on the clicked XY position. Only the height of the slice
#' differs:
#' \itemize{
#'   \item \strong{XY mode} (right click): The tree in the inventory is snapped whatever the 
#'   Z of the clicked point. The Z of the tree is left unchanged.
#'   \item \strong{XYZ mode} (middle click): The Z of the tree is updated to the clicked Z.
#' }
#'
#' @param las A \code{LAS} object (\pkg{lidR}) with the attributes \code{hag}
#' (height above ground) and \code{treeID}.
#' @param inventory An \code{sf} object. See \link{aim_invenotry}.
#' @param pick_radius numeric. Maximum distance, in screen pixels, between the
#' mouse and a tree or a cloud point for it to be picked. Default is 15.
#' @param z_range numeric vector of length 2. Range of height above ground
#' (\code{hag}) of the points displayed in the 3D window. Default is
#' \code{c(0.5, 4)}.
#' @param slice_dz numeric. Half-thickness of the horizontal slice used for the
#' snap refinement. The slice spans \code{z - slice_dz} to \code{z + slice_dz}.
#' Default is 0.10.
#' @param radius_factor numeric. The slice is limited to a disc centred on the
#' clicked XY with radius \code{radius_factor * DHP / 2}. Default is 1.5.
#'
#' @return An \code{sf} object identical to \code{inventory} except for its
#' geometry, which is now \code{XYZ} with the edited tree positions.
#'
#' @export
aim_inventory_editor <- function(las, inventory, pick_radius = 15, z_range = c(0.5, 4), slice_dz = 0.10, radius_factor = 1.5)
{
  snap_fun <- circloid_snaper
  bottom  <- lidR::filter_poi(las, hag > z_range[1], hag < z_range[2])
  
  treeloc = bottom@data[, .(X = mean(X), Y = mean(Y)), by = treeID]
  treeloc = sf::st_as_sf(treeloc, coords = c("X", "Y"))
  near <- sf::st_is_within_distance(treeloc, inventory, dist = 4)
  treeloc <- treeloc[lengths(near) > 0, ]
  
  bottom = filter_poi(bottom, treeID %in% treeloc$treeID)
  
  offsets <- arbor::plot_instance(bottom)
  
  pc_xyz <- bottom@data
  pc_xyz <- as.matrix(pc_xyz[, 1:3])
  pc_view <- pc_xyz
  pc_norm <- sweep(pc_view, 2, c(offsets, 0))
  
  cloud_X <- las@data$X
  cloud_Y <- las@data$Y
  cloud_Z <- las@data$Z
  
  inv_xyz <- sf::st_coordinates(inventory)
  inv_xyz[,1] = inv_xyz[,1] - offsets[1]
  inv_xyz[,2] = inv_xyz[,2] - offsets[2]
  if (ncol(inv_xyz) < 3) inv_xyz <- cbind(inv_xyz, Z = min(pc_xyz[, 3]))
  
  inv_df <- data.frame(
    ID = if ("ID_Arbre" %in% names(inventory)) as.character(inventory$ID_Arbre) else as.character(seq_len(nrow(inventory))),
    X = inv_xyz[, 1], Y = inv_xyz[, 2], Z = inv_xyz[, 3],
    Radius = if ("DHP" %in% names(inventory)) inventory$DHP / 2 else rep(1, nrow(inventory)),
    stringsAsFactors = FALSE
  )
  n_trees <- nrow(inv_df)
  
  # ---------- scene ----------
  
  dev <- rgl::cur3d()
  sub <- rgl::currentSubscene3d(dev)
  
  st <- new.env()
  st$sel <- NA_integer_
  st$moved <- rep(FALSE, n_trees)
  st$ids <- vector("list", n_trees)
  st$status <- NULL
  
  set_status <- function(msg) 
  {
    if (!is.null(st$status)) try(rgl::pop3d(id = st$status), silent = TRUE)
    st$status <- unname(rgl::title3d(main = msg, col = "white", cex = 1))
  }
  
  tree_xyz <- function() 
  {
    cbind(inv_df$X, inv_df$Y, inv_df$Z)
  }
  
  draw_tree <- function(i, state = "normal")
  {
    if (length(st$ids[[i]])) try(rgl::pop3d(id = st$ids[[i]]), silent = TRUE)
    col <- switch(state, selected = "yellow", moved = "#00e676", "orange")
    p <- tree_xyz()[i, ]
    theta <- seq(0, 2 * pi, length.out = 60)
    a <- rgl::points3d(p[1], p[2], p[3], col = col, size = if (state == "selected") 12 else 8)
    b <- rgl::lines3d(p[1] + inv_df$Radius[i] * cos(theta),
                      p[2] + inv_df$Radius[i] * sin(theta),
                      rep(p[3], length(theta)),
                      col = col, lwd = if (state == "selected") 5 else 3)
    d <- rgl::texts3d(p[1], p[2], p[3] + 0.1, texts = inv_df$ID[i], col = col)
    st$ids[[i]] <- unname(c(a, b, d))
  }
  
  redraw <- function(i)
  {
    draw_tree(i, if (!is.na(st$sel) && st$sel == i) "selected" else if (st$moved[i]) "moved" else "normal")
  }
  
  for (i in seq_len(n_trees)) 
  {
    draw_tree(i)
  }
  
  idle_msg <- "CLICK (right/middle) a tree to select it  |  right-drag = pan  |  close window = save & exit"
  set_status(idle_msg)
  
  # ---------- picking (screen space) ----------
  pick <- function(mx, my, xyz, by = c("screen", "depth"))
  {
    by <- match.arg(by)
    vp <- rgl::par3d("viewport", dev = dev, subscene = sub)
    w <- vp[3]; h <- vp[4]
    wc <- rgl::rgl.user2window(xyz[, 1], xyz[, 2], xyz[, 3], projection = rgl::rgl.projection(dev = dev, subscene = sub))
    # window y starts at the bottom, mouse y at the top
    d <- sqrt((wc[, 1] * w - mx)^2 + (wc[, 2] * h - (h - my))^2)
    cand <- which(is.finite(d) & d <= pick_radius & wc[, 3] >= 0 & wc[, 3] <= 1)
    if (!length(cand)) return(NA_integer_)
    if (by == "screen") cand[which.min(d[cand])] else cand[which.min(wc[cand, 3])]
  }
  
  # ---------- RANSAC refinement ----------
  # x, y: clicked XY in the normalised (offset-subtracted) system.
  # z: height at which the slice is extracted (Z is not offset, so it is
  #    the same in the normalised and original systems).
  #    XY mode  -> Z of the tree from the inventory
  #    XYZ mode -> Z of the clicked cloud point
  # Returns the refined centre (normalised system) or NULL if no refinement.
  refine_with_snaper <- function(i, x, y, z)
  {
    # back to original coordinates to query the full cloud
    cx <- x + offsets[1]
    cy <- y + offsets[2]
    r  <- inv_df$Radius[i] * radius_factor
    
    keep <- cloud_Z >= z - slice_dz & cloud_Z <= z + slice_dz &  cloud_X >= cx - r & cloud_X <= cx + r & cloud_Y >= cy - r & cloud_Y <= cy + r
    idx <- which(keep)
    idx <- idx[(cloud_X[idx] - cx)^2 + (cloud_Y[idx] - cy)^2 <= r^2]
    
    # not enough points in the slice: no snap (adjust the minimum to what the fitter needs)
    if (length(idx) < 5)
    {
      message("Slice at Z = ", round(z, 2), " has only ", length(idx), " points; no snap.")
      return(NULL)
    }
    
    xyz <- cbind(X = cloud_X[idx], Y = cloud_Y[idx], Z = cloud_Z[idx])
    
    fit <- tryCatch(snap_fun(xyz), error = function(e) 
    {
      message("Snap function error: ", conditionMessage(e))
      NULL
    })
    
    if (is.null(fit) || length(fit) < 2 || !all(is.finite(fit[1:2]))) return(NULL)
    
    c(fit[[1]] - offsets[1], fit[[2]] - offsets[2])
  }
  
  # with_z = FALSE (right click): snap XY only
  #   -> slice at the clicked XY, at the Z of the tree from the inventory
  # with_z = TRUE  (middle click): snap XY + Z
  #   -> slice at the clicked XY, at the clicked Z
  handle_click <- function(mx, my, with_z) 
  {
    if (is.na(st$sel)) 
    {
      i <- pick(mx, my, tree_xyz(), "screen")
      if (is.na(i)) return(set_status("No tree under cursor. Zoom in or click closer to a tree centre."))
      st$sel <- i
      draw_tree(i, "selected")
      set_status(sprintf("%s selected: RIGHT-CLICK a cloud point = snap XY | MIDDLE-CLICK = snap XYZ (click the tree again to cancel)", inv_df$ID[i]))
    }
    else
    {
      i <- st$sel
      hit_tree <- pick(mx, my, tree_xyz()[i, , drop = FALSE], "screen")
      if (!is.na(hit_tree)) 
      {
        st$sel <- NA_integer_
        redraw(i)
        return(set_status(idle_msg))
      }
      j <- pick(mx, my, pc_norm, "depth")
      if (is.na(j)) return(set_status("No cloud point there. Try again (zoom in for more precision)."))
      
      new_x <- pc_norm[j, 1]
      new_y <- pc_norm[j, 2]
      new_z <- pc_norm[j, 3]
      
      # XY mode : slice at the tree's current (inventory) Z, whatever the clicked Z
      # XYZ mode: slice at the clicked Z, whatever the tree's Z
      slice_z <- if (with_z) new_z else inv_df$Z[i]
      
      refined <- refine_with_snaper(i, new_x, new_y, slice_z)
      if (!is.null(refined))
      {
        new_x <- refined[1]
        new_y <- refined[2]
      }
      
      inv_df$X[i] <<- new_x
      inv_df$Y[i] <<- new_y
      if (with_z) inv_df$Z[i] <<- new_z
      st$moved[i] <- TRUE
      st$sel <- NA_integer_
      redraw(i)
      set_status(sprintf("%s snapped (%s%s).  %s", inv_df$ID[i], if (with_z) "XYZ" else "XY", if (!is.null(refined)) " + Snaper" else "", idle_msg))
    }
  }
  
  # ---------- right button: pan on drag, select/snap(XY) on click ----------
  pan <- new.env()
  begin <- function(x, y) 
  {
    st$x0 <- x; st$y0 <- y; st$x <- x; st$y <- y
    active <- rgl::par3d("activeSubscene", dev = dev)
    pan$listeners <- rgl::par3d("listeners", dev = dev, subscene = active)
    for (s in pan$listeners) 
    {
      init <- rgl::par3d(c("userProjection", "viewport"), dev = dev, subscene = s)
      init$pos <- c(x / init$viewport[3], 1 - y / init$viewport[4], 0.5)
      pan[[as.character(s)]] <- init
    }
  }
  pan_update <- function(x, y) 
  {
    st$x <- x; st$y <- y
    for (s in pan$listeners) 
    {
      init <- pan[[as.character(s)]]
      xlat <- 2 * (c(x / init$viewport[3], 1 - y / init$viewport[4], 0.5) - init$pos)
      M <- rgl::translationMatrix(xlat[1], xlat[2], xlat[3])
      rgl::par3d(userProjection = M %*% init$userProjection, dev = dev, subscene = s)
    }
  }
  make_end <- function(with_z) 
  {
    function(...) 
    {
      if (sqrt((st$x - st$x0)^2 + (st$y - st$y0)^2) < 5) 
      {
        tryCatch(handle_click(st$x0, st$y0, with_z), error = function(e) message("Click error: ", conditionMessage(e)))
      }
    }
  }
  
  # ---------- middle button: click only (select/snap XYZ) ----------
  begin_click <- function(x, y) 
  {
    st$x0 <- x; st$y0 <- y; st$x <- x; st$y <- y
  }
  update_click <- function(x, y) 
  {
    st$x <- x; st$y <- y
  }
  
  mm <- rgl::par3d("mouseMode")
  mm["left"] <- "trackball"; mm["wheel"] <- "pull"
  rgl::par3d(mouseMode = mm)
  rgl::rgl.setMouseCallbacks(2, begin,       pan_update,   make_end(FALSE), dev = dev, subscene = sub)  # right
  rgl::rgl.setMouseCallbacks(3, begin_click, update_click, make_end(TRUE),  dev = dev, subscene = sub)  # middle
  
  cat("\n=== INTERACTIVE SNAPPING ===\n",
      "Left drag: rotate | Right drag: pan | Wheel: zoom\n",
      "Right click: select tree / snap XY\n",
      "Middle click: select tree / snap XYZ\n",
      "Close the 3D window (or press Esc in the console) to save & exit.\n", sep = "")
  
  # ---------- wait until the window is closed ----------
  tryCatch(
  {
    while (dev %in% rgl::rgl.dev.list()) Sys.sleep(0.05)
  }, 
  interrupt = function(e) 
  {
     if (dev %in% rgl::rgl.dev.list()) rgl::close3d(dev = dev)
  })
  
  cat(sprintf("Done. %d tree(s) moved.\n", sum(st$moved)))
  
  # ---------- rebuild sf ----------
  inv_out <- inv_df
  inv_out$X <- inv_out$X + offsets[1]
  inv_out$Y <- inv_out$Y + offsets[2]
  
  updated_sf <- inventory
  new_geom <- sf::st_as_sf(inv_out, coords = c("X", "Y", "Z"), dim = "XYZ", crs = sf::st_crs(inventory))
  sf::st_geometry(updated_sf) <- sf::st_geometry(new_geom)
  updated_sf
}