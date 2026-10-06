#' Disentangle Tree IDs Between Point Cloud Data and Inventory
#'
#' Identifies overlapping tree IDs between a point cloud object and a forest
#' inventory dataset, remapping conflicting tree IDs.
#'
#' If a tree has an ID (e.g., 54) in the point cloud and the inventory also
#' contains an ID = 54, they typically correspond to two different trees (the
#' probability that they represent the same tree is extremely low). In some
#' circumstances, this creates entangled IDs. To avoid complex downstream errors,
#' users can first apply \code{aim_disentangle_ids()}. The function remaps the
#' IDs in the point cloud such that there is no ID overlap between the point
#' cloud and the inventory.
#'
#' @param las A \code{LAS} object containing point cloud data with a \code{treeID} attribute.
#' @param inventory An \code{sf} object containing forest inventory data with an \code{ID_Arbre} column.
#'
#' @seealso [aim_inventory()]
#'
#' @return A \code{LAS} object with updated tree IDs, or the original \code{LAS} object
#'   if no matching IDs are found.
#' @export
#' @md
aim_disentangle_ids <- function(las, inventory) {
  # --- Defensive Programming / Input Validation ---
  if (!inherits(las, "LAS")) {
    stop("Argument 'las' must be a valid 'LAS' object.")
  }
  
  if (!"treeID" %in% names(las)) {
    stop("Argument 'las' does not contain a 'treeID' attribute.")
  }
  
  if (!inherits(inventory, "sf") && !is.data.frame(inventory)) {
    stop("Argument 'inventory' must be an 'sf' spatial object or data frame.")
  }
  
  if (!"ID_Arbre" %in% names(inventory)) {
    stop("Argument 'inventory' does not contain the required 'ID_Arbre' column.")
  }
  
  # --- Extraction and Remapping Logic ---
  las_ids <- unique(las$treeID)
  las_ids <- las_ids[!is.na(las_ids)]  # Exclude NA values if present
  target_ids <- inventory$ID_Arbre
  
  # Identify IDs present in both LAS and inventory
  ids_to_remap <- intersect(las_ids, target_ids)
  
  if (length(ids_to_remap) == 0) {
    message("No matching treeIDs found between LAS data and inventory.")
    return(las)
  }
  
  # Determine new IDs starting from max(las$treeID) + 1
  max_id <- max(las_ids, na.rm = TRUE)
  new_ids <- seq(from = max_id + 1, length.out = length(ids_to_remap))
  
  # Create a lookup mapping from old matching IDs to new IDs
  map_vector <- stats::setNames(new_ids, ids_to_remap)
  
  # Identify points in LAS that need updating
  mask <- las$treeID %in% ids_to_remap
  
  # Perform remapping on LAS points
  las$treeID[mask] <- as.integer(map_vector[as.character(las$treeID[mask])])
  
  return(las)
}

#' Create a Tree Matching Table
#'
#' Matches LiDAR-detected trees with trees from a field inventory.
#'
#' @param las A `LAS` object instance-segmented with the arbor package.
#' @param inventory A sf data frame containing the field inventory, including
#'   `ID_Arbre` and `POM` columns. Use \link{aim_inventory}
#'
#' @return A data frame containing the matching results for each inventory tree.
#'
#' @export
aim_matching_table = function(las, inventory)
{
  required <- c("X", "Y", "Z", "treeID", "hag")
  missing <- setdiff(required, names(las))
  if (length(missing) > 0)
    stop("Missing required column(s): ", paste(missing, collapse = ", "))
  
  treeID = unique(las$treeID)
  targetID = inventory$ID_Arbre
  
  if (any(treeID %in% targetID)) {
    stop("Some tree IDs in the point cloud are also in the inventory tree IDs. This can cause matching issues. Call aim_disentangle_ids() first.")
  }
  
  hag <- NULL
  z   <- max(inventory$POM, na.rm = TRUE) + 2.5
  useful <- lidR::filter_poi(las, hag <= z)
  
  n <- nrow(inventory)
  results <- vector("list", n)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  
  for (i in seq_len(n))
  {
    tree = inventory[i, ]
    res = aim_fit_tree(useful, tree)
    results[[i]] = aim_make_matching_table(res, tree$ID_Arbre)
    
    nfit = 0
    if (length(results[[i]]$treeID) > 1L || !is.na(results[[i]]$treeID)) {
      nfit = length(res)
    }
    results[[i]]$nfit = nfit
    
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # rbindlist is significantly faster for data.frames/data.tables
  if (requireNamespace("data.table", quietly = TRUE)) {
    out <- as.data.frame(data.table::rbindlist(results))
  } else {
    out <- do.call(rbind, results)
  }
  
  class(out) <- c("aim_matching_table", class(out))
  return(out)
}

aim_make_matching_table = function(res, target_id)
{
  to_be_merge_id = lapply(res, function(x) x$to_be_merge_id)
  to_be_merge_id = do.call(c, to_be_merge_id)
  to_be_merge_id = unique(to_be_merge_id)

  if (length(to_be_merge_id) > 0)
    return(data.frame(treeID = to_be_merge_id, targetID = target_id))
  else
    return(data.frame(treeID = NA, targetID = target_id))
}

#' Reassign Tree IDs Using a Matching Table
#'
#' Reassigns tree IDs in a vector according to a matching table.
#'
#' @param las A `LAS` object instance-segmented with the arbor package.
#' @param matching_table A data frame containing treeID and targetID
#' columns defining the ID reassignment.
#'
#' @return A vector with tree IDs reassigned according to matching_table.
#'
#' @export
aim_reassign_ids <- function(las, matching_table)
{
  v = las$treeID
  
  # Check v
  if (!is.numeric(v) || !is.atomic(v) || is.object(v))
    stop("las$treeID must be a numeric vector.")
  if (any(!is.na(v) & !is.finite(v)))
    stop("las$treeID must contain only finite values or NA.")
  if (!is.integer(v))
    stop("las$treeID must be of integer type")
  
  # Check matching_table
  if (!is.data.frame(matching_table))
    stop("matching_table must be a data.frame.")
  required <- c("treeID", "targetID")
  missing <- setdiff(required, names(matching_table))
  if (length(missing) > 0)
    stop("matching_table is missing required column(s): ", paste(missing, collapse = ", "))
  
  # Check ID columns
  for (name in required)
  {
    x <- matching_table[[name]]
    if (!is.numeric(x))
      stop("matching_table$", name, " must be numeric.")
    if (any(!is.na(x) & !is.finite(x)))
      stop("matching_table$", name, " must contain only finite values or NA.")
    if (any(!is.na(x) & x != floor(x)))
      stop("matching_table$", name, " must contain integer values.")
  }
  
  if (anyNA(matching_table$targetID))
    stop("matching_table$targetID must not contain NA values.")
  if (anyNA(matching_table$treeID))
    matching_table <- matching_table[!is.na(matching_table$treeID), , drop = FALSE]
  
  if (length(v) == 0 || nrow(matching_table) == 0)
    return(v)
  
  # treeID must be unique: each source maps to exactly one target
  if (anyDuplicated(matching_table$treeID))
  {
    warning("Matching table contains duplicated 'treeID': use aim_resolve_multi_matching().", call. = FALSE)
    matching_table <- matching_table[
      !duplicated(matching_table$treeID) &
        !duplicated(matching_table$treeID, fromLast = TRUE),
      ,
      drop = FALSE
    ]
  }
  
  # Identity rows (e.g. 54 -> 54) change nothing. They must NOT be treated as
  # real rules: otherwise the ID is considered "handled" by the temp-routing
  # and is not protected from the points reassigned to it (e.g. 105 -> 54),
  # which would silently merge the two trees.
  matching_table <- matching_table[matching_table$treeID != matching_table$targetID, , drop = FALSE]
  
  if (nrow(matching_table) == 0)
    return(v)
  
  src_ids <- as.integer(matching_table$treeID)
  tgt_ids <- as.integer(matching_table$targetID)
  n_src   <- length(src_ids)
  
  # Identify pre-existing values in v that would collide with a targetID.
  # A target collides if it is NOT itself a source: a source ID is moved to a
  # temporary ID in Phase 1 and leaves its original value, so it frees the slot.
  # Identity rows were removed above, so a target that is also a source is
  # guaranteed to be moved elsewhere.
  colliding_targets <- setdiff(unique(tgt_ids), unique(src_ids))
  
  # Only the VALUES actually present in v matter, and each distinct value must
  # map to one single new id (all its occurrences must stay together as a group).
  bump_values <- intersect(unique(v), colliding_targets)
  n_bump      <- length(bump_values)
  
  # Build collision-free temporary / replacement IDs.
  # Must exceed EVERYTHING that could appear: v, all treeIDs, all targetIDs.
  max_val  <- max(c(v, src_ids, tgt_ids), na.rm = TRUE)
  new_pool <- max_val + seq_len(n_src + n_bump)
  temp_ids <- new_pool[seq_len(n_src)]            # one temp slot per mapping row
  bump_ids <- new_pool[n_src + seq_len(n_bump)]   # one fresh id per bumped VALUE
  
  # Phase 0: move pre-existing colliding values out of the way first
  if (n_bump > 0)
  {
    idx0 <- match(v, bump_values)
    hit0 <- !is.na(idx0)
    v[hit0] <- bump_ids[idx0[hit0]]
  }
  
  # Phase 1: source -> temporary
  idx1 <- match(v, src_ids)
  hit1 <- !is.na(idx1)
  v[hit1] <- temp_ids[idx1[hit1]]
  
  # Phase 2: temporary -> target
  idx2 <- match(v, temp_ids)
  hit2 <- !is.na(idx2)
  v[hit2] <- tgt_ids[idx2[hit2]]
  
  return(v)
}

#' Re-segment Tree IDs Using a Matching Table and an Inventory
#'
#' Re-segment Tree IDs Using a Matching Table and an Inventory.
#'
#' @param las A `LAS` object instance-segmented with the arbor package.
#' @param matching_table A data frame containing treeID and targetID
#' columns defining the ID reassignment.
#' @param inventory An sf data.frame produced by \link{aim_inventory}.
#'
#' @return A vector with tree IDs reassigned according to matching_table.
#'
#' @export
aim_resolve_multi_matching = function(las, matching_table, inventory)
{
  treeID <- NULL
  
  cnt    <- table(matching_table$treeID)
  dups   <- which(cnt > 1)
  dupIDs <- as.integer(names(dups))
  dupTargetID <- matching_table$targetID[matching_table$treeID %in% dupIDs]
  dubInventory <- inventory[inventory$ID_Arbre %in% dupTargetID,]

  #x = arbor::plot_instance(lidR::filter_poi(las, treeID %in% dupIDs))
  #aim_add_inventory3d(x, dubInventory)

  seeds = lapply(1:nrow(dubInventory), function(i)
  {
    res = aim_fit_tree(las, dubInventory[i,])
    seeds = lapply(res, function(x) x$seeds)
    seeds = do.call(rbind, seeds)
    if (is.matrix(seeds)) seeds = as.data.frame(seeds)
    seeds$treeID = dubInventory[i,]$ID_Arbre
    seeds
  })
  seeds = do.call(rbind, seeds)
  header = lidR::header(las)
  header@VLR = list()
  seeds = lidR::LAS(seeds, header)
  #lidR::plot(seeds, add = x, color = "treeID")

  las@data$PID = 1:lidR::npoints(las)
  dupTrees = lidR::filter_poi(las, treeID %in% dupIDs)

  dupTreesFixed = arbor::segment_instance(dupTrees, seeds)

  #x = arbor::plot_instance(dupTreesFixed)
  #aim_add_inventory3d(x, dubInventory)

  las$treeID[dupTreesFixed$PID] = dupTreesFixed$treeID

  #inventory$ID_Arbre %in% las$treeID

  las$treeID
}

#' @method print aim_matching_table
#' @export
print.aim_matching_table <- function(x, ...)
{
  stopifnot(
    is.data.frame(x),
    all(c("treeID", "targetID") %in% names(x))
  )

  mappings <- split(x$targetID, x$treeID)

  conflicts <- mappings[
    vapply(mappings, function(z) length(unique(z)) > 1L, logical(1))
  ]

  valid <- mappings[
    !vapply(mappings, function(z) length(unique(z)) > 1L, logical(1))
  ]

  cat(
    "Matching table: ",
    nrow(x), " mappings | ",
    length(mappings), " treeIDs",
    if (length(conflicts))
      paste0(" | \u26a0 ", length(conflicts), " conflicts"),
    "\n\n",
    sep = ""
  )

  if (length(valid))
  {
    txt <- vapply(
      names(valid),
      function(id)
        paste0(id, "\u2192", unique(valid[[id]])[1L]),
      character(1)
    )

    # Wrap according to console width
    cat(
      paste(
        strwrap(
          paste(txt, collapse = "  "),
          width = getOption("width") - 1
        ),
        collapse = "\n"
      ),
      "\n",
      sep = ""
    )
  }

  if (length(conflicts))
  {
    cat("\n\u26a0 Conflicts: ")

    txt <- vapply(
      names(conflicts),
      function(id)
        paste0(
          id, "\u2192",
          paste(unique(conflicts[[id]]), collapse = ",")
        ),
      character(1)
    )

    cat(
      paste(
        strwrap(
          paste(txt, collapse = "  "),
          width = getOption("width") - 15,
          initial = "",
          prefix = "              "
        ),
        collapse = "\n"
      ),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

