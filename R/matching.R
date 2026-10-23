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
#' @param inventory An \code{sf} object produced by \link{aim_inventory}.
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
  
  if (!any(unique(las$treeID) %in% inventory$ID_Arbre)) {
    message("No matching treeIDs found between LAS data and inventory.")
    return(las)
  }

  las$treeID <- arbor:::aim_disentangle_ids_cpp(las@data, aim_inventory_df(inventory))
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
  
  out <- arbor:::aim_matching_table_cpp(las@data, aim_inventory_df(inventory))
  class(out) <- c("aim_matching_table", class(out))
  return(out)
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
  
  # Duplicated treeIDs are left untouched (use aim_resolve_multi_matching()),
  # identity rows are ignored and pre-existing IDs colliding with a target are moved
  arbor:::aim_reassign_ids_cpp(las@data, matching_table)
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
  arbor:::aim_resolve_multi_matching_cpp(las@data, matching_table, aim_inventory_df(inventory))
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

