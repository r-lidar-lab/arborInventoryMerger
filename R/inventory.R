#' Convert a Field Inventory to an XYZ sf Object
#'
#' Validates and standardizes a field inventory data frame, converts measurements to
#' metres when required, extracts ground elevations from a Digital Terrain Model (DTM),
#' and generates a 3D (`XYZ`) `sf` point object representing tree stem positions.
#'
#' @details
#' The function automatically detects and standardizes non-standard column names
#' using case-insensitive regular expression matching. Columns are internally mapped to
#' canonical names (`ID_Arbre`, `X_correg`, `Y_correg`, `DHP`, `POM`, `Z`).\cr\cr
#' Validation checks performed:
#' \itemize{
#'   \item \strong{Unit Conversion}: Measurements (`DHP`, `POM`) are checked for scale.
#'         Values consistently exceeding unit thresholds are converted from centimetres to metres.
#'   \item \strong{Spatial Extraction & Elevation Handling}:
#'         If \code{Z} is not present in the input table, elevation (\code{Z}) is derived by sampling
#'         the DTM raster at \code{(X_correg, Y_correg)} coordinates and adding \code{POM}.
#'         If \code{Z} is present in the input table, \code{Z} is kept as is and \code{POM} is calculated as \code{Z - DTM}.
#'         Rows with missing spatial or raster coverage are discarded with a warning.
#'   \item \strong{Topology Check}: Overlapping tree stem buffers (calculated from \code{DHP})
#'         are evaluated and flagged if substantial overlap exists between tree stems.
#' }
#'
#' @param data A data frame containing tree inventory records. The following columns
#'   and their authorized variations (case-insensitive) are accepted:
#'   \describe{
#'     \item{\code{ID_Arbre}}{Tree unique identifier. Authorized names: \code{ID_Arbre}, \code{Tree_ID}, \code{TreeID}, \code{ID}, \code{IDtree}.}
#'     \item{\code{X_correg}}{X coordinate (metres). Authorized names: \code{X_correg}, \code{X_corrected}, \code{X_coord}, \code{X}, \code{CenterX}.}
#'     \item{\code{Y_correg}}{Y coordinate (metres). Authorized names: \code{Y_correg}, \code{Y_corrected}, \code{Y_coord}, \code{Y}, \code{CenterY}.}
#'     \item{\code{DHP}}{Diameter at breast height. Authorized names: \code{DHP}, \code{DBH}.}
#'     \item{\code{POM}}{Point of measurement height. Authorized names: \code{POM}, \code{HT_POM}, \code{Point_of_measurement}. (Optional if \code{Z} is provided).}
#'     \item{\code{Z}}{Absolute elevation above sea level. Authorized names: \code{Z}, \code{Z_coord}, \code{Elevation}, \code{Alt}, \code{Altitude}. (Optional if \code{POM} is provided).}
#'   }
#' @param dtm A \code{\link[terra]{SpatRaster}} object representing the Digital Terrain Model (DTM) used to sample ground surface elevation.
#'
#' @return A 3D simple feature object (\code{\link[sf]{sf}}) with \code{POINT} geometry
#'   and \code{XYZ} dimensions containing standardized attribute names (\code{ID_Arbre},
#'   \code{X_correg}, \code{Y_correg}, \code{DHP}, \code{POM}, \code{Z}).
#'
#' @importFrom terra extract
#' @importFrom sf st_as_sf st_buffer st_intersection st_area st_agr<-
#'
#' @export
aim_inventory <- function(data, dtm)
{
  fix_names_encoding <- function(x)
  {
    stopifnot(is.character(x))
    
    # Keep valid UTF-8 strings unchanged
    ok <- !is.na(iconv(x, from = "UTF-8", to = "UTF-8", sub = NA))
    
    # Replace invalid strings with a safe ASCII representation
    x[!ok] <- iconv(x[!ok], from = "", to = "ASCII//TRANSLIT", sub = "_")
    
    x
  }
  
  names(data) <- fix_names_encoding(names(data))
  
  # 1. Column Detection and Standardization
  
  column_specs <- list(
    ID_Arbre = "^(ID_Arbre|Tree_ID|TreeID|ID|IDtree)$",
    X_correg = "^(X_correg|X_corrected|X_coord|X|CenterX)$",
    Y_correg = "^(Y_correg|Y_corrected|Y_coord|Y|CenterY)$",
    DHP      = "^(DHP|DBH)$",
    POM      = "^(POM|HT_POM|Point_of_measurement)$",
    Z        = "^(Z|Z_coord|Elevation|Alt|Altitude)$"
  )
  
  detected_cols <- character(0)
  required_cols <- c("ID_Arbre", "X_correg", "Y_correg", "DHP")
  missing_cols  <- c()
  
  for (canonical_name in names(column_specs))
  {
    pattern <- column_specs[[canonical_name]]
    match   <- suppressWarnings(grep(pattern, names(data), ignore.case = TRUE, value = TRUE))
    
    if (length(match) == 1)
    {
      detected_cols[canonical_name] <- match
    }
    else if (length(match) > 1)
    {
      stop(sprintf("Ambiguous match for '%s': found multiple columns matching pattern: %s", canonical_name, paste(match, collapse = ", ")))
    }
    else if (canonical_name %in% required_cols)
    {
      missing_cols <- c(missing_cols, canonical_name)
    }
  }
  
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  has_z   <- "Z" %in% names(detected_cols)
  has_pom <- "POM" %in% names(detected_cols)
  
  if (!has_z && !has_pom) {
    stop("Input data must contain at least one height/elevation column: 'POM' or 'Z'.")
  }
  
  message("Column Detection Log")
  for (canonical_name in names(detected_cols))
  {
    message(sprintf("  Canonical: %-10s -> Detected: '%s'", canonical_name, detected_cols[canonical_name]))
  }
  
  for (canonical_name in names(detected_cols))
  {
    orig_name <- detected_cols[canonical_name]
    if (orig_name != canonical_name)
    {
      names(data)[names(data) == orig_name] <- canonical_name
    }
  }
  
  
  # 2. ID Validation
  
  id <- data$ID_Arbre
  
  if (!is.numeric(id))
  {
    stop(paste0("ID_Arbre must contain numbers not ", typeof(id)))
  }
  
  if (any(!is.na(id) & (!is.finite(id) | id != floor(id))))
  {
    stop("ID_Arbre must contain integral numbers.")
  }
  data$ID_Arbre <- as.integer(id)
  
  
  # 3. Measurement Validation & Unit Conversion Log
  
  message("Data Units & Quality Log")
  
  check_measurement <- function(x, name, th)
  {
    if (!is.numeric(x)) {
      stop(name, " must be numeric.")
    }
    
    has_issues <- FALSE
    
    # Check non-finite values
    non_finite_idx <- which(!is.finite(x))
    if (length(non_finite_idx) > 0)
    {
      message(sprintf("  [Issue] %s contains %d non-finite values (NA/NaN/Inf) at indices: %s",
                      name, length(non_finite_idx), paste(utils::head(non_finite_idx, 5), collapse = ", ")))
      has_issues <- TRUE
    }
    
    x_valid <- x[is.finite(x)]
    
    if (length(x_valid) == 0)
    {
      message(sprintf("  [Critical] %s contains NO valid numeric values.", name))
      warning(sprintf("%s contains no valid values. Check input data.", name), call. = FALSE)
      return(x)
    }
    
    # Unit Conversion checks
    if (all(x_valid > th))
    {
      message(sprintf("  [Info] %s values all > %g. Converting from cm to metres (dividing by 100).", name, th))
      x <- x / 100
      x_valid <- x[is.finite(x)]
    }
    else if (any(x_valid > th))
    {
      suspect_count <- sum(x_valid > th)
      message(sprintf("  [Issue] %s contains %d value(s) > %g m. Mixed units suspected; no conversion applied.", name, suspect_count, th))
      has_issues <- TRUE
    }
    
    # Negative check
    neg_count <- sum(x_valid < 0)
    if (neg_count > 0)
    {
      message(sprintf("  [Issue] %s contains %d negative value(s).", name, neg_count))
      has_issues <- TRUE
    }
    
    if (has_issues)
    {
      warning(sprintf("Quality issues detected in '%s'. See log details above.", name), call. = FALSE)
    }
    
    return(x)
  }
  
  data$DHP <- check_measurement(data$DHP, "DHP", 1)
  
  if (!has_z)
  {
    data$POM <- check_measurement(data$POM, "POM", 10)
  }
  else
  {
    if (!is.numeric(data$Z))
    {
      stop("Z must be numeric.")
    }
    if (any(!is.finite(data$Z)))
    {
      warning("Z contains non-finite values (NA/NaN/Inf).", call. = FALSE)
    }
  }
  
  
  # 4. Spatial Extraction & Coordinate Validation
  
  if (!is.numeric(data$X_correg) || !is.numeric(data$Y_correg))
  {
    stop("X_correg and Y_correg must be numeric.")
  }
  
  if (any(!is.finite(data$X_correg)) || any(!is.finite(data$Y_correg)))
  {
    stop("X_correg and Y_correg must contain only finite values.")
  }
  
  dtm_vals <- terra::extract(dtm, cbind(data$X_correg, data$Y_correg))[[1]]
  
  if (has_z)
  {
    message("  [Info] 'Z' detected in input table. Keeping Z as is and computing POM = Z - DTM elevation.")
    data$POM <- data$Z - dtm_vals
    
    neg_pom_count <- sum(data$POM < 0, na.rm = TRUE)
    if (neg_pom_count > 0)
    {
      message(sprintf("  [Issue] Calculated POM is negative for %d tree(s) (Z is lower than DTM surface).", neg_pom_count))
      warning(sprintf("%d tree(s) have calculated POM < 0 m.", neg_pom_count), call. = FALSE)
    }
  }
  else
  {
    data$Z <- dtm_vals + data$POM
  }
  
  na_spatial_count <- sum(is.na(data$Z) | is.na(data$POM))
  if (na_spatial_count > 0)
  {
    message(sprintf("  [Warning] %d trees outside DTM extent or missing spatial/elevation data. These rows were removed.", na_spatial_count))
    warning(sprintf("%d observations were removed due to missing DTM elevation or height values.", na_spatial_count), call. = FALSE)
    data <- data[!is.na(data$Z) & !is.na(data$POM), ]
  }
  
  inventory <- sf::st_as_sf(
    data,
    coords = c("X_correg", "Y_correg", "Z"),
    dim = "XYZ"
  )
  
  
  # 5. Topology / Overlap Check
  
  message("Inventory Topology Check")
  
  circles <- sf::st_buffer(inventory, inventory$DHP / 2)
  sf::st_agr(circles) <- "constant"
  
  intersections <- sf::st_intersection(circles, circles)
  intersections <- intersections[intersections$ID_Arbre != intersections$ID_Arbre.1, ]
  
  if (nrow(intersections) > 0)
  {
    intersections <- intersections[intersections$ID_Arbre < intersections$ID_Arbre.1, ]
    
    area_inter <- sf::st_area(intersections)
    area_1      <- sf::st_area(circles[match(intersections$ID_Arbre, circles$ID_Arbre), ])
    area_2      <- sf::st_area(circles[match(intersections$ID_Arbre.1, circles$ID_Arbre), ])
    
    min_area    <- pmin(area_1, area_2)
    overlap_pct <- as.numeric((area_inter / min_area) * 100)
    
    moderate_or_major_issues <- 0
    
    for (k in seq_along(overlap_pct))
    {
      id_a <- intersections$ID_Arbre[k]
      id_b <- intersections$ID_Arbre.1[k]
      pct  <- round(overlap_pct[k], 2)
      
      if (pct >= 5 && pct <= 80)
      {
        moderate_or_major_issues <- moderate_or_major_issues + 1
        message(sprintf("  [Warning] Circles %s and %s overlap by %s%%", id_a, id_b, pct))
      }
      else if (pct > 80)
      {
        moderate_or_major_issues <- moderate_or_major_issues + 1
        message(sprintf("  [Critical] Circles %s and %s overlap by %s%%", id_a, id_b, pct))
      }
      else if (pct >= 1)
      {
        message(sprintf("  [Message] Circles %s and %s overlap by %s%%", id_a, id_b, pct))
      }
    }
    
    if (moderate_or_major_issues > 0)
    {
      warning(sprintf("%d moderate/critical topological overlap issues detected in inventory. See logs above.", moderate_or_major_issues), call. = FALSE)
    }
    else
    {
      message("  No significant topological overlaps detected.")
    }
  }
  else
  {
    message("  No overlapping geometries found.")
  }
  
  return(inventory)
}