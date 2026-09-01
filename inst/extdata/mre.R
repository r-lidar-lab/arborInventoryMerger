library(lidR)
library(arbor)
library(arborInventoryMerger)
library(ggplot2)

# Path
finventory <- "Mbalmayo_Cameroun/inventory_Mbal09.txt"
fscene     <- "quadrat_4_buf15m_segmented.laz"
fdtm       <- "quadrat_4_buf15m_dtm.tif"

# Load data
data <- data.table::fread(file.path(dir, finventory))
las  <- readTLS(file.path(dir, fscene))
dtm  <- terra::rast(file.path(dir, fdtm))

# Convert inventory data to spatial data with internal validation
hull      <- sf::st_convex_hull(lidR::filter_poi(las, hag > 0.5, hag < 2))
inventory <- aim_inventory(data, dtm)
inventory <- sf::st_filter(inventory, hull)

# Quick 3D inventory plot
gnd <- lidR::filter_poi(las, hag > 0.25 & hag < 4)
x   <- arbor::plot_instance(gnd)
aim_add_inventory3d(x, inventory)

# Quick 2D inventory plot
ggplot(inventory, aes(color = as.factor(ID_Arbre))) + geom_sf() +
  geom_sf(data = hull, inherit.aes = FALSE, alpha = 0.5) +
  guides(color="none") +
  theme_bw()

# Quick visualization of what is happening for tree 79
tree <- dplyr::filter(inventory, ID_Arbre == 79)
res  <- arborInventoryMerger:::aim_fit_tree(las, tree)
arborInventoryMerger:::aim_view(las, tree, res)

# Merge inventory and segmented point cloud
matching_table <- aim_matching_table(las, inventory)
las$treeID <- aim_reassign_ids(las, matching_table)
las$treeID <- aim_resolve_multi_matching(las, matching_table, inventory)

# Render low points
gnd <- lidR::filter_poi(las, hag < 4)
x   <- arbor::plot_instance(gnd)
aim_add_inventory3d(x, inventory)

# Extract and render inventory trees
lasinventory <- lidR::filter_poi(las, treeID %in% inventory$ID_Arbre)
x <- arbor::plot_instance(lasinventory)
aim_add_inventory3d(x, inventory)
lidR::add_dtm3d(x, dtm)
