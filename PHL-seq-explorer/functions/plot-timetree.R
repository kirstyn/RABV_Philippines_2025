plot_timetree <- function(metadata, region_palette) {
    library(pacman)
   pacman::p_load(
      ggplot2,
      dplyr,
      here
    )

    # Then handle Bioconductor packages separately
    for (pkg in c("ggtree", "treeio")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = TRUE)
      }
      library(pkg, character.only = TRUE)
    }
  # Load tree
  timetree_path <- here("processed_data", "trees", "timetree.nexus")
  timetree <- read.nexus(timetree_path)
  
  # Ensure seq_name is first column
  metadata <- metadata %>% select(seq_name, everything())
  
  # Convert dates
  metadata$Preferred_date <- as.Date(metadata$Preferred_date, format = "%d-%b-%Y")
  max_date <- max(metadata$Preferred_date, na.rm = TRUE)
  
  # Ensure Region_short exists and is factor with consistent levels
  metadata$Region_short <- factor(metadata$Region_short, levels = names(region_palette))
  
  # Filter metadata to only tips present in tree
  metadata <- metadata %>% filter(seq_name %in% timetree$tip.label)
  
  # Dummy plot to extract x range
  dummy_plot <- ggtree(timetree)
  range_x <- range(dummy_plot$data$x)
  
  # Plot tree
  p <- ggtree(timetree) %<+% metadata +
    geom_tippoint(aes(colour = Region_short), size = 2) +
    theme_tree2() +
    scale_x_continuous(
      name = "Year",
      breaks = seq(floor(range_x[1]), ceiling(range_x[2]), by = 5),
      labels = rev(round(as.numeric(format(max_date, "%Y")) -
                           seq(floor(range_x[1]), ceiling(range_x[2]), by = 5)))
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_colour_manual(values = region_palette,na.value = alpha("grey", 0.5))
  
  return(p)
}
