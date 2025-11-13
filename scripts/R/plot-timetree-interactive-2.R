plot_timetree_interactive <- function(metadata, region_palette, tree_path = NULL, n_tip_sample = NULL) {
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(here)
  library(plotly)
  
  # 1. Load tree
  if (is.null(tree_path)) {
    tree_path <- here("analysis", "Timetree", "PHL_n786", "timetree.nexus")
  }
  timetree <- read.nexus(tree_path)
  
  # 2. Filter metadata to only tips in tree
  metadata <- metadata %>%
    filter(seq_name %in% timetree$tip.label) %>%
    mutate(
      Preferred_date = as.Date(Preferred_date, format = "%d-%b-%Y"),
      Region_short = factor(Region_short, levels = names(region_palette))
    ) %>%
    select(seq_name, everything())
  
  # 3. Optionally subset to first n tips
  if (!is.null(n_tip_sample)) {
    keep_tips <- head(intersect(timetree$tip.label, metadata$seq_name), n_tip_sample)
    metadata <- metadata %>% filter(seq_name %in% keep_tips)
    timetree <- ape::drop.tip(timetree, setdiff(timetree$tip.label, keep_tips))
  }
  
  # 4. Merge metadata with tree tips using ggtree 4.x helper
  ggtree_data <- as_tibble(timetree) %>%
    left_join(metadata, by = c("label" = "seq_name"))
  
  # 5. Build ggtree plot
  p <- ggtree(timetree, branch.length = "none") %<+% ggtree_data +
    geom_tippoint(aes(colour = Region_short), size = 2, alpha = 0.7) +
    scale_colour_manual(values = region_palette, na.value = alpha("grey", 0.5)) +
    theme_tree2() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # 6. Convert to interactive Plotly object
  ggplotly(p, tooltip = c("label", "Region_short"))
}