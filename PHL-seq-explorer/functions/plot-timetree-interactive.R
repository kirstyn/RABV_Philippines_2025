plot_timetree_interactive <- function(metadata, region_palette, tree_path = NULL,
                                      n_tip_sample = NULL, colour_by = "Region_short") {
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(plotly)
  library(ape)
  
  # 1. Load tree
  if (is.null(tree_path)) {
      tree_path <- here::here("processed_data", "trees", "timetree.nexus")  }
  timetree <- read.nexus(tree_path)
  
  # 2. Filter metadata to only tips in tree
  metadata <- metadata %>%
    filter(seq_name %in% timetree$tip.label) %>%
    mutate(Preferred_date = as.Date(Preferred_date))
  
  # 3. Optional subset of tips
  if (!is.null(n_tip_sample)) {
    keep_tips <- head(intersect(timetree$tip.label, metadata$seq_name), n_tip_sample)
    metadata <- metadata %>% filter(seq_name %in% keep_tips)
    timetree <- drop.tip(timetree, setdiff(timetree$tip.label, keep_tips))
  }
  
  # 4. Prepare dummy plot to extract coordinates
  tree_data <- ggtree(timetree)$data
  metat <- tree_data %>%
    inner_join(metadata, by = c("label" = "seq_name"))
  
  # 5. Fix factor levels to match the shared palette
  # (Use palette names directly from region_palette argument)
  if (!is.null(region_palette) && length(region_palette) > 0) {
    # limit levels to those present in palette
    if (colour_by %in% names(metadata)) {
      metat[[colour_by]] <- factor(metat[[colour_by]],
                                   levels = names(region_palette))
    }
  }
  
  # 6. Define colour scale using the provided palette for ALL variables
  colour_scale <- scale_colour_manual(
    values = region_palette,
    na.value = alpha("grey70", 0.7)
  )
  
  # 7. Build base tree plot
  dummy_plot <- ggtree(timetree)
  range_x <- range(dummy_plot$data$x)
  max_date <- max(metadata$Preferred_date, na.rm = TRUE)
  
  p <- ggplot(metat, aes(x = x, y = y, colour = .data[[colour_by]])) +
    geom_tree(data = ggtree(timetree)$data, aes(x = x, y = y),
              colour = "grey70", linewidth = 0.3) +
    geom_point(
      aes(
        text = paste0(
          "<b>Sample:</b> ", label, "<br>",
          "<b>Date:</b> ", Preferred_date, "<br>",
          Barangay, "-", Province_std, "-", Region_short
        )
      ),
      size = 2
    ) +
    theme_tree2() +
    colour_scale +
    labs(colour = colour_by) +
    scale_x_continuous(
      name = "Year",
      breaks = seq(floor(range_x[1]), ceiling(range_x[2]), by = 5),
      labels = rev(round(as.numeric(format(max_date, "%Y")) -
                           seq(floor(range_x[1]), ceiling(range_x[2]), by = 5)))
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 8. Remove interactive segment layers
  p$layers <- lapply(p$layers, function(layer) {
    if (inherits(layer$geom, "GeomInteractiveSegment")) {
      layer$geom <- ggplot2::GeomSegment
    }
    layer
  })
  
  # 9. Convert to interactive plotly
  ggplotly(p, tooltip = "text") %>%
    layout(legend = list(title = list(text = colour_by)))
}
