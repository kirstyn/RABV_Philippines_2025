plot_timetree_interactive <- function(metadata, region_palette, tree_path = NULL, n_tip_sample = NULL) {
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(plotly)
  library(ape)
  library(here)
  library(ggraph)

  # Load tree
  if (is.null(tree_path)) {
    tree_path <- here::here("processed_data", "trees", "timetree.nexus")
  }
  timetree <- read.nexus(tree_path)

  # Filter metadata to only tips in tree
  metadata <- metadata %>%
    filter(seq_name %in% timetree$tip.label) %>%
    mutate(
      Preferred_date = as.Date(Preferred_date),
      Region_short = factor(Region_short, levels = names(region_palette))
    )

  # Optional subset of tips
  if (!is.null(n_tip_sample)) {
    keep_tips <- head(intersect(timetree$tip.label, metadata$seq_name), n_tip_sample)
    metadata <- metadata %>% filter(seq_name %in% keep_tips)
    timetree <- drop.tip(timetree, setdiff(timetree$tip.label, keep_tips))
  }

  # Dummy plot to get x-range
  dummy_plot <- ggtree(timetree)
  range_x <- range(dummy_plot$data$x)
  max_date <- max(metadata$Preferred_date, na.rm = TRUE)
  
  # 5. Merge metadata with tree tips
  metat <- dummy_plot$data %>%
    dplyr::inner_join(metadata, by = c("label" = "seq_name"))
  
  # 4. Merge metadata with tree tip positions
  tree_data <- ggtree(timetree)$data
  metat <- tree_data %>%
    inner_join(metadata, by = c("label" = "seq_name"))

  # Build ggtree plot
  p <- ggtree(timetree) +
    geom_point(
      data = metat,
      aes(
        x = x,
        y = y,
        colour = Region_short,
        text = paste0(
          "<b>Sample:</b> ", label, "<br>",
          "<b>Date:</b> ", Preferred_date, "<br>",
          Barangay,"-",Province_std, "-",Region_short
        )
      ),
      size = 2
    ) +
    scale_x_continuous(
        name = "Year",
        breaks = seq(floor(range_x[1]), ceiling(range_x[2]), by = 5),
        labels = rev(round(as.numeric(format(max_date, "%Y")) -
                             seq(floor(range_x[1]), ceiling(range_x[2]), by = 5)))
      ) +
    scale_colour_manual(values = region_palette, na.value = alpha("grey70",0.7)) +
    theme_tree2()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

  # Remove interactive segment layers
  p$layers <- lapply(p$layers, function(layer) {
   if (inherits(layer$geom, "GeomInteractiveSegment")) {
      layer$geom <- ggplot2::GeomSegment
   }
    layer
  })

  # Convert to interactive plotly plot
  # 6. Convert to plotly interactive plot
  ggplotly(p, tooltip = "text") %>%
    layout(showlegend = FALSE)
}
