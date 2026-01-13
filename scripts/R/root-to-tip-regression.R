# ---- Load packages ----
library(ape)
library(ggplot2)
library(dplyr)
library(lubridate)
library(phytools)
library(adephylo)
library(ggtree)

# ---- 1. Input files ----
tree_file <- "analysis/ML_trees/241025_PHL_all_n794_withSEA2outgroups.ft.newick"
meta_file <- "processed_data/processed_metadata/gathered_metadata/final/PHL_metadata_matchedToSeq_2026-01-05_n794.csv"

# ---- 2. Read in tree and metadata ----
tree <- read.tree(tree_file)
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Trim whitespace from tree tip labels
tree$tip.label <- str_trim(tree$tip.label)

# Trim whitespace from key metadata columns
meta <- meta %>%
  mutate(
    Sample_ID = str_trim(Sample_ID),
    Accession = str_trim(Accession),
    Source = str_trim(Source)
  )
# ---- 3. Parse dates ----
meta <- meta %>%
  mutate(
    Preferred_date = parse_date_time(Preferred_date, orders ="%d-%b-%Y"),
    decimal_date = decimal_date(Preferred_date)
  )


# Optional: check how many have matches
cat(sum(!is.na(meta$seq_name)), "metadata entries matched to tree tips\n")
# ---- Check which metadata entries didn't match tree tips ----
unmatched_meta_tree <- meta %>%
  filter(is.na(seq_name)) %>%
  select(Sample_ID, seq_name,Accession, Source)

# View summary
cat(nrow(unmatched_meta_tree), "metadata entries not matched to tree tips\n")

# Optional: inspect them
unmatched_meta_tree

# ---- 5. Root tree ----
# Root using outgroup (specify tip labels of outgroup)
tree$tip.label[tree$tip.label == "EU643590.1"] <- "EU643590_SEA2a"
# Identify tips containing "_SEA2a" or "_SEA2b"
outgroup_tips <- tree$tip.label[grep("_SEA2[ab]", tree$tip.label)]

# Root the tree using these tips
tree_rooted <- root(tree, outgroup = outgroup_tips, resolve.root = TRUE)

# # Root using midpoint (less preferred)
# tree_rooted <- midpoint.root(tree)
# Drop the outgroup tips from the tree
tree_final <- drop.tip(tree_rooted, outgroup_tips)

ggtree(tree_final) +
  geom_tiplab(size = 3) +
  theme_tree2() +
  ggtitle("Rooted Phylogenetic Tree (Outgroup Removed)")


# ---- 6. Calculate root-to-tip distances ----
rtt_data <- data.frame(
  tip = tree_final$tip.label,
  distance = distRoot(tree_final)
)

# ---- 7. Merge with metadata ----
rtt_meta <- rtt_data %>%
  left_join(meta, by = c("tip" = "seq_name")) %>%
  filter(!is.na(decimal_date))
rtt_missing_dates <- rtt_data %>%
  left_join(meta, by = c("tip" = "seq_name")) %>%
  filter(is.na(decimal_date)) %>%
  select(tip, Sample_ID, Accession, Preferred_date, Source, Region_std)

# Summary
cat(nrow(rtt_missing_dates), "tips have no associated date\n")

# View them
rtt_missing_dates

# ---- 8. Root-to-tip regression ----
fit <- lm(distance ~ decimal_date, data = rtt_meta)
summary(fit)

# ---- Estimate root date ----
intercept <- coef(fit)[1]
slope     <- coef(fit)[2]

root_decimal_year <- -intercept / slope
root_decimal_year

# ---- 9. Plot root-to-tip regression ----
p <- ggplot(rtt_meta, aes(x = decimal_date, y = distance)) +
  geom_point(aes(color = Major_Island), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Root-to-tip Regression",
    x = "Sampling Date (Decimal Year)",
    y = "Root-to-tip Distance",
    color = "Major_Island"
  ) +
  annotate("text",
           x = min(rtt_meta$decimal_date, na.rm = TRUE),
           y = max(rtt_meta$distance, na.rm = TRUE),
           hjust = 0,
           label = paste0("R² = ", round(summary(fit)$r.squared, 3),
                          "\nRate = ", signif(coef(fit)[2], 3), " subs/site/year")
  )

print(p)

# ---- 10. Save outputs ----
write.tree(tree_rooted, file = "analysis/Temporal_signal/20251023_PHL_all_filtered_withSEA2outgroups.OutgpRooted.ft.newick")
ggsave("analysis/Temporal_signal/20251023_PHL_all_filtered_withSEA2outgroups_RTTplot.png", p, width = 8, height = 6, dpi = 300)
write.csv(rtt_meta, "analysis/Temporal_signal/27Oct25_gathered_metadata_n797_raddl_and_manual_Corrected_stdGeo_RTTmetadata.csv", row.names = FALSE)
write.csv(meta,  "processed_data/processed_metadata/gathered_metadata/27Oct25_gathered_metadata_n782_raddl_and_manual_Corrected_stdGeo_seqNames.csv", row.names = FALSE)

## Compare island groups
# Root-to-tip regression
rtt_meta %>%
  filter(!is.na(Major_Island)) %>%
  group_by(Major_Island) %>%
  summarise(
    n = n(),
    min_year = min(decimal_date, na.rm = TRUE),
    max_year = max(decimal_date, na.rm = TRUE),
    year_span = max_year - min_year,
    rate = coef(lm(distance ~ decimal_date))[2],
    R2 = summary(lm(distance ~ decimal_date))$r.squared,
    .groups = "drop"
  )


p_pub <- ggplot(
  rtt_meta %>% filter(!is.na(Major_Island)),
  aes(x = decimal_date, y = distance, colour = Major_Island)
) +
  geom_point(size = 1.8, alpha = 0.75) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1
  ) +
  scale_colour_viridis_d(option = "D", end = 0.9) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    plot.title   = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Root-to-tip Regression by Major Island Group",
    x = "Sampling date (decimal year)",
    y = "Root-to-tip genetic distance",
    colour = "Major island"
  )

print(p_pub)

ggsave(
  filename = "analysis/Temporal_signal/241025_PHL_all_n794_RTT_by_island_publication.svg",
  plot = p_pub,
  width = 180,
  height = 140,
  units = "mm"
)




# 1️⃣ Calculate regression stats per Region
reg_stats <- rtt_meta %>%
  filter(!is.na(Region_short)) %>%
  group_by(Region_short) %>%
  summarise(
    n         = n(),
    min_year  = min(decimal_date, na.rm = TRUE),
    max_year  = max(decimal_date, na.rm = TRUE),
    year_span = max_year - min_year,
    rate      = coef(lm(distance ~ decimal_date))[2],
    R2        = summary(lm(distance ~ decimal_date))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    # Label for top-left of each facet
    label = paste0(
      "\nRate = ", signif(rate, 3),
      "\nR² = ", round(R2, 3)
    )
  )

# 2️⃣ Faceted plot by Region with top-left labels
p_faceted_region <- ggplot(
  rtt_meta %>% filter(!is.na(Region_short)),
  aes(x = decimal_date, y = distance)
) +
  geom_point(size = 1.8, alpha = 0.75, colour = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, colour = "darkred", linewidth = 1) +
  facet_wrap(~ Region_short, scales = "free_y") +
  geom_text(
    data = reg_stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,  # slightly inset from left
    vjust = 0.8,   # slightly inset from top
    size = 3.5
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    plot.title   = element_text(size = 16, face = "bold")
  ) +
  labs(
    title = "Root-to-tip Regression by Region",
    x = "Sampling date (decimal year)",
    y = "Root-to-tip genetic distance"
  )

print(p_faceted_region)

# Save as high-resolution PNG
ggsave(
  filename = "analysis/Temporal_signal/241025_PHL_all_n794_RTT_by_Region.png",
  plot = p_faceted_region,
  width = 12,       # figure width in inches
  height = 8,       # figure height in inches
  dpi = 600         # high resolution for print
)
