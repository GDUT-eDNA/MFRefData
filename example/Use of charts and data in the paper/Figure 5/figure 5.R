REQUIRED_PACKAGES <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "stringr",
  "forcats",
  "scales",
  "tidyr",
  "grid",
  "officer",
  "rvg",
  "cowplot"
)

MISSING_PACKAGES <- REQUIRED_PACKAGES[
  !vapply(REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]

if (length(MISSING_PACKAGES) > 0) {
  stop(
    "Missing R packages: ",
    paste(MISSING_PACKAGES, collapse = ", "),
    "\nPlease install them first:\n",
    "install.packages(c(",
    paste(sprintf('\"%s\"', MISSING_PACKAGES), collapse = ", "),
    "))"
  )
}

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)
library(scales)
library(tidyr)
library(grid)
library(officer)
library(rvg)
library(cowplot)

INPUT_BUBBLE_DATA <- "D:/figure 5.csv"

OUT_DIR <- "D:/Figure 5"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PLOT_DATA <- file.path(OUT_DIR, "bubble_plot_data_5class_top_hits.csv")
OUT_SUMMARY_TABLE <- file.path(OUT_DIR, "bubble_upper_summary_5class.csv")

OUT_MAIN_PREFIX <- file.path(OUT_DIR, "bubble_5class_species_tool_primer")

FONT_FAMILY <- "Arial"

TEXT_COLOR <- "#000000"

if (.Platform$OS.type == "windows") {
  try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
}

OUTPUT_FORMATS <- c("png", "tiff", "pdf", "svg")

FONT_SIZE_AXIS_X <- 13.4
FONT_SIZE_AXIS_Y <- 14.4
FONT_SIZE_AXIS_TITLE_X <- 18
FONT_SIZE_AXIS_TITLE_Y <- 18
FONT_SIZE_TITLE <- 18
FONT_SIZE_SUBTITLE <- 12
FONT_SIZE_CAPTION <- 10

FONT_SIZE_X_MARK <- 4.0

LEGEND_TITLE_SIZE <- 4.1
LEGEND_TEXT_SIZE <- 3.45
LEGEND_SYMBOL_TEXT_SIZE <- 3.60

LEGEND_TOOLKIT_POINT_SIZE <- 4.8
LEGEND_PRIMER_POINT_SIZE <- 4.8
LEGEND_TOPHIT_SIZE_1 <- 1.4
LEGEND_TOPHIT_SIZE_10 <- 3.0
LEGEND_TOPHIT_SIZE_40 <- 5.2

LEGEND_ITEM_STEP <- 0.018
LEGEND_GROUP_GAP <- 0.08

LEGEND_TITLE_ITEM_GAP <- LEGEND_ITEM_STEP * 1.30

ANNOTATION_RESULT_TITLE <- "Annotation result"

ANNOTATION_RESULT_SYMBOLS <- c(
  "Marked",
  "Blank",
  "x",
  "xx",
  "xxx"
)

ANNOTATION_RESULT_LABELS <- c(
  "correct species",
  "unassigned",
  "same genus, different species",
  "same family, different genus/species",
  "beyond-family mismatch"
)

ANNOTATION_RESULT_WRAP <- TRUE

ANNOTATION_RESULT_WRAP_WIDTH <- 26

ANNOTATION_RESULT_SYMBOL_X <- 0.003
ANNOTATION_RESULT_TEXT_X <- 0.38

ANNOTATION_RESULT_LINEHEIGHT <- 0.92

ANNOTATION_RESULT_EXTRA_LINE_STEP <- LEGEND_ITEM_STEP * 0.65

TOOL_LEVELS <- c("MFRefData", "refdb", "CRABS")

TOOL_COLORS <- c(
  "MFRefData" = "#D73027",
  "refdb"     = "#1A9850",
  "CRABS"     = "#4575B4"
)

TOOL_BLOCK_COLORS <- c(
  "MFRefData" = "#FEE5D9",
  "refdb"     = "#E5F5E0",
  "CRABS"     = "#DEEBF7"
)

PRIMER_LEVELS <- c("Tele02", "Riaz-12S", "12S-V5")

PRIMER_SHAPES <- c(
  "Tele02"   = 24,
  "Riaz-12S" = 21,
  "12S-V5"   = 22
)

TOP_HIT_COUNT_CAP <- 40
SIZE_RANGE <- c(1.8, 8.0)

SIZE_BREAKS <- c(1, 10, 40)
SIZE_LABELS <- c("1", "10", "40")

HIDE_UNASSIGNED_POINTS <- TRUE
KEEP_ALL_SPECIES_ON_Y_AXIS <- TRUE

EXACT_BORDER_COLOR <- "black"
EXACT_POINT_STROKE <- 1.35

X_MARK_NUDGE_X <- 0
X_MARK_NUDGE_Y <- 0

CLASS_LABEL_EXACT <- paste0(ANNOTATION_RESULT_SYMBOLS[1], ": ", ANNOTATION_RESULT_LABELS[1])
CLASS_LABEL_UNASSIGNED <- paste0(ANNOTATION_RESULT_SYMBOLS[2], ": ", ANNOTATION_RESULT_LABELS[2])
CLASS_LABEL_SAME_GENUS <- paste0(ANNOTATION_RESULT_SYMBOLS[3], ": ", ANNOTATION_RESULT_LABELS[3])
CLASS_LABEL_SAME_FAMILY <- paste0(ANNOTATION_RESULT_SYMBOLS[4], ": ", ANNOTATION_RESULT_LABELS[4])
CLASS_LABEL_BEYOND_FAMILY <- paste0(ANNOTATION_RESULT_SYMBOLS[5], ": ", ANNOTATION_RESULT_LABELS[5])

CLASS_LABEL_LEVELS <- c(
  CLASS_LABEL_EXACT,
  CLASS_LABEL_UNASSIGNED,
  CLASS_LABEL_SAME_GENUS,
  CLASS_LABEL_SAME_FAMILY,
  CLASS_LABEL_BEYOND_FAMILY
)

MAIN_WIDTH <- 13
MAIN_HEIGHT <- 11
MAIN_DPI <- 800

FINAL_WIDTH <- 15.0
FINAL_HEIGHT <- 11.0
FINAL_DPI <- 800

MAIN_REL_WIDTH <- 1.00
LEGEND_REL_WIDTH <- 0.20

MAIN_TITLE <- "Mock Community BLAST Annotation Results"
MAIN_X_TITLE <- ""
MAIN_Y_TITLE <- ""

REVERSE_SPECIES_ORDER <- TRUE
X_LABEL_ANGLE <- 0

clean_taxon <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "na", "nan", "NULL", "null", "None", "none")] <- NA_character_
  x
}

split_taxon_list <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || trimws(x) == "") {
    return(character(0))
  }
  y <- unlist(strsplit(x, ";", fixed = TRUE))
  y <- clean_taxon(y)
  unique(y[!is.na(y) & y != ""])
}

count_semicolon_items <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  ifelse(x == "", 0L, stringr::str_count(x, ";") + 1L)
}

all_equal_to_target <- function(values, target) {
  values <- clean_taxon(values)
  target <- clean_taxon(target)
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0 || is.na(target) || target == "") {
    return(FALSE)
  }
  all(values == target)
}

make_tool_primer_position <- function(df) {
  df %>%
    mutate(
      tool = factor(tool, levels = TOOL_LEVELS),
      primer = factor(primer, levels = PRIMER_LEVELS),
      tool_order = as.integer(tool),
      primer_order = as.integer(primer),
      x_pos = (tool_order - 1) * length(PRIMER_LEVELS) + primer_order,
      x_label = paste0(as.character(tool), "\n", as.character(primer))
    )
}

make_tool_blocks <- function() {
  tibble(
    tool = factor(TOOL_LEVELS, levels = TOOL_LEVELS),
    tool_order = seq_along(TOOL_LEVELS),
    xmin = (tool_order - 1) * length(PRIMER_LEVELS) + 0.5,
    xmax = tool_order * length(PRIMER_LEVELS) + 0.5,
    ymin = -Inf,
    ymax = Inf
  )
}

classify_one_query_5class <- function(
    assignment_rank,
    n_hits_after_99_95,
    true_species,
    true_genus,
    true_family,
    top_hit_species_list,
    top_hit_genus_list,
    top_hit_family_list
) {

  if (!is.na(assignment_rank) && assignment_rank == "unassigned") {
    return("unassigned")
  }

  if (!is.na(n_hits_after_99_95) && as.numeric(n_hits_after_99_95) <= 0) {
    return("unassigned")
  }

  top_species <- split_taxon_list(top_hit_species_list)
  top_genus <- split_taxon_list(top_hit_genus_list)
  top_family <- split_taxon_list(top_hit_family_list)

  if (length(top_species) == 0 && length(top_genus) == 0 && length(top_family) == 0) {
    return("different_family_different_taxa")
  }

  if (all_equal_to_target(top_species, true_species)) {
    return("exact_species")
  }

  if (all_equal_to_target(top_genus, true_genus)) {
    return("same_genus_different_species")
  }

  if (all_equal_to_target(top_family, true_family)) {
    return("same_family_different_genus_species")
  }

  return("different_family_different_taxa")
}

make_depth_mark <- function(class_5) {
  dplyr::case_when(
    class_5 == "exact_species" ~ "",
    class_5 == "same_genus_different_species" ~ "x",
    class_5 == "same_family_different_genus_species" ~ "xx",
    class_5 == "different_family_different_taxa" ~ "xxx",
    TRUE ~ ""
  )
}

theme_bubble <- function() {
  theme_bw(base_size = FONT_SIZE_AXIS_X, base_family = FONT_FAMILY) +
    theme(
      text = element_text(
        family = FONT_FAMILY,
        color = TEXT_COLOR
      ),

      plot.title = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_TITLE,
        face = "bold",
        hjust = 0.5,
        color = TEXT_COLOR
      ),

      plot.subtitle = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_SUBTITLE,
        hjust = 0.5,
        color = TEXT_COLOR
      ),

      plot.caption = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_CAPTION,
        hjust = 0,
        color = TEXT_COLOR
      ),

      axis.text.x = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_AXIS_X,
        angle = X_LABEL_ANGLE,
        hjust = 0.5,
        vjust = 0.5,
        color = TEXT_COLOR
      ),

      axis.text.y = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_AXIS_Y,
        face = "italic",
        color = TEXT_COLOR
      ),

      axis.title.x = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_AXIS_TITLE_X,
        face = "bold",
        color = TEXT_COLOR
      ),

      axis.title.y = element_text(
        family = FONT_FAMILY,
        size = FONT_SIZE_AXIS_TITLE_Y,
        face = "bold",
        color = TEXT_COLOR
      ),

      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

make_annotation_y_positions <- function(title_y, labels) {

  wrapped_labels <- labels

  if (ANNOTATION_RESULT_WRAP) {
    wrapped_labels <- stringr::str_wrap(
      labels,
      width = ANNOTATION_RESULT_WRAP_WIDTH
    )
  }

  line_counts <- stringr::str_count(wrapped_labels, "\n") + 1

  y <- numeric(length(labels))
  current_y <- title_y - LEGEND_TITLE_ITEM_GAP

  for (i in seq_along(labels)) {
    y[i] <- current_y
    current_y <- current_y -
      LEGEND_ITEM_STEP -
      (line_counts[i] - 1) * ANNOTATION_RESULT_EXTRA_LINE_STEP
  }

  tibble::tibble(
    label = wrapped_labels,
    line_count = line_counts,
    y = y,
    next_title_y = current_y - LEGEND_GROUP_GAP
  )
}

make_unified_legend_panel <- function() {

  x_title <- 0.03
  x_symbol <- 0.12
  x_text <- 0.29

  y_start <- 0.95

  next_group_title_y <- function(last_item_y) {
    last_item_y - LEGEND_GROUP_GAP
  }

  make_item_y <- function(title_y, n_items) {
    title_y - LEGEND_TITLE_ITEM_GAP - LEGEND_ITEM_STEP * seq(0, n_items - 1)
  }

  y_toolkit_title <- y_start
  y_toolkit_items <- make_item_y(y_toolkit_title, 3)

  y_primer_title <- next_group_title_y(min(y_toolkit_items))
  y_primer_items <- make_item_y(y_primer_title, 3)

  y_tophit_title <- next_group_title_y(min(y_primer_items))
  y_tophit_items <- make_item_y(y_tophit_title, 3)

  y_anno_title <- next_group_title_y(min(y_tophit_items))
  anno_y_info <- make_annotation_y_positions(
    title_y = y_anno_title,
    labels = ANNOTATION_RESULT_LABELS
  )

  toolkit_df <- tibble::tibble(
    label = TOOL_LEVELS,
    x_symbol = x_symbol,
    x_text = x_text,
    y = y_toolkit_items,
    tool = factor(TOOL_LEVELS, levels = TOOL_LEVELS)
  )

  primer_df <- tibble::tibble(
    label = PRIMER_LEVELS,
    x_symbol = x_symbol,
    x_text = x_text,
    y = y_primer_items,
    shape_value = unname(PRIMER_SHAPES[PRIMER_LEVELS])
  )

  tophit_df <- tibble::tibble(
    label = c("1", "10", "40"),
    x_symbol = x_symbol,
    x_text = x_text,
    y = y_tophit_items,
    size_value = c(
      LEGEND_TOPHIT_SIZE_1,
      LEGEND_TOPHIT_SIZE_10,
      LEGEND_TOPHIT_SIZE_40
    )
  )

  annotation_df <- tibble::tibble(
    symbol = ANNOTATION_RESULT_SYMBOLS,
    label = anno_y_info$label,
    x_symbol = ANNOTATION_RESULT_SYMBOL_X,
    x_text = ANNOTATION_RESULT_TEXT_X,
    y = anno_y_info$y
  )

  y_min <- min(annotation_df$y) -
    LEGEND_ITEM_STEP * 1.2 -
    max(0, max(anno_y_info$line_count) - 1) * ANNOTATION_RESULT_EXTRA_LINE_STEP

  y_max <- y_start + LEGEND_ITEM_STEP * 1.2

  ggplot() +

  annotate(
    "text",
    x = x_title,
    y = y_toolkit_title,
    label = "Toolkit",
    hjust = 0,
    vjust = 0.5,
    family = FONT_FAMILY,
    fontface = "bold",
    size = LEGEND_TITLE_SIZE,
    color = TEXT_COLOR
  ) +
    geom_point(
      data = toolkit_df,
      aes(x = x_symbol, y = y, fill = tool),
      shape = 21,
      size = LEGEND_TOOLKIT_POINT_SIZE,
      stroke = 1.0,
      color = "grey20",
      show.legend = FALSE
    ) +
    geom_text(
      data = toolkit_df,
      aes(x = x_text, y = y, label = label),
      hjust = 0,
      vjust = 0.5,
      family = FONT_FAMILY,
      size = LEGEND_TEXT_SIZE,
      color = TEXT_COLOR,
      show.legend = FALSE
    ) +

  annotate(
    "text",
    x = x_title,
    y = y_primer_title,
    label = "Primer",
    hjust = 0,
    vjust = 0.5,
    family = FONT_FAMILY,
    fontface = "bold",
    size = LEGEND_TITLE_SIZE,
    color = TEXT_COLOR
  ) +
    geom_point(
      data = primer_df,
      aes(x = x_symbol, y = y, shape = shape_value),
      size = LEGEND_PRIMER_POINT_SIZE,
      fill = "grey70",
      color = "grey20",
      stroke = 1.0,
      show.legend = FALSE
    ) +
    geom_text(
      data = primer_df,
      aes(x = x_text, y = y, label = label),
      hjust = 0,
      vjust = 0.5,
      family = FONT_FAMILY,
      size = LEGEND_TEXT_SIZE,
      color = TEXT_COLOR,
      show.legend = FALSE
    ) +

  annotate(
    "text",
    x = x_title,
    y = y_tophit_title,
    label = "Top-hit count",
    hjust = 0,
    vjust = 0.5,
    family = FONT_FAMILY,
    fontface = "bold",
    size = LEGEND_TITLE_SIZE,
    color = TEXT_COLOR
  ) +
    geom_point(
      data = tophit_df,
      aes(x = x_symbol, y = y, size = size_value),
      shape = 16,
      color = "grey50",
      show.legend = FALSE
    ) +
    geom_text(
      data = tophit_df,
      aes(x = x_text, y = y, label = label),
      hjust = 0,
      vjust = 0.5,
      family = FONT_FAMILY,
      size = LEGEND_TEXT_SIZE,
      color = TEXT_COLOR,
      show.legend = FALSE
    ) +

  annotate(
    "text",
    x = x_title,
    y = y_anno_title,
    label = ANNOTATION_RESULT_TITLE,
    hjust = 0,
    vjust = 0.5,
    family = FONT_FAMILY,
    fontface = "bold",
    size = LEGEND_TITLE_SIZE,
    color = TEXT_COLOR
  ) +
    geom_text(
      data = annotation_df,
      aes(x = x_symbol, y = y, label = symbol),
      hjust = 0,
      vjust = 0.5,
      family = FONT_FAMILY,
      fontface = "bold",
      size = LEGEND_SYMBOL_TEXT_SIZE,
      color = TEXT_COLOR,
      show.legend = FALSE
    ) +
    geom_text(
      data = annotation_df,
      aes(x = x_text, y = y, label = label),
      hjust = 0,
      vjust = 0.5,
      family = FONT_FAMILY,
      size = LEGEND_TEXT_SIZE,
      lineheight = ANNOTATION_RESULT_LINEHEIGHT,
      color = TEXT_COLOR,
      show.legend = FALSE
    ) +

    scale_fill_manual(values = TOOL_COLORS, guide = "none") +
    scale_shape_identity(guide = "none") +
    scale_size_identity(guide = "none") +

    coord_cartesian(
      xlim = c(0, 1.40),
      ylim = c(y_min, y_max),
      clip = "off"
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      plot.margin = margin(42, 0, 30, 0),
      text = element_text(
        family = FONT_FAMILY,
        color = TEXT_COLOR
      )
    )
}

save_plot_all_formats <- function(plot_obj, file_prefix, width, height, dpi) {

  png_file <- paste0(file_prefix, ".png")
  tiff_file <- paste0(file_prefix, ".tiff")
  pdf_file <- paste0(file_prefix, ".pdf")
  svg_file <- paste0(file_prefix, ".svg")

  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "png",
    bg = "white"
  )

  ggsave(
    filename = tiff_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )

  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = svg_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = grDevices::svg,
    bg = "white"
  )

  cat("Saved raster/vector files:\n")
  cat("  ", png_file, "\n")
  cat("  ", tiff_file, "\n")
  cat("  ", pdf_file, "\n")
  cat("  ", svg_file, "\n")
}

save_plot_editable_office <- function(plot_obj, file_prefix, width, height, title_text) {

  pptx_file <- paste0(file_prefix, "_editable.pptx")
  docx_file <- paste0(file_prefix, "_word.docx")
  png_file <- paste0(file_prefix, ".png")

  ppt <- officer::read_pptx()
  ppt <- officer::add_slide(
    ppt,
    layout = "Blank",
    master = "Office Theme"
  )

  slide_size <- officer::slide_size(ppt)
  margin <- 0.25

  scale_factor <- min(
    (slide_size$width - 2 * margin) / width,
    (slide_size$height - 2 * margin) / height
  )

  fig_w <- width * scale_factor
  fig_h <- height * scale_factor
  fig_left <- (slide_size$width - fig_w) / 2
  fig_top <- (slide_size$height - fig_h) / 2

  ppt <- officer::ph_with(
    ppt,
    rvg::dml(ggobj = plot_obj),
    location = officer::ph_location(
      left = fig_left,
      top = fig_top,
      width = fig_w,
      height = fig_h
    )
  )

  print(ppt, target = pptx_file)

  doc <- officer::read_docx()

  landscape_sec <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape"),
    page_margins = officer::page_mar(
      top = 0.35,
      bottom = 0.35,
      left = 0.35,
      right = 0.35
    )
  )

  doc <- officer::body_set_default_section(doc, value = landscape_sec)

  doc <- officer::body_add_par(
    doc,
    title_text,
    style = "heading 1"
  )

  doc_max_w <- 10.6
  doc_max_h <- 7.1

  doc_scale <- min(doc_max_w / width, doc_max_h / height)
  doc_fig_w <- width * doc_scale
  doc_fig_h <- height * doc_scale

  doc <- officer::body_add_img(
    x = doc,
    src = png_file,
    width = doc_fig_w,
    height = doc_fig_h
  )

  doc <- officer::body_add_par(
    doc,
    "Figure note. Bubble indicates correct species-level assignment. Blank cells indicate unassigned queries. x, xx and xxx indicate same-genus mismatch, same-family mismatch and beyond-family mismatch, respectively. Bubble colour indicates toolkit, bubble shape indicates primer, and bubble size indicates top-hit count.",
    style = "Normal"
  )

  print(doc, target = docx_file)

  cat("Saved Office files:\n")
  cat("  Editable PPTX: ", pptx_file, "\n")
  cat("  Word file:     ", docx_file, "\n")
  cat("  Word source:   ", png_file, "\n")
}

if (!file.exists(INPUT_BUBBLE_DATA)) {
  stop(
    "INPUT_BUBBLE_DATA does not exist:\n",
    INPUT_BUBBLE_DATA,
    "\nPlease edit INPUT_BUBBLE_DATA at the top of this script."
  )
}

bubble_raw <- readr::read_csv(INPUT_BUBBLE_DATA, show_col_types = FALSE)

required_cols <- c(
  "tool",
  "primer",
  "true_species",
  "true_genus",
  "true_family",
  "assignment_rank",
  "n_hits_after_99_95",
  "top_hit_species_list",
  "top_hit_genus_list",
  "top_hit_family_list",
  "top_hit_subject_list"
)

missing_cols <- setdiff(required_cols, names(bubble_raw))
if (length(missing_cols) > 0) {
  stop(
    "Input bubble data is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

unknown_tools <- setdiff(unique(bubble_raw$tool), TOOL_LEVELS)
missing_color_names <- setdiff(TOOL_LEVELS, names(TOOL_COLORS))
missing_block_color_names <- setdiff(TOOL_LEVELS, names(TOOL_BLOCK_COLORS))

cat("\nToolkit names in input data:\n")
print(sort(unique(bubble_raw$tool)))

cat("\nToolkit levels expected by script:\n")
print(TOOL_LEVELS)

if (length(unknown_tools) > 0) {
  stop(
    "These toolkit names are present in the input data but not in TOOL_LEVELS: ",
    paste(unknown_tools, collapse = ", "),
    "\nPlease make TOOL_LEVELS exactly match the CSV tool column."
  )
}

if (length(missing_color_names) > 0) {
  stop(
    "These TOOL_LEVELS are missing from TOOL_COLORS: ",
    paste(missing_color_names, collapse = ", ")
  )
}

if (length(missing_block_color_names) > 0) {
  stop(
    "These TOOL_LEVELS are missing from TOOL_BLOCK_COLORS: ",
    paste(missing_block_color_names, collapse = ", ")
  )
}

plot_data <- bubble_raw %>%
  mutate(
    true_species = clean_taxon(true_species),
    true_genus = clean_taxon(true_genus),
    true_family = clean_taxon(true_family),
    n_hits_after_99_95 = as.numeric(n_hits_after_99_95)
  ) %>%
  rowwise() %>%
  mutate(
    class_5 = classify_one_query_5class(
      assignment_rank = assignment_rank,
      n_hits_after_99_95 = n_hits_after_99_95,
      true_species = true_species,
      true_genus = true_genus,
      true_family = true_family,
      top_hit_species_list = top_hit_species_list,
      top_hit_genus_list = top_hit_genus_list,
      top_hit_family_list = top_hit_family_list
    )
  ) %>%
  ungroup() %>%
  make_tool_primer_position() %>%
  mutate(
    is_annotated = class_5 != "unassigned",
    is_exact_species = class_5 == "exact_species",
    top_hit_count = count_semicolon_items(top_hit_subject_list),
    top_hit_count = ifelse(is_annotated, pmax(top_hit_count, 1L), 0L),
    top_hit_count_capped = pmin(top_hit_count, TOP_HIT_COUNT_CAP),
    x_mark = make_depth_mark(class_5),
    class_5_label = case_when(
      class_5 == "exact_species" ~ CLASS_LABEL_EXACT,
      class_5 == "unassigned" ~ CLASS_LABEL_UNASSIGNED,
      class_5 == "same_genus_different_species" ~ CLASS_LABEL_SAME_GENUS,
      class_5 == "same_family_different_genus_species" ~ CLASS_LABEL_SAME_FAMILY,
      class_5 == "different_family_different_taxa" ~ CLASS_LABEL_BEYOND_FAMILY,
      TRUE ~ "Other"
    ),
    class_5_label = factor(
      class_5_label,
      levels = CLASS_LABEL_LEVELS
    )
  )

if ("species_order" %in% names(plot_data)) {
  species_levels <- plot_data %>%
    distinct(true_species, species_order) %>%
    arrange(as.numeric(species_order)) %>%
    pull(true_species)
} else {
  species_levels <- unique(plot_data$true_species)
}

if (!KEEP_ALL_SPECIES_ON_Y_AXIS) {
  species_levels <- plot_data %>%
    filter(is_annotated) %>%
    pull(true_species) %>%
    unique()
}

if (REVERSE_SPECIES_ORDER) {
  species_levels <- rev(species_levels)
}

plot_data <- plot_data %>%
  mutate(
    true_species = factor(true_species, levels = species_levels),
    class_5 = factor(
      class_5,
      levels = c(
        "exact_species",
        "same_genus_different_species",
        "same_family_different_genus_species",
        "different_family_different_taxa",
        "unassigned"
      )
    )
  )

point_data <- plot_data %>%
  filter(if (HIDE_UNASSIGNED_POINTS) is_annotated else TRUE)

x_axis_labels <- plot_data %>%
  distinct(x_pos, x_label) %>%
  arrange(x_pos)

tool_blocks <- make_tool_blocks()

summary_tbl <- plot_data %>%
  group_by(tool, primer, x_pos, x_label) %>%
  summarise(
    C = sum(class_5 == "exact_species", na.rm = TRUE),
    U = sum(class_5 == "unassigned", na.rm = TRUE),
    x = sum(class_5 == "same_genus_different_species", na.rm = TRUE),
    xx = sum(class_5 == "same_family_different_genus_species", na.rm = TRUE),
    xxx = sum(class_5 == "different_family_different_taxa", na.rm = TRUE),
    annotated = sum(class_5 != "unassigned", na.rm = TRUE),
    mean_top_hits = mean(top_hit_count[class_5 != "unassigned"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tool = factor(tool, levels = TOOL_LEVELS),
    primer = factor(primer, levels = PRIMER_LEVELS),
    total = C + U + x + xx + xxx
  )

readr::write_csv(plot_data, OUT_PLOT_DATA)
readr::write_csv(summary_tbl, OUT_SUMMARY_TABLE)

p_main_core <- ggplot() +
  geom_rect(
    data = tool_blocks,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = tool),
    alpha = 0.13,
    inherit.aes = FALSE
  ) +

  geom_hline(
    yintercept = seq_along(species_levels),
    color = "grey92",
    linewidth = 0.25
  ) +

  geom_point(
    data = point_data %>% filter(is_exact_species),
    aes(
      x = x_pos,
      y = true_species,
      size = top_hit_count_capped,
      shape = primer,
      fill = tool
    ),
    color = EXACT_BORDER_COLOR,
    stroke = EXACT_POINT_STROKE,
    alpha = 0.95,
    show.legend = FALSE
  ) +

  geom_text(
    data = point_data %>% filter(x_mark != ""),
    aes(
      x = x_pos,
      y = true_species,
      label = x_mark
    ),
    family = FONT_FAMILY,
    size = FONT_SIZE_X_MARK,
    fontface = "bold",
    color = TEXT_COLOR,
    nudge_x = X_MARK_NUDGE_X,
    nudge_y = X_MARK_NUDGE_Y,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +

  scale_x_continuous(
    breaks = x_axis_labels$x_pos,
    labels = x_axis_labels$x_label,
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  scale_y_discrete(
    expand = expansion(add = 0.6),
    drop = !KEEP_ALL_SPECIES_ON_Y_AXIS
  ) +

  scale_fill_manual(
    values = TOOL_COLORS,
    guide = "none"
  ) +

  scale_shape_manual(
    values = PRIMER_SHAPES,
    guide = "none"
  ) +

  scale_size_continuous(
    range = SIZE_RANGE,
    breaks = SIZE_BREAKS,
    limits = c(1, TOP_HIT_COUNT_CAP),
    labels = SIZE_LABELS,
    guide = "none"
  ) +

  labs(
    x = MAIN_X_TITLE,
    y = MAIN_Y_TITLE,
    title = MAIN_TITLE
  ) +

  theme_bubble() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_legend <- make_unified_legend_panel()

p_main <- cowplot::plot_grid(
  p_main_core,
  p_legend,
  ncol = 2,
  rel_widths = c(MAIN_REL_WIDTH, LEGEND_REL_WIDTH),
  align = "h",
  axis = "tb"
)

save_plot_all_formats(
  plot_obj = p_main,
  file_prefix = OUT_MAIN_PREFIX,
  width = FINAL_WIDTH,
  height = FINAL_HEIGHT,
  dpi = FINAL_DPI
)

save_plot_editable_office(
  plot_obj = p_main,
  file_prefix = OUT_MAIN_PREFIX,
  width = FINAL_WIDTH,
  height = FINAL_HEIGHT,
  title_text = MAIN_TITLE
)

cat("\n============================================================\n")
cat("Five-class bubble plot finished.\n")
cat("Input data:\n")
cat("  ", INPUT_BUBBLE_DATA, "\n")

cat("\nOutput figure formats:\n")
cat("  PNG:  high-resolution bitmap, 800 dpi\n")
cat("  TIFF: high-resolution bitmap, 800 dpi\n")
cat("  PDF:  vector graphic\n")
cat("  SVG:  vector graphic\n")
cat("  PPTX: editable Office DrawingML\n")
cat("  DOCX: Word file with high-resolution PNG\n")

cat("\nOutput figure prefix:\n")
cat("  Main plot:      ", OUT_MAIN_PREFIX, "\n")

cat("\nOutput data:\n")
cat("  Plot data:      ", OUT_PLOT_DATA, "\n")
cat("  Summary table:  ", OUT_SUMMARY_TABLE, "\n")
cat("============================================================\n\n")

print(summary_tbl)
