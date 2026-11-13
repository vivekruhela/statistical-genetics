# Figure helpers for the posts (R track); lib/figures.py is the Python twin.
# Figures are drawn in the light theme; lib/figure_theme.js switches them to dark where the
# reader's site or system is dark. With pandoc the page is self-contained; without it the
# JavaScript goes in a side folder next to the page.
suppressPackageStartupMessages({ library(plotly); library(htmlwidgets) })

COLORS <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
            magenta = "#e87ba4", green = "#008300")
INK <- "#0b0b0b"; INK2 <- "#52514e"; LINE <- "#c9c8c2"; SURFACE <- "#fcfcfb"
SEQ <- list(c(0, SURFACE), c(0.15, "#cde2fb"), c(0.5, "#6da7ec"), c(1, "#104281"))
FONT <- list(family = "system-ui, -apple-system, Segoe UI, sans-serif", size = 12, color = INK)

base_layout <- function(fig, ...) {
  layout(fig, paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)", font = FONT,
         dragmode = FALSE, hoverlabel = list(align = "left"),
         margin = list(l = 60, r = 20, t = 40, b = 50), ...)
}

save_figure <- function(fig, path, title, min_width = 640) {
  js <- paste(readLines(file.path("lib", "figure_theme.js")), collapse = "\n")
  fig <- config(fig, displayModeBar = FALSE, responsive = TRUE)
  fig <- onRender(fig, sprintf("function (gd) {\n%s\nstatgenTheme(gd, %d);\n}", js, min_width))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, mustWork = FALSE)
  saveWidget(fig, path, title = title, background = "transparent",
             selfcontained = requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available())
  cat("wrote", path, "\n")
}
