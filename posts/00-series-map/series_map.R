#!/usr/bin/env Rscript
# Post 0 figure: the clickable series map (R version; series_map.py is the Python twin).
#
# Reads posts/series.tsv and writes a Plotly HTML page. Each post is a node; the GWAS lane runs
# above the shared trunk and the expression lane below, and both meet in fine-mapping.
# Published posts are filled and open when clicked; upcoming posts are rings.
#
# Usage: Rscript posts/00-series-map/series_map.R posts/series.tsv series-map.html [YYYY-MM-DD]
# The optional date decides which posts count as published (default: today). Without pandoc
# the page is written with its JavaScript in a side folder (series-map_files/).
suppressPackageStartupMessages({ library(data.table); library(plotly); library(htmlwidgets) })
args <- commandArgs(trailingOnly = TRUE)
tsv <- args[1]; out <- args[2]
today <- if (length(args) > 2) args[3] else format(Sys.Date())
site <- "https://vivekruhela.github.io/blog/"

# colours: the validated three-slot categorical set, light and dark steps
lanes <- data.table(lane = c("shared", "gwas", "expression"),
                    label = c("Shared trunk", "Genotype lane: GWAS and ROH", "Expression lane"),
                    light = c("#2a78d6", "#eb6834", "#1baf7a"), dark = c("#3987e5", "#d95926", "#199e70"))
ink <- list(light = c("#0b0b0b", "#52514e", "#c9c8c2", "#fcfcfb"), dark = c("#f0efec", "#c3c2b7", "#4a4a47", "#1a1a19"))

p <- fread(tsv, colClasses = "character")
p[, `:=`(x = as.numeric(x), y = as.numeric(y))]
p[, published := publish <= today]
p[, url := ifelse(published, paste0(site, publish, "-", slug, "/"), "")]
pos <- p[, .(id, x, y)]
setkey(pos, id)
after <- strsplit(p$after, ";", fixed = TRUE)
# a label above a node would sit on an edge that leaves straight up; put those below
up <- unique(unlist(Map(function(a, x, y) a[pos[a, x] == x & y > pos[a, y]], after, p$x, p$y)))
p[, textpos := ifelse(id %in% up, "bottom center", "top center")]

edges <- rbindlist(Map(function(a, x, y) if (length(a)) data.table(x = c(rbind(pos[a, x], x, NA)), y = c(rbind(pos[a, y], y, NA))),
                       after, p$x, p$y))
fig <- plot_ly(height = 470) |>
  add_trace(data = edges, x = ~x, y = ~y, type = "scatter", mode = "lines", hoverinfo = "skip",
            line = list(width = 2, color = ink$light[3]), showlegend = FALSE, name = "edges")
for (k in seq_len(nrow(lanes))) {
  d <- p[lane == lanes$lane[k]]
  hover <- sprintf("<b>%s  %s</b><br>%s<br>%s%s%s", d$id, d$title, d$question,
                   ifelse(d$published, "Published ", "Coming "), format(as.Date(d$publish), "%d %b %Y"),
                   ifelse(d$published, "<br><i>Click to open</i>", ""))
  fig <- fig |> add_trace(
    x = d$x, y = d$y, type = "scatter", mode = "markers+text", name = lanes$label[k], text = d$id,
    textposition = d$textpos, textfont = list(size = 13, color = ink$light[1]),
    marker = list(size = 26, color = ifelse(d$published, lanes$light[k], ink$light[4]),
                  line = list(width = 3, color = lanes$light[k])),
    customdata = d$url, hovertext = hover, hoverinfo = "text")
}

heads <- p[, .(x0 = min(x), x1 = max(x), lane = lane[1]), by = part]
ann <- c(lapply(seq_len(nrow(heads)), function(i) list(
  x = (heads$x0[i] + heads$x1[i]) / 2, y = if (heads$lane[i] == "expression") -1.55 else 2.75,
  text = sprintf("<b>%s</b>", heads$part[i]), showarrow = FALSE, font = list(size = 12, color = ink$light[2]))),
  list(list(x = 0, y = -0.11, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top", showarrow = FALSE,
            text = "Filled: published (click to open) \u00b7 Ring: coming soon", font = list(size = 12, color = ink$light[2]))))

fig <- fig |>
  layout(margin = list(l = 10, r = 10, t = 10, b = 75),
         paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
         font = list(family = "system-ui, -apple-system, Segoe UI, sans-serif", color = ink$light[1]),
         xaxis = list(visible = FALSE, range = c(-0.6, 13.6), fixedrange = TRUE),
         yaxis = list(visible = FALSE, range = c(-1.8, 3.0), fixedrange = TRUE),
         legend = list(orientation = "h", x = 0, y = -0.02, yanchor = "top", font = list(size = 12),
                       itemclick = FALSE, itemdoubleclick = FALSE),
         hoverlabel = list(align = "left"), dragmode = FALSE, annotations = ann) |>
  config(displayModeBar = FALSE, responsive = TRUE)

# click opens a published post in the page around the frame; dark mode follows the site's
# theme toggle when the frame is on the same site, else the OS setting (same as the Python page)
js <- sprintf("function (gd) {
  gd.on('plotly_click', function (e) { var u = e.points[0].customdata; if (u) { (window.top || window).location.href = u; } });
  function isDark() {
    try { var t = window.parent.document.documentElement.dataset.theme; if (t) return t === 'dark'; } catch (err) {}
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
  }
  function paint() {
    var d = isDark(), ink = d ? %s : %s, cols = d ? %s : %s;
    document.body.style.background = 'transparent';
    document.documentElement.style.colorScheme = d ? 'dark' : 'light';
    for (var k = 1; k < gd.data.length; k++) {
      var c = cols[k - 1], cd = [].concat(gd.data[k].customdata);
      Plotly.restyle(gd, {'marker.color': [cd.map(function (u) { return u ? c : ink[3]; })],
                          'marker.line.color': c, 'textfont.color': ink[0]}, [k]);
    }
    Plotly.restyle(gd, {'line.color': ink[2]}, [0]);
    var ann = (gd.layout.annotations || []).map(function (a) { a.font.color = ink[1]; return a; });
    Plotly.relayout(gd, {'font.color': ink[0], 'legend.font.color': ink[0], annotations: ann});
  }
  paint();
  if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', paint);
  try { new MutationObserver(paint).observe(window.parent.document.documentElement, {attributes: true, attributeFilter: ['data-theme']}); } catch (err) {}
}", jsonlite::toJSON(ink$dark), jsonlite::toJSON(ink$light), jsonlite::toJSON(lanes$dark), jsonlite::toJSON(lanes$light))
fig <- onRender(fig, js)

saveWidget(fig, out, selfcontained = requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available(),
           title = "Statistical genetics series map", background = "transparent")
cat(sprintf("wrote %s: %d posts, %d published as of %s\n", out, nrow(p), sum(p$published), today))
