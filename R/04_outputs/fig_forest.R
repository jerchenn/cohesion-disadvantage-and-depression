# fig_forest.R -- JAMA-style forest plots: log HR axis, estimate and 95% CI column at left, drawn
# from aggregate estimate files only (runs off the Workbench). One panel, so text and bars align.
# Usage: Rscript R/04_outputs/fig_forest.R <fig_estimates.csv> <out_dir> [fig_items_data.csv]
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(TRUE); est <- read.csv(args[1], check.names = FALSE); out <- args[2]

forest <- function(rows, file, xlab, breaks, lim, w = 7.2, cols = c("No.", "Events"), title = NULL, save = TRUE, mhead = "Measure", cihead = "HR (95% CI)"){
  rows$y <- rev(seq_len(nrow(rows))); top <- max(rows$y) + 1
  X0 <- 6.4; X1 <- 10                                  # plot region in panel units
  px <- function(h) X0 + (log(h) - log(lim[1])) / (log(lim[2]) - log(lim[1])) * (X1 - X0)
  rows$ci <- ifelse(rows$head, "", sprintf("%.2f (%.2f-%.2f)", rows$hr, rows$lo, rows$hi))
  f <- function(v) ifelse(is.na(v), "", formatC(v, format = "d", big.mark = " "))
  pts <- rows[!rows$head, ]
  g <- ggplot() +
    annotate("text", x = 0, y = top, label = mhead, hjust = 0, fontface = "bold", size = 2.9) +
    annotate("text", x = 3.55, y = top, label = if (is.null(cols)) "" else cols[1], hjust = 1, fontface = "bold", size = 2.9) +
    annotate("text", x = 4.35, y = top, label = if (is.null(cols)) "" else cols[2], hjust = 1, fontface = "bold", size = 2.9) +
    annotate("text", x = 6.1, y = top, label = cihead, hjust = 1, fontface = "bold", size = 2.9) +
    geom_text(data = rows, aes(0, y, label = label, fontface = ifelse(head, "bold", "plain")), hjust = 0, size = 2.9) +
    geom_text(data = rows, aes(3.55, y, label = f(n)), hjust = 1, size = 2.9) +
    geom_text(data = rows, aes(4.35, y, label = f(events)), hjust = 1, size = 2.9) +
    geom_text(data = rows, aes(6.1, y, label = ci), hjust = 1, size = 2.9) +
    annotate("segment", x = px(1), xend = px(1), y = 0.5, yend = top - 0.5, linetype = 2, colour = "grey50", linewidth = 0.35) +
    geom_segment(data = pts, aes(x = px(lo), xend = px(hi), y = y, yend = y), linewidth = 0.45) +
    geom_point(data = pts, aes(px(hr), y), shape = 15, size = 2) +
    annotate("segment", x = X0, xend = X1, y = 0.4, yend = 0.4, linewidth = 0.35) +
    annotate("segment", x = px(breaks), xend = px(breaks), y = 0.4, yend = 0.25, linewidth = 0.35) +
    annotate("text", x = px(breaks), y = 0.05, label = formatC(breaks, format = "f", digits = max(nchar(sub("^[^.]*\\.?", "", format(breaks, scientific = FALSE, drop0trailing = TRUE))))), size = 2.6, vjust = 1) +
    annotate("text", x = (X0 + X1) / 2, y = -0.55, label = xlab, size = 2.8) +
    annotate("segment", x = 0, xend = X1, y = top - 0.5, yend = top - 0.5, linewidth = 0.3, colour = "grey60") +
    scale_x_continuous(limits = c(0, X1 + 0.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.9, top + 0.4), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(6, 8, 4, 6), plot.title = element_text(face = "bold", size = 9, margin = margin(0, 0, 4, 0)))
  if (!is.null(title)) g <- g + ggtitle(title)
  h <- 0.24 * (nrow(rows) + 3) + 0.3 + if (is.null(title)) 0 else 0.2
  if (!save) return(list(g = g, h = h))
  ggsave(file.path(out, paste0(file, ".pdf")), g, width = w, height = h, device = "pdf")
  ggsave(file.path(out, paste0(file, ".eps")), g, width = w, height = h, device = "eps")
  ggsave(file.path(out, paste0(file, ".png")), g, width = w, height = h, dpi = 300)
  invisible(g)
}
grouped <- function(d, gcol){ do.call(rbind, lapply(unique(d[[gcol]]), function(k)
  rbind(data.frame(label = k, hr = NA, lo = NA, hi = NA, n = NA, events = NA, head = TRUE),
        data.frame(label = paste0("   ", d$label[d[[gcol]] == k]), hr = d$hr[d[[gcol]] == k], lo = d$lo[d[[gcol]] == k],
                   hi = d$hi[d[[gcol]] == k], n = d$n[d[[gcol]] == k], events = d$events[d[[gcol]] == k], head = FALSE)))) }

## Figure 3
con <- est[est$panel == "concordance", ]
grp <- c("Cohesion (exposure)" = "Exposure", "Self-report: happiness" = "Self-reported affect",
         "Self-report: meaning" = "Self-reported affect", "Self-report: social connection" = "Self-reported affect",
         "Device: steps" = "Device-measured activity", "Device: MVPA" = "Device-measured activity",
         "Device: sleep efficiency" = "Device-measured sleep", "Device: sleep WASO" = "Device-measured sleep")
nice <- c("Cohesion (exposure)" = "Neighborhood social cohesion", "Self-report: happiness" = "Happiness",
          "Self-report: meaning" = "Meaning in life", "Self-report: social connection" = "Social connection",
          "Device: steps" = "Daily steps", "Device: MVPA" = "Moderate to vigorous activity",
          "Device: sleep efficiency" = "Sleep efficiency", "Device: sleep WASO" = "Wake after sleep onset")
con$g <- grp[con$label]; con$label <- nice[con$label]
forest(grouped(con, "g"), "Figure3_forest", "HR per 1-SD higher (healthier) level (log scale)",
       c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1), c(0.5, 1.12))
cat("Wrote Figure3_forest.{pdf,eps,png}\n")

## Figure 2: item-level associations (A) and their interaction with disadvantage (B)
if (length(args) >= 3) {
  suppressPackageStartupMessages(library(gridExtra))
  it <- read.csv(args[3]); it$n <- NA; it$events <- NA
  ord <- c("Willing to help", "Can be trusted", "Get along", "Share values", "Watch out for each other")
  lab <- c("Marginal" = "Item alone", "Joint (independent)" = "All 5 items jointly",
           "x area deprivation" = "Item x area deprivation", "x lower income" = "Item x lower income")
  mk <- function(pn){ d <- it[it$panel == pn, ]; d$item <- factor(d$item, ord); d <- d[order(d$item), ]
    d$label <- lab[d$series]; d$g <- as.character(d$item); grouped(d, "g") }
  A <- forest(mk("A"), NULL, "HR per 1-SD higher item score (log scale)", c(0.85, 0.9, 0.95, 1.0, 1.05), c(0.84, 1.06),
              cols = NULL, mhead = "Item", title = "A  Association with incident depression", save = FALSE)
  B <- forest(mk("B"), NULL, "Interaction HR (log scale)", c(0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1),
              c(0.98, 1.10), cols = NULL, mhead = "Item", cihead = "Interaction HR (95% CI)", title = "B  Interaction with socioeconomic disadvantage", save = FALSE)
  g2 <- arrangeGrob(A$g, B$g, ncol = 1, heights = c(A$h, B$h)); h2 <- A$h + B$h
  ggsave(file.path(out, "Figure2_forest.pdf"), g2, width = 7.2, height = h2, device = "pdf")
  ggsave(file.path(out, "Figure2_forest.eps"), g2, width = 7.2, height = h2, device = "eps")
  ggsave(file.path(out, "Figure2_forest.png"), g2, width = 7.2, height = h2, dpi = 300)
  cat("Wrote Figure2_forest.{pdf,eps,png}\n")
}
