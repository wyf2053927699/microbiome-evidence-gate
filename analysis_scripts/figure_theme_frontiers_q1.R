# Frontiers/Q1 final-size theme for base-R manuscript figures.
# This file is preloaded before the frozen figure scripts. It changes only
# typography and rendering defaults; it does not alter source data or statistics.

options(stringsAsFactors = FALSE)

.q1_cex_floor <- 0.82       # 6.6 pt absolute floor at final size
.q1_title_cex <- 1.00       # 8 pt panel title
.q1_panel_cex <- 1.20       # 9.6 pt panel label

.q1_cex <- function(x, floor = .q1_cex_floor) {
  if (is.null(x) || !length(x) || !is.finite(x[1])) return(floor)
  max(as.numeric(x[1]), floor)
}

text <- function(x, y = NULL, labels = seq_along(x), adj = NULL, pos = NULL,
                 offset = 0.5, vfont = NULL, cex = 1, col = NULL,
                 font = NULL, ...) {
  graphics::text(x = x, y = y, labels = labels, adj = adj, pos = pos,
                 offset = offset, vfont = vfont, cex = .q1_cex(cex),
                 col = col, font = font, ...)
}

axis <- function(side, at = NULL, labels = TRUE, tick = TRUE, line = NA,
                 pos = NA, outer = FALSE, font = NA, lty = "solid",
                 lwd = 1, lwd.ticks = lwd, col = NULL, col.ticks = NULL,
                 hadj = NA, padj = NA, gap.axis = NA, cex.axis = 1, ...) {
  graphics::axis(side = side, at = at, labels = labels, tick = tick, line = line,
                 pos = pos, outer = outer, font = font, lty = lty, lwd = lwd,
                 lwd.ticks = lwd.ticks, col = col, col.ticks = col.ticks,
                 hadj = hadj, padj = padj, gap.axis = gap.axis,
                 cex.axis = .q1_cex(cex.axis), ...)
}

legend <- function(..., cex = 1, pt.cex = cex) {
  graphics::legend(..., cex = .q1_cex(cex), pt.cex = .q1_cex(pt.cex))
}

title <- function(main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                  line = NA, outer = FALSE, ...) {
  dots <- list(...)
  if (!is.null(main)) dots$cex.main <- max(dots$cex.main %||% 1, .q1_title_cex)
  dots$cex.lab <- max(dots$cex.lab %||% 1, .q1_cex_floor)
  do.call(graphics::title, c(list(main = main, sub = sub, xlab = xlab,
                                  ylab = ylab, line = line, outer = outer), dots))
}

mtext <- function(text, side = 3, line = 0, outer = FALSE, at = NA,
                  adj = NA, padj = NA, cex = NA, col = NA, font = NA, ...) {
  is_panel <- length(text) == 1 && grepl("^[A-Ha-h]$", text)
  target <- if (is_panel) .q1_panel_cex else .q1_cex_floor
  graphics::mtext(text = text, side = side, line = line, outer = outer, at = at,
                  adj = adj, padj = padj, cex = .q1_cex(cex, target),
                  col = col, font = font, ...)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
