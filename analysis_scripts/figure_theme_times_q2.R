if (.Platform$OS.type == "windows") {
  windowsFonts(`Times New Roman` = windowsFont("Times New Roman"))
}

Q2_FONT <- "Times New Roman"
Q2_COL <- c(
  ink = "#242424",
  neutral = "#7A7A7A",
  pale = "#D9D9D9",
  blue = "#3F6C8E",
  blue_light = "#DCE6ED",
  teal = "#4F8A7A",
  teal_light = "#DDEBE6",
  amber = "#C58A3A",
  amber_light = "#F3E7D3",
  red = "#B5524E",
  red_light = "#F3DEDC"
)

q2_par <- function(mar = c(3.0, 3.2, 1.8, 0.8)) {
  par(
    family = Q2_FONT,
    fg = Q2_COL[["ink"]],
    col.axis = Q2_COL[["ink"]],
    col.lab = Q2_COL[["ink"]],
    mar = mar,
    mgp = c(1.75, 0.42, 0),
    tcl = -0.20,
    cex.axis = 0.78,
    cex.lab = 0.82,
    xaxs = "i",
    yaxs = "i"
  )
}

q2_panel <- function(label, title, label_adj = 0, title_adj = 0.105) {
  mtext(label, side = 3, line = 0.16, adj = label_adj, font = 2, cex = 0.88, family = Q2_FONT)
  mtext(title, side = 3, line = 0.18, adj = title_adj, font = 2, cex = 0.76, family = Q2_FONT)
}

q2_arrow <- function(x0, y0, x1, y1, col = Q2_COL[["neutral"]], lwd = 1.2) {
  arrows(x0, y0, x1, y1, length = 0.055, angle = 22, code = 2, col = col, lwd = lwd, xpd = NA)
}

q2_card <- function(x0, y0, x1, y1, fill = "white", border = Q2_COL[["neutral"]], lwd = 1) {
  rect(x0, y0, x1, y1, col = fill, border = border, lwd = lwd, xpd = NA)
}

q2_save_base <- function(draw, base, width_mm = 180, height_mm = 125, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(base, ".svg"), width = w, height = h,
                     system_fonts = list(sans = Q2_FONT), pointsize = 8)
  } else {
    svg(paste0(base, ".svg"), width = w, height = h, family = Q2_FONT, pointsize = 8)
  }
  draw(); dev.off()
  cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = Q2_FONT, pointsize = 8)
  draw(); dev.off()
  png(paste0(base, ".png"), width = width_mm, height = height_mm, units = "mm",
      res = dpi, family = Q2_FONT, pointsize = 8, bg = "white")
  draw(); dev.off()
  tiff(paste0(base, ".tiff"), width = width_mm, height = height_mm, units = "mm",
       res = dpi, compression = "lzw", family = Q2_FONT, pointsize = 8, bg = "white")
  draw(); dev.off()
}
