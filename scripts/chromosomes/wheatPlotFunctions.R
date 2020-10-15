#!/applications/R/R-3.5.0/bin/Rscript

# Function to plot chromosome profiles of coverage (type = "h")
chrPlotCov <- function(xplot, title, cenStart, cenEnd,
                       dat1, col1, Ylab1, min1, max1) {
  plot(xplot, dat1, col = col1, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = title, cex.main = 2.5)
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = col1)
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 0.75, col = "black")
  box(lwd = 2.0)
}

# Function to plot chromosome profiles of coverage (type = "h")
chrPlotFeatureFreq <- function(xplot, title, cenStart, cenEnd,
                               dat1, col1, Ylab1, min1, max1) {
  plot(xplot, dat1, col = col1, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = title, cex.main = 2.5)
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = col1)
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 0.75, col = "black")
  box(lwd = 2.0)
}

# Function to plot chromosome profiles of coverage (type = "h")
# overlaid with feature frequency (type = "l")
chrPlotCovFeatureFreq <- function(xplot, title, cenStart, cenEnd,
                                  dat1, col1, Ylab1, min1, max1,
                                  dat2, col2, Ylab2, min2, max2) {
  plot(xplot, dat1, col = col1, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)*": "*italic("r"[s])*" = "*.(round(cor(dat1,
                                                                    dat2,
                                                                    method = "spearman"),
                                                                digits = 2))),
       cex.main = 2.5)
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = col1)
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  par(new = T)
  plot(xplot, dat2, col = col2, type = "l", lwd = 2,
       ylim = c(min2,
                max2),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 2.3, adj = c(0.5, -1.8), labels = Ylab2, xpd = NA, srt = -90, col = col2)
  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 0.75, col = "black")
  box(lwd = 2.0)
}

## Function to plot chromosome profiles of coverage (type = "h")
## overlaid with cM/Mb (type = "l")
#chrPlotCov_cMMb <- function(title, cenStart, cenEnd,
#                            xplot1, dat1, col1, Ylab1, min1, max1,
#                            xplot2, dat2, col2, Ylab2, min2, max2,
#                            legendLoc, legendLabs) {
#  plot(xplot1, dat1, col = col1, type = "h", lwd = 0.5,
#       ylim = c(min1,
#                max1),
#       xlim = c(0, max(chrLens)),
#       xlab = "", ylab = "",
#       xaxt = "n", yaxt = "n",
#       main = title,
#       cex.main = 2.5)
#  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
#       labels = c("0", "200", "400", "600", "800"),
#       at = c(0, 2e+08, 4e+08, 6e+08, 8e+08))
#  mtext(side = 2, line = 2.5, cex = 1.5, text = Ylab1, col = "black")
#  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
#  par(new = T)
#  plot(xplot2, dat2, col = col2, type = "l", lwd = 2,
#       ylim = c(min2,
#                max2),
#       xlim = c(0, max(chrLens)),
#       xlab = "", ylab = "",
#       xaxt = "n", yaxt = "n")
#  legend(legendLoc,
#         legend = legendLabs,
#         col = c(col1),
#         text.col = c(col1),
#         text.font = c(1),
#         ncol = 1, cex = 2.0, lwd = 1.5, bty = "n")
#  p <- par('usr')
#  text(p[2], mean(p[3:4]), cex = 2.3, adj = c(0.5, -1.8), labels = Ylab2, xpd = NA, srt = -90, col = col2)
#  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
#  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 0.75, col = "black")
#  box(lwd = 2.0)
#}

# Function to plot profile for one chromosome of one dataset (type = "h")
# overlaid with feature frequency or cM/Mb (type = "l")
chrPartitionPlotCov1_feature <- function(chrx, title, cenStart, cenEnd, rug1, rug1Col, 
                                         xplot1, dat1A, col1A, Ylab1, min1, max1,
                                         legendLoc, legendLabs,
                                         xplot2, dat2, col2, Ylab2, min2, max2) {
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)), font.main = 1,
       cex.main = 2.5)
  axis(side = 1, cex.axis = 2.0, lwd.tick = 1.5,
       labels = as.character(seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*100),
       at = seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*1e+08)
  mtext(side = 2, line = 2.50, cex = 2.0, text = Ylab1, col = "black")
  mtext(side = 1, line = 3.25, cex = 2.0, text = "Coordinates (Mb)", col = "black")
  axis(side = 2, cex.axis = 2.0, lwd.tick = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c("white"),
         text.col = c(col1A),
         text.font = c(1, 1),
         ncol = 1, cex = 2.0, lwd = 1.5, bty = "n")
  par(new = T)
  plot(xplot2, dat2, col = col2, type = "l", lwd = 2,
       ylim = c(min2,
                max2),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -1.4), labels = Ylab2, xpd = NA, srt = -90, col = col2)
  axis(side = 4, cex.axis = 2.0, lwd.tick = 1.5)
  xblocks(x = xplot1,
          y = (xplot1 < chrPartitions[x,2] |
               xplot1 > chrPartitions[x,5]),
          col = "midnightblue", height = 6.0)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,2] &
               xplot1 < chrPartitions[x,3]) |
              (xplot1 > chrPartitions[x,4] &
               xplot1 <= chrPartitions[x,5]),
          col = "turquoise3", height = 6.0)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,3] &
               xplot1 <= chrPartitions[x,4]),
          col = "lemonchiffon", height = 6.0)
  rug(x = rug1, ticksize = 0.03, side = 3, lwd = 0.75, col = rug1Col)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  box(lwd = 2.0)
}
#chrPartitionPlotCov1_feature <- function(chrx, title, cenStart, cenEnd, rug1, rug1Col, 
#                                         xplot1, dat1A, col1A, Ylab1, min1, max1, legendLoc, legendLabs,
#                                         xplot2, dat2, col2, Ylab2, min2, max2) {
#  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 0.5,
#       ylim = c(min1,
#                max1),
#       xlim = c(0, max(chrLens[chrx])),
#       xlab = "", ylab = "",
#       xaxt = "n", yaxt = "n",
#       main = bquote(.(title)), font.main = 1,
#       cex.main = 2.5)
#  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
#       labels = as.character(seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*100),
#       at = seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*1e+08)
#  mtext(side = 2, line = 2.25, cex = 2.0, text = Ylab1, col = "black")
#  mtext(side = 1, line = 3.25, cex = 2.0, text = "Coordinates (Mb)", col = "black")
#  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
#  legend(legendLoc,
#         legend = legendLabs,
#         col = c("white"),
#         text.col = c(col1A),
#         text.font = c(1, 1),
#         ncol = 1, cex = 2.0, lwd = 1.5, bty = "n")
#  par(new = T)
#  plot(xplot2, dat2, col = col2, type = "l", lwd = 2,
#       ylim = c(min2,
#                max2),
#       xlim = c(0, max(chrLens[chrx])),
#       xlab = "", ylab = "",
#       xaxt = "n", yaxt = "n")
#  p <- par('usr')
#  text(p[2], mean(p[3:4]), cex = 2.3, adj = c(0.5, -1.8), labels = Ylab2, xpd = NA, srt = -90, col = col2)
#  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
#  xblocks(x = xplot1,
#          y = (xplot1 < chrPartitions[x,2] |
#               xplot1 > chrPartitions[x,5]),
#          col = "midnightblue", height = 6.0)
#  xblocks(x = xplot1,
#          y = (xplot1 >= chrPartitions[x,2] &
#               xplot1 < chrPartitions[x,3]) |
#              (xplot1 > chrPartitions[x,4] &
#               xplot1 <= chrPartitions[x,5]),
#          col = "turquoise3", height = 6.0)
#  xblocks(x = xplot1,
#          y = (xplot1 >= chrPartitions[x,3] &
#               xplot1 <= chrPartitions[x,4]),
#          col = "lemonchiffon", height = 6.0)
#  rug(x = rug1, ticksize = 0.03, side = 3, lwd = 0.75, col = rug1Col)
#  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
#  box(lwd = 2.0)
#}

# Function to plot profile for one chromosome of two datasets (type = "h")
# overlaid with feature frequency or cM/Mb (type = "l")
chrPartitionPlotCov2_feature <- function(chrx, title, cenStart, cenEnd, rug1, rug1Col, 
                                         xplot1, dat1A, col1A, Ylab1, min1, max1,
                                         dat1B, col1B,
                                         legendLoc, legendLabs,
                                         xplot2, dat2, col2, Ylab2, min2, max2) {
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)), font.main = 1,
       cex.main = 2.5)
  lines(xplot1, dat1B, col = col1B, type = "h", lwd = 0.5)
  axis(side = 1, cex.axis = 2.0, lwd.tick = 1.5,
       labels = as.character(seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*100),
       at = seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*1e+08)
  mtext(side = 2, line = 2.50, cex = 2.0, text = Ylab1, col = "black")
  mtext(side = 1, line = 3.25, cex = 2.0, text = "Coordinates (Mb)", col = "black")
  axis(side = 2, cex.axis = 2.0, lwd.tick = 1.5)
  par(new = T)
  plot(xplot2, dat2, col = col2, type = "l", lwd = 2,
       ylim = c(min2,
                max2),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -1.4), labels = Ylab2, xpd = NA, srt = -90, col = col2)
  axis(side = 4, cex.axis = 2.0, lwd.tick = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c("white"),
         text.col = c(col1A, col1B),
         text.font = c(1, 1),
         ncol = 1, cex = 2.0, lwd = 1.5, bty = "n")
  xblocks(x = xplot1,
          y = (xplot1 < chrPartitions[x,2] |
               xplot1 > chrPartitions[x,5]),
          col = "midnightblue", height = 6.0)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,2] &
               xplot1 < chrPartitions[x,3]) |
              (xplot1 > chrPartitions[x,4] &
               xplot1 <= chrPartitions[x,5]),
          col = "turquoise3", height = 6.0)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,3] &
               xplot1 <= chrPartitions[x,4]),
          col = "lemonchiffon", height = 6.0)
  rug(x = rug1, ticksize = 0.03, side = 3, lwd = 0.75, col = rug1Col)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  box(lwd = 2.0)
}

# Function to plot profile for one chromosome of two log2(ChIP/input) datasets (type = "h")
# overlaid with feature frequency for two sets of features (type = "l")
chrPartitionPlotCov2_feature2 <- function(chrx, title, cenStart, cenEnd, rug1, rug1Col, 
                                          xplot1, dat1A, col1A, Ylab1, min1A, max1A, min1B, max1B,
                                          dat1B, col1B,
                                          legendLoc, legendLabs,
                                          xplot2, dat2A, col2A, Ylab2, min2A, max2A, min2B, max2B,
                                          dat2B, col2B) {
  par(mgp = c(3, 1, 0))
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 1.0,
       ylim = c(min1A,
                max1A),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)), font.main = 1,
       cex.main = 2.5)
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col1A, col.axis = col1A, line = 0.2)
#  lines(xplot1, dat1B, col = col1B, type = "h", lwd = 0.5)
  par(new = T, mgp = c(3, 1, 0))
  plot(xplot1, dat1B, col = col1B, type = "h", lwd = 1.0,
       ylim = c(min1B,
                max1B),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col1B, col.axis = col1B, line = 3.0)
  mtext(side = 2, line = 5.5, cex = 2.0, text = Ylab1, col = "black")
  par(mgp = c(3, 1.25, 0))
  axis(side = 1, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0,
       labels = as.character(seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*100),
       at = seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*1e+08)
  mtext(side = 1, line = 3.5, cex = 2.0, text = "Coordinates (Mb)", col = "black")
  par(new = T, mgp = c(3, 1.5, 0))
  plot(xplot2, dat2A, col = col2A, type = "l", lwd = 2,
       ylim = c(min2A,
                max2A),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
#  lines(xplot2, dat2B, col = col2B, type = "l", lwd = 2)
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col2A, col.axis = col2A, line = 0.2)
  par(new = T, mgp = c(3, 1.5, 0))
  plot(xplot2, dat2B, col = col2B, type = "l", lwd = 2,
       ylim = c(min2B,
                max2B),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col2B, col.axis = col2B, line = 3.0)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -3.5), labels = Ylab2, xpd = NA, srt = -90, col = "black")
  legend(legendLoc,
         inset = c(0.35, 0.08),
         legend = legendLabs,
         col = c("white"),
         text.col = c(col1A, col1B, col2A, col2B),
         text.font = c(1, 1, 1, 1),
         ncol = 1, cex = 1.25, lwd = 1.5, bty = "n")
  xblocks(x = xplot1,
          y = (xplot1 < chrPartitions[x,2] |
               xplot1 > chrPartitions[x,5]),
          col = "midnightblue", height = 6.0)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,2] &
               xplot1 < chrPartitions[x,3]) |
              (xplot1 > chrPartitions[x,4] &
               xplot1 <= chrPartitions[x,5]),
          col = "turquoise3", height = 6.0)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,3] &
               xplot1 <= chrPartitions[x,4]),
          col = "lemonchiffon", height = 6.0)
  rug(x = rug1, ticksize = 0.03, side = 3, lwd = 0.75, col = rug1Col)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  box(lwd = 2.0)
}

# Function to plot profile for one chromosome of two datasets (e.g., MNase and mCG; type = "h")
# overlaid with feature frequency for two sets of features (type = "l")
chrPartitionPlotCovMeth_feature2 <- function(chrx, title, cenStart, cenEnd,
#                                             rug1, rug1Col, 
                                             xplot1, dat1A, col1A, Ylab1, min1A, max1A, min1B, max1B,
                                             dat1B, col1B,
                                             legendLoc, legendLabs,
                                             xplot2, dat2A, col2A, Ylab2, min2A, max2A, min2B, max2B,
                                             dat2B, col2B) {
  par(mgp = c(3, 1, 0))
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 1.0,
       ylim = c(min1A,
                max1A),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)), font.main = 1,
       cex.main = 2.5)
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col1A, col.axis = col1A, line = 0.2)
#  lines(xplot1, dat1B, col = col1B, type = "h", lwd = 0.5)
  par(new = T, mgp = c(3, 1, 0))
  plot(xplot1, dat1B, col = col1B, type = "h", lwd = 1.0,
       ylim = c(min1B,
                max1B),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col1B, col.axis = col1B, line = 3.0)
  mtext(side = 2, line = 5.5, cex = 2.0, text = Ylab1, col = "black")
  par(mgp = c(3, 1.25, 0))
  axis(side = 1, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0,
       labels = as.character(seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*100),
       at = seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*1e+08)
  mtext(side = 1, line = 3.5, cex = 2.0, text = "Coordinates (Mb)", col = "black")
  par(new = T, mgp = c(3, 1.5, 0))
  plot(xplot2, dat2A, col = col2A, type = "l", lwd = 2,
       ylim = c(min2A,
                max2A),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
#  lines(xplot2, dat2B, col = col2B, type = "l", lwd = 2)
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col2A, col.axis = col2A, line = 0.2)
  par(new = T, mgp = c(3, 1.5, 0))
  plot(xplot2, dat2B, col = col2B, type = "l", lwd = 2,
       ylim = c(min2B,
                max2B),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col2B, col.axis = col2B, line = 3.0)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -3.5), labels = Ylab2, xpd = NA, srt = -90, col = "black")
  legend(legendLoc,
         inset = c(0.00, 0.00),
         legend = legendLabs,
         col = c("white"),
         text.col = c(col1A, col1B, col2A, col2B),
         text.font = c(1, 1, 1, 1),
         ncol = 1, cex = 1.25, lwd = 1.5, bty = "n")
  xblocks(x = xplot1,
          y = (xplot1 < chrPartitions[x,2] |
               xplot1 > chrPartitions[x,5]),
          col = "midnightblue", ybottom = 213)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,2] &
               xplot1 < chrPartitions[x,3]) |
              (xplot1 > chrPartitions[x,4] &
               xplot1 <= chrPartitions[x,5]),
          col = "turquoise3", ybottom = 213)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,3] &
               xplot1 <= chrPartitions[x,4]),
          col = "lemonchiffon", ybottom = 213)
#  rug(x = rug1, ticksize = 0.03, side = 3, lwd = 0.75, col = rug1Col)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  box(lwd = 2.0)
}

# Function to plot profile for one chromosome of two log2(ChIP/input) datasets (type = "h")
# overlaid with feature frequency for two sets of features (type = "l")
chrPartitionPlotCov2_featureQuantiles <- function(chrx, title, cenStart, cenEnd, rug1, rug1Col, 
                                                  xplot1, dat1A, col1A, Ylab1, min1A, max1A, min1B, max1B,
                                                  dat1B, col1B,
                                                  legendLoc, legendLabs,
                                                  xplot2, dat2A, col2A, Ylab2, min2A, max2A, min2B, max2B,
                                                  dat2B, col2B) {
  par(mgp = c(3, 1, 0))
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 1.0,
       ylim = c(min1A,
                max1A),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)), font.main = 1,
       cex.main = 2.5)
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col1A, col.axis = col1A, line = 0.2)
#  lines(xplot1, dat1B, col = col1B, type = "h", lwd = 0.5)
  par(new = T, mgp = c(3, 1, 0))
  plot(xplot1, dat1B, col = col1B, type = "h", lwd = 1.0,
       ylim = c(min1B,
                max1B),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col1B, col.axis = col1B, line = 3.0)
  mtext(side = 2, line = 5.5, cex = 2.0, text = Ylab1, col = "black")
  par(mgp = c(3, 1.25, 0))
  axis(side = 1, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0,
       labels = as.character(seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*100),
       at = seq(0, round(max(chrLens[chrx])/1e+08), by = 2)*1e+08)
  mtext(side = 1, line = 3.5, cex = 2.0, text = "Coordinates (Mb)", col = "black")
  par(new = T, mgp = c(3, 1.5, 0))
  plot(xplot2, dat2A, col = col2A, type = "l", lwd = 2,
       ylim = c(min2A,
                max2A),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
#  lines(xplot2, dat2B, col = col2B, type = "l", lwd = 2)
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = col2A, col.axis = col2A, line = 0.2)
  par(new = T, mgp = c(3, 1.5, 0))
  plot(xplot2, dat2B[,3], col = col2B[1], type = "l", lwd = 2,
       ylim = c(min2B,
                max2B),
       xlim = c(0, max(chrLens[chrx])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 4, cex.axis = 2.0, lwd = 2.0, lwd.tick = 2.0, col = "black", col.axis = "black", line = 3.0)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.0, adj = c(0.5, -3.5), labels = Ylab2, xpd = NA, srt = -90, col = "black")
  lines(xplot2, dat2B[,4], col = col2B[2], type = "l", lwd = 2.0)
  lines(xplot2, dat2B[,5], col = col2B[3], type = "l", lwd = 2.0)
  lines(xplot2, dat2B[,6], col = col2B[4], type = "l", lwd = 2.0)
  legend(legendLoc,
         inset = c(0.15, 0.08),
         legend = legendLabs,
         col = c("white"),
         text.col = c(col2B, col1A, col1B, col2A),
         text.font = c(1, 1, 1, 1),
         ncol = 2, cex = 1.25, lwd = 1.5, bty = "n")
  xblocks(x = xplot1,
          y = (xplot1 < chrPartitions[x,2] |
               xplot1 > chrPartitions[x,5]),
          col = "midnightblue", height = 0.16)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,2] &
               xplot1 < chrPartitions[x,3]) |
              (xplot1 > chrPartitions[x,4] &
               xplot1 <= chrPartitions[x,5]),
          col = "turquoise3", height = 0.16)
  xblocks(x = xplot1,
          y = (xplot1 >= chrPartitions[x,3] &
               xplot1 <= chrPartitions[x,4]),
          col = "lemonchiffon", height = 0.16)
  rug(x = rug1, ticksize = 0.03, side = 3, lwd = 0.75, col = rug1Col)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  box(lwd = 2.0)
}

# Function to plot chromosome profiles of coverage of one dataset (type = "h")
# overlaid with feature frequency or cM/Mb (type = "l")
chrPlotCov1_feature <- function(title, cenStart, cenEnd, R1End, R3Start, rug1, rug2, rug3,
                                regionCol, rug1Col, rug2Col, rug3Col,
                                xplot1, dat1A, col1A, Ylab1, min1, max1,
                                xplot2, dat2, col2, Ylab2, min2, max2) {
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)*": "*italic("r"[s])*" = "*.(round(cor(dat1A,
                                                                    dat2,
                                                                    method = "spearman",
                                                                    use = "pairwise.complete.obs"),
                                                                digits = 2))),
       cex.main = 2.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       labels = c("0", "200", "400", "600", "800"),
       at = c(0, 2e+08, 4e+08, 6e+08, 8e+08))
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = col1A)
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  par(new = T)
  plot(xplot2, dat2, col = col2, type = "l", lwd = 2,
       ylim = c(min2,
                max2),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 2.3, adj = c(0.5, -1.8), labels = Ylab2, xpd = NA, srt = -90, col = col2)
  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  abline(v = c(R1End, R3Start), lty = 5, lwd = 1, col = regionCol)
  rug(x = rug1, ticksize = 0.03, side = 1, lwd = 0.75, col = rug1Col)
  rug(x = rug2, ticksize = 0.03, side = 3, lwd = 0.75, col = rug2Col)
  rug(x = rug3, ticksize = 0.03, side = 3, lwd = 0.75, col = rug3Col)
  box(lwd = 2.0)
}

# Function to plot chromosome profiles of coverage of one dataset (type = "h")
# overlaid with feature quantiles (type = "l")
chrPlotCov1_featureQuantiles <- function(title, cenStart, cenEnd, R1End, R3Start, rug1, rug2, rug3,
                                         regionCol, rug1Col, rug2Col, rug3Col,
                                         xplot1, dat1A, col1A, Ylab1, min1, max1,
                                         xplot2, dat2, col2, Ylab2, min2, max2,
                                         legendLoc, legendLabs) {
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)), font.main = 1,
       cex.main = 2)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       labels = c("0", "200", "400", "600", "800"),
       at = c(0, 2e+08, 4e+08, 6e+08, 8e+08))
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = col1A)
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  par(new = T)
  plot(xplot2, dat2[[1]]$filt_feature, col = col2[[1]], type = "l", lwd = 2,
       ylim = c(min2,
                max2),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot2, dat2[[2]]$filt_feature, col = col2[[2]], type = "l", lwd = 2)
  lines(xplot2, dat2[[3]]$filt_feature, col = col2[[3]], type = "l", lwd = 2)
  lines(xplot2, dat2[[4]]$filt_feature, col = col2[[4]], type = "l", lwd = 2)
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 2.3, adj = c(0.5, -1.8), labels = Ylab2, xpd = NA, srt = -90, col = "black")
  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  abline(v = c(R1End, R3Start), lty = 5, lwd = 1, col = regionCol)
  rug(x = rug1, ticksize = 0.03, side = 1, lwd = 0.75, col = rug1Col)
  rug(x = rug2, ticksize = 0.03, side = 3, lwd = 0.75, col = rug2Col)
  rug(x = rug3, ticksize = 0.03, side = 3, lwd = 0.75, col = rug3Col)
  legend(legendLoc,
         legend = legendLabs,
         col = c(col2),
         text.col = c(col2),
         text.font = c(1),
         ncol = 1, cex = 1.5, lwd = 2, bty = "n")
  box(lwd = 2.0)
}

# Function to plot chromosome profiles of coverage of two datasets (type = "h")
# overlaid with cM/Mb (type = "l")
chrPlotCov2_feature <- function(title, cenStart, cenEnd, R1End, R3Start, rug1, rug2, rug3,
                                regionCol, rug1Col, rug2Col, rug3Col,
                                xplot1, dat1A, col1A, Ylab1, min1, max1,
                                dat1B, col1B,
                                legendLoc, legendLabs,
                                xplot2, dat2, col2, Ylab2, min2, max2) {
  plot(xplot1, dat1A, col = col1A, type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = bquote(.(title)*": "*italic("r"[s])*" = "*.(round(cor(dat1A,
                                                                    dat1B,
                                                                    method = "spearman"),
                                                                digits = 2))),
       cex.main = 2.5)
  lines(xplot1, dat1B, col = col1B, type = "h", lwd = 0.5)
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5,
       labels = c("0", "200", "400", "600", "800"),
       at = c(0, 2e+08, 4e+08, 6e+08, 8e+08))
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = "black")
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = c("white"),
         text.col = c(col1A, col1B),
         text.font = c(1, 1),
         ncol = 1, cex = 2.0, lwd = 1.5, bty = "n")
  par(new = T)
  plot(xplot2, dat2, col = col2, type = "l", lwd = 2,
       ylim = c(min2,
                max2),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 2.3, adj = c(0.5, -1.8), labels = Ylab2, xpd = NA, srt = -90, col = col2)
  axis(side = 4, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 1, col = "black")
  abline(v = c(R1End, R3Start), lty = 5, lwd = 1, col = regionCol)
  rug(x = rug1, ticksize = 0.03, side = 1, lwd = 0.75, col = rug1Col)
  rug(x = rug2, ticksize = 0.03, side = 3, lwd = 0.75, col = rug2Col)
  rug(x = rug3, ticksize = 0.03, side = 3, lwd = 0.75, col = rug3Col)
  box(lwd = 2.0)
}

# Function to plot chromosome profiles of coverage of 6 datasets (type = "h")
chrPlotCov6 <- function(title, cenStart, cenEnd,
                        xplot1, Ylab1, min1, max1, cols1,
                        dat1A,
                        dat1B,
                        dat1C,
                        dat1D,
                        dat1E,
                        dat1F,
                        legendLoc, legendLabs) {
  plot(xplot1, dat1A, col = cols1[1], type = "h", lwd = 0.5,
       ylim = c(min1,
                max1),
       xlim = c(0, max(chrLens)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = title,
       cex.main = 2.5)
  lines(xplot1, dat1B, col = cols1[2], type = "h", lwd = 0.5)
  lines(xplot1, dat1C, col = cols1[3], type = "h", lwd = 0.5)
  lines(xplot1, dat1D, col = cols1[4], type = "h", lwd = 0.5)
  lines(xplot1, dat1E, col = cols1[5], type = "h", lwd = 0.5)
  lines(xplot1, dat1F, col = cols1[6], type = "h", lwd = 0.5)
  mtext(side = 2, line = 2.25, cex = 1.5, text = Ylab1, col = "black")
  axis(side = 2, cex.axis = 1.5, lwd.tick = 1.5)
  legend(legendLoc,
         legend = legendLabs,
         col = cols1,
         text.col = c(cols1),
         text.font = c(1),
         ncol = 1, cex = 2.0, lwd = 1.5, bty = "n")
  axis(side = 1, cex.axis = 1.5, lwd.tick = 1.5)
  abline(v = c(cenStart, cenEnd), lty = 5, lwd = 0.75, col = "black")
  box(lwd = 2.0)
}

# Function to plot mean coverage profile of one chromatin mark around feature and random loci
plotAvgCov <- function(xplot,
                       dat1, ranDat1, col1, Ylab1,
                       flankSize, binSize,
                       flankLabL, flankLabR,
                       featureStartLab, featureEndLab,
                       ranLocStartLab, ranLocEndLab) {
  # Feature loci
  plot(xplot, dat1,
       ylim = c(min(c(dat1, ranDat1), na.rm = T),
                max(c(dat1, ranDat1), na.rm = T)),
       type = "l", lwd = 3, col = col1, ann = F,
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5,  at = pretty(c(dat1, ranDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylab1, col = col1)
  axis(side = 1, labels = c("", "", "", ""), lwd.tick = 1.5,
       at = c(1,
              (flankSize/binSize)+1,
              length(dat1)-(flankSize/binSize),
              length(dat1)))
  mtext(side = 1, line = 1, cex = 0.7,
        text = c(flankLabL,
                 featureStartLab,
                 featureEndLab,
                 flankLabR),
        at = c(1,
               (flankSize/binSize)+1,
               length(dat1)-(flankSize/binSize),
               length(dat1)))
  abline(v = c((flankSize/binSize)+1,
               length(dat1)-(flankSize/binSize)), lty = 3, lwd = 2)
  box(lwd = 2.0)

  # Random loci
  plot(xplot, ranDat1,
       ylim = c(min(c(dat1, ranDat1), na.rm = T),
                max(c(dat1, ranDat1), na.rm = T)),
       type = "l", lwd = 3, col = col1, ann = F,
       xaxt = "n", yaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5,  at = pretty(c(dat1, ranDat1)))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylab1, col = col1)
  axis(side = 1, labels = c("", "", "", ""), lwd.tick = 1.5,
       at = c(1,
              (flankSize/binSize)+1,
              length(dat1)-(flankSize/binSize),
              length(dat1)))
  mtext(side = 1, line = 1, cex = 0.7,
        text = c(flankLabL,
                 ranLocStartLab,
                 ranLocEndLab,
                 flankLabR),
        at = c(1,
               (flankSize/binSize)+1,
               length(ranDat1)-(flankSize/binSize),
               length(ranDat1)))
  abline(v = c((flankSize/binSize)+1,
               length(ranDat1)-(flankSize/binSize)), lty = 3, lwd = 2)
  box(lwd = 2.0)
}

# Function to plot mean coverage profile of one chromatin mark around feature and random loci
plotAvgCovYlim <- function(xplot,
                           dat1, ranDat1, col1, Ylab1, Ylim1,
                           flankSize, binSize,
                           flankLabL, flankLabR,
                           featureStartLab, featureEndLab,
                           ranLocStartLab, ranLocEndLab) {
  # Feature loci
  plot(xplot, dat1,
       ylim = Ylim1,
       type = "l", lwd = 3, col = col1, ann = F,
       xaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylab1, col = col1)
  axis(side = 1, labels = c("", "", "", ""), lwd.tick = 1.5,
       at = c(1,
              (flankSize/binSize)+1,
              length(dat1)-(flankSize/binSize),
              length(dat1)))
  mtext(side = 1, line = 1, cex = 0.7,
        text = c(flankLabL,
                 featureStartLab,
                 featureEndLab,
                 flankLabR),
        at = c(1,
               (flankSize/binSize)+1,
               length(dat1)-(flankSize/binSize),
               length(dat1)))
  abline(v = c((flankSize/binSize)+1,
               length(dat1)-(flankSize/binSize)), lty = 3, lwd = 2)
  box(lwd = 2.0)

  # Random loci
  plot(xplot, ranDat1,
       ylim = Ylim1,
       type = "l", lwd = 3, col = col1, ann = F,
       xaxt = "n")
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylab1, col = col1)
  axis(side = 1, labels = c("", "", "", ""), lwd.tick = 1.5,
       at = c(1,
              (flankSize/binSize)+1,
              length(dat1)-(flankSize/binSize),
              length(dat1)))
  mtext(side = 1, line = 1, cex = 0.7,
        text = c(flankLabL,
                 ranLocStartLab,
                 ranLocEndLab,
                 flankLabR),
        at = c(1,
               (flankSize/binSize)+1,
               length(ranDat1)-(flankSize/binSize),
               length(ranDat1)))
  abline(v = c((flankSize/binSize)+1,
               length(ranDat1)-(flankSize/binSize)), lty = 3, lwd = 2)
  box(lwd = 2.0)
}

# Function to plot mean coverage profile of one chromatin mark around genes in dominant, suppressed and balanced triads
plotAvgCovDomSupBal <- function(xplot,
                                balanced,
                                dominant_dominant, 
                                dominant_nondominant, 
                                suppressed_suppressed, 
                                suppressed_nonsuppressed, 
                                Ylab, colours,
                                flankSize, binSize,
                                flankLabL, flankLabR,
                                featureStartLab, featureEndLab,
                                legendLoc, legendLabs) {
  # Feature loci
  plot(xplot, balanced,
       ylim = c(min(balanced, dominant_dominant, dominant_nondominant,
                    suppressed_suppressed, suppressed_nonsuppressed),
                max(balanced, dominant_dominant, dominant_nondominant,
                    suppressed_suppressed, suppressed_nonsuppressed)),
       type = "l", lwd = 3, col = colours[1], ann = F,
       xaxt = "n")
  lines(xplot, dominant_dominant, type = "l", lwd = 3, col = colours[2])
  lines(xplot, dominant_nondominant, type = "l", lwd = 3, col = colours[3])
  lines(xplot, suppressed_suppressed, type = "l", lwd = 3, col = colours[4])
  lines(xplot, suppressed_nonsuppressed, type = "l", lwd = 3, col = colours[5])
  axis(side = 2, cex.axis = 1, lwd.tick = 1.5)
  mtext(side = 2, line = 2, cex = 0.8, text = Ylab, col = "black")
  axis(side = 1, labels = c("", "", "", ""), lwd.tick = 1.5,
       at = c(1,
              (flankSize/binSize)+1,
              length(dominant_dominant)-(flankSize/binSize),
              length(dominant_dominant)))
  mtext(side = 1, line = 1, cex = 0.7,
        text = c(flankLabL,
                 featureStartLab,
                 featureEndLab,
                 flankLabR),
        at = c(1,
               (flankSize/binSize)+1,
               length(dominant_dominant)-(flankSize/binSize),
               length(dominant_dominant)))
  abline(v = c((flankSize/binSize)+1,
               length(dominant_dominant)-(flankSize/binSize)), lty = 3, lwd = 2)
  legend(legendLoc,
         legend = legendLabs,
         col = colours,
         text.col = colours,
         text.font = c(1),
         ncol = 1, cex = 0.65, lwd = 1.5, bty = "n")
  box(lwd = 2.0)
}
