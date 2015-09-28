# FIGURE: all four results plots with a common legend

library(ggplot2)
library(gridExtra)

results.final = {}

for (file in list.files("/Users/Beth/Desktop/ddi-results", full.names=TRUE)) {
    split.name = strsplit(file, "/")
    size = as.numeric(strsplit(split.name[[1]][6], "\\.")[[1]][1])
    d = read.delim(file, head = F)
    aucs = d[,1]
    results = data.frame(size, mean(aucs > 0.7), mean(aucs > 0.8), mean(aucs > 0.9), median(aucs), min(aucs), max(aucs), "ebc3d")
    names(results) = c("S", "fracg07", "fracg08", "fracg09", "median", "min", "max", "method")
    results.final = rbind(results.final, results)
}

colors = c("#629df5", "#3bbcba", "#9cebe9", "#efb331", "#f47433")


p1 <- ggplot(results.final, aes(x = S, y = fracg07, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("ebc3d" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))

p2 <- ggplot(results.final, aes(x = S, y = fracg08, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.8") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("ebc3d" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))

p3 <- ggplot(results.final, aes(x = S, y = fracg09, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.9") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("ebc3d" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))


g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend <- g_legend(p1)

png("/Users/Beth/Desktop/figure-ddis.png", height = 10, width = 10, units="in", res=300)

plot_full <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), mylegend, nrow = 2))

dev.off()

