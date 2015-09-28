# FIGURE: all four results plots with a common legend

library(ggplot2)
library(gridExtra)

results.final = read.delim("~/Desktop/results-final.tsv")

results.final = subset(results.final, method %in% c("ebcranksum"))

results.pgx = results.final[results.final$datatype == "gene.2d",]
results.target = results.final[results.final$datatype == "target.2d",]
results.pgx.full = results.final[results.final$datatype == "gene.3d",]
results.target.full = results.final[results.final$datatype == "target.3d",]

colors = c("#629df5", "#3bbcba", "#9cebe9", "#efb331", "#f47433")

png("/Users/Beth/Desktop/figure-comparison-2d-3d-pgx-target.png", height = 10, width = 10, units="in", res=300)

p1 <- ggplot(results.pgx, aes(x = S, y = fracg07, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + labs(title = "(a)") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("avgcosine" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))

p2 <- ggplot(results.target, aes(x = S, y = fracg07, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + labs(title = "(b)") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("avgcosine" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))

p3 <- ggplot(results.pgx.full, aes(x = S, y = fracg07, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + labs(title = "(c)") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("avgcosine" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))

p4 <- ggplot(results.target.full, aes(x = S, y = fracg07, colour = method)) + geom_line(size = 0.75) + geom_point(aes(colour = method), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + labs(title = "(d)") + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("avgcosine" = colors[1], "bestcosine" = colors[2], "ranksum" = colors[3], "ebctotal" = colors[4], "ebcranksum" = colors[5]))

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

names(results.pgx)[9] = "Method"

mylegend <- g_legend(p1)

plot_full <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), nrow = 2), mylegend, nrow = 2, heights = c(10, 1))

dev.off()

# FIGURE: comparing different iteration numbers

results.pgx = results.final[grep("gene", results.final$datatype),]

p1 <- ggplot(results.pgx, aes(x = S, y = fracg07, colour = datatype)) + geom_line(size = 0.75) + geom_point(aes(colour = datatype), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("gene.2d.iterations.0" = colors[1], "gene.2d.iterations.1" = colors[2], "gene.2d.iterations.2" = colors[3]), labels = c("iterations = 0", "iterations = 1", "iterations = 2"))

p2 <- ggplot(results.pgx, aes(x = S, y = fracg08, colour = datatype)) + geom_line(size = 0.75) + geom_point(aes(colour = datatype), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.8") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("gene.2d.iterations.0" = colors[1], "gene.2d.iterations.1" = colors[2], "gene.2d.iterations.2" = colors[3]), labels = c("iterations = 0", "iterations = 1", "iterations = 2"))

p3 <- ggplot(results.pgx, aes(x = S, y = fracg09, colour = datatype)) + geom_line(size = 0.75) + geom_point(aes(colour = datatype), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.9") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("gene.2d.iterations.0" = colors[1], "gene.2d.iterations.1" = colors[2], "gene.2d.iterations.2" = colors[3]), labels = c("iterations = 0", "iterations = 1", "iterations = 2"))

mylegend <- g_legend(p1)

plot_full <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), mylegend, nrow = 2))


results.target = results.final[grep("target", results.final$datatype),]

p1 <- ggplot(results.target, aes(x = S, y = fracg07, colour = datatype)) + geom_line(size = 0.75) + geom_point(aes(colour = datatype), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.7") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("target.2d.iterations.0" = colors[1], "target.2d.iterations.1" = colors[2], "target.2d.iterations.2" = colors[3]), labels = c("iterations = 0", "iterations = 1", "iterations = 2"))

p2 <- ggplot(results.target, aes(x = S, y = fracg08, colour = datatype)) + geom_line(size = 0.75) + geom_point(aes(colour = datatype), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.8") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("target.2d.iterations.0" = colors[1], "target.2d.iterations.1" = colors[2], "target.2d.iterations.2" = colors[3]), labels = c("iterations = 0", "iterations = 1", "iterations = 2"))

p3 <- ggplot(results.target, aes(x = S, y = fracg09, colour = datatype)) + geom_line(size = 0.75) + geom_point(aes(colour = datatype), size = 3) + xlab("Seed Set Size") + ylab("Fraction of Trials with AUC > 0.9") + xlim(0,100) + ylim(0,1) + theme(plot.title = element_text(hjust = 0)) + theme(plot.title = element_text(lineheight=.8, face="bold")) + theme(panel.background = element_rect(fill = "gray95")) + scale_colour_manual(values = c("target.2d.iterations.0" = colors[1], "target.2d.iterations.1" = colors[2], "target.2d.iterations.2" = colors[3]), labels = c("iterations = 0", "iterations = 1", "iterations = 2"))

mylegend <- g_legend(p1)

plot_full <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), mylegend, nrow = 2))

