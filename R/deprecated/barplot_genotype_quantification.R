# create data.frame for plotting

gt_quant <- data.frame(
  seq_method = c("ngs", "ngs", "cloning", "cloning"),
  mendelian = c("Y", "N", "Y", "N"),
  count_in = c(5, 6, 17, 3),
  total_pairs = c(11, 11, 20, 20)
)

# quickly test whether groups differ statistically
a <- as.logical(c(rep("T", 17), rep("F", 3)))
t.test(a)

b <- as.logical(c(rep("T", 5), rep("F", 6)))
t.test(b)

# plot groups and results
ggplot(data = gt_quant, aes(x    = seq_method, 
                            y    = count_in, 
                            fill = mendelian)) +
  geom_bar(stat="identity", position = position_dodge()) +
  geom_text(aes(label = count_in), vjust=1.6, color="white",
            position = position_dodge(0.9), size=4) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,20),
                     name = "Total count\n") +
  scale_x_discrete(name = "\nSequencing method",
                   label = c("Cloning", "MiSeq")) +
  scale_fill_discrete(name = "Mendelian\ninheritance", label = c("No", "Yes"),
                      guide = guide_legend(reverse=TRUE)) +
  annotate(geom = "text", x = 1, y = 18, label = "***", size = 8) +
  annotate(geom = "text", x = 2, y =  7, label = "*", size = 8) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 14),
    plot.margin = unit(c(0.5,0.5,0.5,0,5), "cm"),
    legend.position = "right"
  )

ggsave(filename = "graphics/gtypes_quantification.png",
       plot = last_plot(),
       dpi = 400)

