# vectorize distance matrices

msats_gen.abs[upper.tri(msats_gen.abs, diag = T)] <- NA
a <- msats_gen.abs %>% as.vector() %>% na.omit()

clone_gen.ufrac[upper.tri(clone_gen.ufrac)] <- NA
diag(clone_gen.ufrac) <- NA
b <- clone_gen.ufrac %>% as.vector() %>% na.omit()

df <- cbind(a,b) %>% as.data.frame()

model <- lm(b ~ a)
m<- summary(model)
m


(
  mantel_plot <- ggplot(df, aes(b,a)) +
  geom_point() +
    # geom_jitter(height = 0.002, width = 0.07) +
  geom_smooth(method = "lm", 
              se = T, 
              color = "orange") +
    scale_y_continuous(name   = "Pairwise microsatellite\nallele sharing",
                       limits = c(20, 70),
                       expand = c(0,0)) +
    scale_x_continuous(name   = "Pairwise MHC distances") +
    annotate(geom = "text", x = 0.03, y = 69,
             label ="italic(R²)==0.01", parse = T) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = unit(c(0.5,0.5,0.5,0,5), "cm"),
      axis.title = element_text(color = "black", 
                                margin = margin(10,10,20,10))
    )
)

ggsave(filename = "graphics/isolation_by_distance.png", dpi = 400)
