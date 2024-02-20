ggpairs_methyl <- function(txt, pdf) {
    library("ggplot2")
    library("Metrics")
    library("dplyr")
    library("tidyr")
    library("GGally")
    library("RColorBrewer")
    rf <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    r <- rf(32)
    bin2d_fn <- function(data, mapping, ...) {
        p <- ggplot(data = data, mapping = mapping) +
            geom_bin2d(bins = 25) +
            scale_fill_gradientn(colors = r, trans = "log10") + theme_bw()
        p
    }
    corr_fn <- function(x, y) {m<-mae(x,y); mae <- format(m, digits=4); ct <- cor.test(x,y);r <- unname(ct$estimate); corr <- format(r, digits=4); rt<-paste(mae,"\n", "Corr:", corr, sep=" ")}

    df <- read.table(txt, col.names = c("chr", "start", "perc", "depth", "sample"))
    df %>% filter(depth >= 5) %>% select(-c(depth)) %>% pivot_wider(id_cols = c(chr, start), names_from = sample, values_from = perc, names_sort = TRUE) %>% drop_na() -> pivot_df
    df_fig <- pivot_df[-c(1, 2)]
    graph <- ggpairs(df_fig,upper = list(continuous = wrap(ggally_statistic,text_fn=corr_fn, title="MAE")), lower = list(continuous = bin2d_fn)) + ggtitle(txt)
    ggsave(pdf, plot = graph, width = 9, height = 9, units = "in")
}
args <- commandArgs(trailingOnly = TRUE)
txt <- args[1]
pdf <- args[2]
ggpairs_methyl(txt, pdf)
