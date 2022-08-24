library("tidyverse")

#colname <- c("chrm", "pos_start", "pos_end", "read_count")
#bed_out <- read_tsv("test_data/bedtools.out", col_names = colname)
#pos_start <- min(bed_out$pos_start)
#pos_end <- max(bed_out$pos_end)
#percnt_pos <- ceiling(quantile(seq(pos_end - pos_start), probs = seq(0, 1, 0.01))) + pos_start

args <- commandArgs(TRUE)

coverage_file <- args[1]
output_prefix <- args[2]

stopifnot(file.exists(coverage_file))

col_name <- c("Chr", "locus", "depth")

coverage <-
  read_tsv(coverage_file, col_names = col_name) %>% as_tibble()
coverage_df <-
  coverage %>% filter(locus >= 125) %>% filter(locus <= 712)


ggplot(coverage_df, aes(x = locus, y = depth)) +
  geom_line(
    colour = "red",
    size = 1,
    shape = 20,
    alpha = 1 / 3
  ) +
  ggtitle("Exon 4 to 7 of ABL1 - coverage plot") +
  theme_bw() +
  geom_vline(xintercept = c(274, 359, 537, 714)) +
  geom_text(aes(x = 274, label = "Exon4", y = 1), text = element_text(size =
                                                                        11)) +
  geom_text(aes(x = 359, label = "Exon5", y = 1), text = element_text(size =
                                                                        11)) +
  geom_text(aes(x = 537, label = "Exon6", y = 1), text = element_text(size =
                                                                        11)) +
  geom_text(aes(x = 713, label = "Exon7", y = 1), text = element_text(size =
                                                                        11))


ggsave(paste0(output_prefix, ".pdf"))
ggsave(paste0(output_prefix, ".png"))