library(data.table)
library(dplyr)
library(R.utils)
library(ggplot2)
library(purrr)


setwd("~/Documents/SDZWA/alala_genomics/analysis/mosdepth/data")

files <- list.files(pattern = "\\.regions\\.bed\\.gz$")

cov <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  setnames(dt, c("chr", "start", "end", "coverage"))
  dt[, sample := sub("\\.regions\\.bed\\.gz$", "", f)]
  dt
}))

cov <- cov[chr == "NC_063222.1" & start >= 14e6 & end <= 24e6]



depth <- ggplot(cov, aes(x = (start + end)/2, y = coverage, group = sample)) +
  geom_line(alpha = 0.2) +
  labs(
    x = "Position (Mb)",
    y = "Coverage (×)"
  ) +
  theme_bw()

ggsave("~/Documents/SDZWA/alala_genomics/analysis/mosdepth/chr10_depth.png",
       depth, width = 8, height = 2)




walk(unique(cov$sample), function(s) {
  p <- ggplot(filter(cov, sample == s),
              aes(x = (start + end)/2, y = coverage)) +
    geom_line() +
    labs(
      title = s,
      x = "Genomic position (chr10)",
      y = "Coverage (×)"
    ) +
    theme_bw()
  
  ggsave(paste0("coverage_", s, ".png"), p, width = 6, height = 3, dpi = 300)
})

