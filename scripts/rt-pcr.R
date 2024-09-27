require(dplyr)
require(ggplot2)
require(readr)
require(tidyr)
require(data.table)
require(openxlsx)
require(RColorBrewer)
require(plotly)
require(htmlwidgets)
require(sicegar)



## change
Cq <-
  read_csv(
    "patrycja/płytka 1_admin_2024-07-10 16-12-37_CT048023 -  Quantification Cq Results.csv",
    na = "NaN"
  )
amp <-
  read_csv(
    "patrycja/płytka 1_admin_2024-07-10 16-12-37_CT048023 -  Quantification Amplification Results_SYBR.csv"
  )

description <- read_delim(
  "plate1-description.tsv",
  delim = "\t",
  escape_double = FALSE,
  col_types = cols(primer = col_integer()),
  trim_ws = TRUE
)


basename <- "plate1-bis" ## output files prefix

Cq <- Cq %>%
  select(Well, Cq) %>%
  mutate(Well = gsub("(^|[^0-9])0+", "\\1", Well, perl = TRUE))

amp <- amp %>%
  select(matches("^[A-Z]")) %>%
  pivot_longer(!Cycle, names_to = "Well") %>%
  left_join(., Cq, by = "Well") %>%
  #left_join(., end_point, by = "Well") %>%
  filter(!is.na(Cq))


data <- amp %>%
  arrange(Well, Cycle) %>%
  group_by(Well) %>%
  mutate("baseline" = mean(value[5:10]), "max" = max(value)) %>%  ## to change cycles taken as baseline change
  ungroup() %>%                                                            ## "baseline" = mean(value[10:15]) to "baseline" = mean(value[from:to])
  mutate("scaled_value" = (value - baseline) / (max))

data.in.range <- data %>%
  filter(scaled_value >= 0.2) %>%
  group_by(Well) %>%
  mutate("first_cycle" = min(Cycle)) %>%
  ungroup() %>%
  filter(Cycle <= first_cycle + 2) %>%
  group_by(Well) %>%
  add_count() %>%
  ungroup() %>%
  as.data.table()

## sigmoidal model fitting
sig <- data %>%
  select(Cycle, scaled_value, Well) %>%
  rename("time" = Cycle, "intensity" = scaled_value) %>%
  as.data.frame()

sig_results <- data.frame(matrix(ncol = 12, nrow = 0))
wells <- unique(sig$Well)

for (well in wells) {
  print(well)
  sig_data <- sig %>%
    filter(Well == well) %>%
    select(time, intensity) %>%
    as.data.frame()
  
  fitObj <- sicegar::fitAndCategorize(dataInput = sig_data)
  summary_fit <- fitObj$summaryVector
  n <-  length(summary_fit)
  summary_fit <- t(summary_fit)
  
  print(n)
  if (n == 12) {
    summary_fit[1] <- well
    sig_results <- rbind(sig_results, summary_fit)
  }
  rm(fitObj, sum_sig)
}

sig_results <- sig_results %>%
  rename("Well" = dataInputName) %>%
  mutate_at(
    vars(
      "maximum_y",
      "midPoint_x",
      "midPoint_y",
      "slope",
      "incrementTime",
      "startPoint_x",
      "startPoint_y",
      "reachMaximum_x",
      "reachMaximum_y"
    ),
    as.numeric
  ) %>%
  mutate("E" = exp((4 * slope) / maximum_y)) %>%
  select(-maximum_x) %>%
  mutate(Well = as.character(Well)) %>%
  left_join(., description, by = "Well") %>%
  arrange(primer, strain, medium) %>%
  as.data.frame()

## log growth
slopes <-
  data.in.range[, coef(lm(log(scaled_value, 2) ~ Cycle))[2], by = Well] %>% as.data.frame()
colnames(slopes) <- c("Well", "slope")

r2 <-
  data.in.range[, summary(lm(log(scaled_value, 2) ~ Cycle))$r.squared, by = Well] %>% as.data.frame()
colnames(r2) <- c("Well", "r2")

intercepts <-
  data.in.range[, coef(lm(log(scaled_value, 2) ~ Cycle))[1], by = Well] %>% as.data.frame()
colnames(intercepts) <- c("Well", "intercept")

n <- data.in.range %>%
  select(Cq, n, Well) %>%
  unique() %>%
  as.data.frame()

result <- slopes %>%
  left_join(., r2, by = "Well") %>%
  left_join(., n, by = "Well") %>%
  left_join(., intercepts, by = "Well") %>%
  mutate("E" = 2 ^ slope) %>%
  mutate("C30" = (log(0.3, 2) - intercept) / slope) %>%
  left_join(., description, by = "Well") %>%
  arrange(primer, strain, medium) %>%
  as.data.frame()

excel <-
  buildWorkbook(result,
                sheetName = "statistics",
                colNames = TRUE,
                rowNames = FALSE)

addWorksheet(excel,
             sheetName = "statistics-sigmoidal-model")

writeData(
  excel,
  "statistics-sigmoidal-model",
  sig_results,
  colNames = TRUE,
  rowNames = FALSE
)

saveWorkbook(excel, paste(basename, "-results.xlsx", sep = ""), overwrite = TRUE)

### raw and scaled plots
ncolors <- amp %>%
  select(Well) %>%
  unique() %>%
  nrow()

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plot <-
  data %>% ggplot(aes(x = Cycle, y = scaled_value, color = Well)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  theme_bw() +
  xlab("Cycle") + ylab("Fluorescence (SYBR)") +
  theme(legend.background = element_rect(
    linewidth = 0.5,
    linetype = "solid",
    colour = "grey"
  )) +
  theme(text = element_text(size = 18)) +
  theme(legend.title = element_blank()) + ggtitle(paste("Scaled data: ", basename, sep =
                                                          "")) +
  theme(plot.title = element_text(hjust = 0.5, 18)) +
  scale_color_manual(values = getPalette(ncolors))

plot
plot <- ggplotly(plot)
plot1 <- plotly_build(plot)
plot1$x$layout$legend$title[[1]] <- ""
saveWidget(partial_bundle(plot1),
           paste(basename, ".html", sep = ""))

plot_raw <-
  amp %>% ggplot(aes(x = Cycle, y = value, color = Well)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  theme_bw() +
  xlab("Cycle") + ylab("Fluorescence (SYBR)") +
  theme(legend.background = element_rect(
    linewidth = 0.5,
    linetype = "solid",
    colour = "grey"
  )) +
  theme(text = element_text(size = 18)) +
  theme(legend.title = element_blank()) + ggtitle(paste("Raw data: ", basename, sep =
                                                          "")) +
  theme(plot.title = element_text(hjust = 0.5, 18)) +
  scale_color_manual(values = getPalette(ncolors))
plot_raw
plot_raw <- ggplotly(plot_raw)
plot1_raw <- plotly_build(plot_raw)
plot1_raw$x$layout$legend$title[[1]] <- ""
saveWidget(partial_bundle(plot1_raw),
           paste(basename, "-raw.html", sep = ""))
