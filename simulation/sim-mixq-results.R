library(tidyverse)
library(gt)

load("output/out-sim-mixq.RData")

out <- map2(cases, out, \(x, y) c(x, "out" = list(y))) |>
  list_transpose() |>
  as_tibble() |>
  unnest(cols = "out") |>
  pivot_longer(cols = starts_with("p"), names_prefix = "p_", names_to = "p", values_to = "rmse") |>
  mutate(method = factor(method,
                         levels = c("mixq", "qr", "qrcm", "gev", "egp", "evgam_gev", "evgam_gpd"),
                         labels = c("mixQ", "QR", "QRCM", "GEV-res", "EGP-res", "GEV-evgam", "GPD-evgam")))

# count the reps with NA for egp
out |>
  group_by(method) |>
  summarise(sum(is.na(rmse)))
out |>
  filter(method == "EGP-res", p == 0.95) |>
  group_by(dgp, n) |>
  summarise(sum(is.na(rmse)))


dgp <- c("pareto","logn", "egp")

p_median <- map(dgp, \(x) out |>
               filter(dgp == x,
                      method != "mixq2") |>
               group_by(n, p) |>
               mutate(rmse = rmse / median(rmse, na.rm = TRUE)) |>
               ggplot(aes(x = method, y = rmse))+
               geom_boxplot(fill = 2, outlier.shape = 1, width = 0.5, staplewidth  = 0.5)+
               facet_grid(rows = vars(p),
                          cols = vars(n),
                          labeller = \(x) label_both(labels = x, sep = " = "),
                          scales = "free")+
               labs(y = "RMSE / median(RMSE)", x = NULL)+
               coord_cartesian(ylim = c(0, 5))+
               theme_bw()+
               theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 9, hjust = 1),
                     strip.background = element_rect(fill = NA),
                     strip.text.y = element_text(angle = 0),
                     axis.text.x = element_text(angle = -30, hjust = 0)))

map2(p_median, dgp, \(x, y) ggsave(plot = x,
                                   filename = paste0("output/sim-mixq-rmse-median-",y,".pdf"),
                                   width = 10, height = 6))

# tables
gt_out <- purrr::map(dgp,
           \(x) out |>
             filter(dgp == x) |>
             group_by(n, p, method) |>
             summarise("rmse" = mean(rmse, na.rm = TRUE), .groups = "drop") |>
             pivot_wider(names_from = method, values_from = rmse) |>
             arrange(desc(p), n) |>
             gt() |>
             fmt_number(columns = -c(n, p), decimals = 2))

map2(dgp, gt_out, \(x, y) gtsave(data = y, filename = paste0("output/sim-mixq-rmse-tab-",x,".html")))
walk(gt_out, \(x) x |>
      as_latex() |>
      as.character() |>
      writeLines()
)

# not used
# p_diff <- map(dgp, \(x) out |>
#                 filter(dgp == x,
#                        method != "mixq2") |>
#                 pivot_wider(names_from = method, values_from = rmse) |>
#                 mutate(across(-c(dgp, n, p, rep), \(x) (x - mixQ))) |>
#                 select(-mixQ) |>
#                 pivot_longer(cols = -c(dgp, n, p, rep), values_to = "rmse", names_to = "method") |>
#                 mutate(method = factor(method, levels = c("QR", "QRCM", "GEV"))) |>
#                 ggplot(aes(x = method, y = rmse))+
#                 geom_hline(yintercept = 0, lty = 3)+
#                 geom_boxplot(fill = 2, outlier.shape = 1, width = 0.5, staplewidth  = 0.5)+
#                 facet_grid(rows = vars(p),
#                            cols = vars(n),
#                            labeller = \(x) label_both(labels = x, sep = " = "),
#                            scales = "free")+
#                 labs(y = expression(RMSE - RMSE["mixQ"]), x = NULL)+
#                 theme_bw()+
#                 theme(axis.title.y = element_text(angle = 0, vjust = 0.5, size = 9, hjust = 1),
#                       strip.background = element_rect(fill = NA),
#                       strip.text.y = element_text(angle = 0)))

# map2(p_diff, dgp, \(x, y) ggsave(plot = x,
#                                  filename = paste0("output/sim-mixq-rmse-diff-",y,".pdf"),
#                                  width = 7, height = 6))
