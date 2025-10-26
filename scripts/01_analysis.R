# requirements ====
library(data.table)
library(ggplot2)
library(plotly)


# file paths ====
fp_dcm <- file.path(Sys.getenv("DCM_PGS_DIR"), "data", "anon_data.csv")
fp_gen <- NULL
N <- 300000
fp_gen <- "/mnt/project/dcm/score/score_ukb_pgs.sscore"
fp_pheno <- "/mnt/project/hermes3_data/hermes3_phenotypes.tsv.gz"


# read ====
dat <- rbindlist(list(
  dcm = fread(fp_dcm), 
  no_dcm = if (is.null(fp_gen)) {
    data.table(
      earlier_age = sample(0:100, size=N, replace=TRUE),
      tier_status = "no_dcm",
      Zheng_PRS   = rnorm(N, mean = -0.5, sd = 1),
      reported_sex= sample(c("F","M"), size=N, replace=TRUE)
    )[, genetic_sex := reported_sex]
  } else {
    d <- fread(fp_gen)
    pheno <- fread(fp_pheno)
  }
), idcol = "dcm_status")


# clean ====
setnames(dat, c("earlier_age","tier_status","Zheng_PRS"), c("age","gene_status","prs"))
dat[, age  := suppressWarnings(as.numeric(age))]
dat[, gene_status  := factor(gene_status,  
                             levels = c("no_dcm","gene_neg","gene_vus","gene_pos"),
                             labels = c("No DCM","Gene -ve","Gene VUS","Gene +ve"))]
dat[, reported_sex := factor(reported_sex, levels = c("M","F"), c("male","female"))]
dat[, genetic_sex  := factor(genetic_sex,  levels = c("M","F"), c("male","female"))]


# descriptives ====
summary(dat)


# ==== compute medians ====
meds <- dat[, .(median_prs = median(prs, na.rm = TRUE)), by = gene_status]


# plot PRS distributions ====
ggplot(dat, aes(x = prs, fill = gene_status)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = meds, aes(xintercept = median_prs, color = gene_status),
             linetype = "dashed", linewidth = 1) +
  labs(x = "Polygenic Risk Score (PRS)",
       y = "Density",
       fill = "Gene status",
       color = "Median PRS") +
  theme_minimal(base_size = 14)



# plot pie chart ====
dat_pie <- dat[gene_status != "No DCM", .N, by = gene_status][, prop := N / sum(N)]
ggplot(dat_pie, aes(x = "", y = prop, fill = gene_status)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  scale_fill_viridis_d() +
  labs(fill = "Gene status") +
  theme_void(base_size = 22)


# define thresholds ====
prs_thresholds <- data.table(
  or_thresh  = c(1.5, 5.0), 
  thresh_lab = c("intermediate", "pathogenic") 
)


# calculate ORs by PGS quantile ====
quantiles   <- c(0, 1, 5, 10, 20, 40, 60, 80, 90, 95, 99, 100)
dat[, prs_quantile := cut(prs, 
                          breaks = quantile(prs, probs = quantiles/100, na.rm = TRUE), 
                          labels = c(paste0("<",quantiles[2]), 
                                     paste0(quantiles[2:(length(quantiles)-2)], "-", quantiles[3:(length(quantiles)-1)]), 
                                     paste0(">", quantiles[length(quantiles)-1])), 
                          include.lowest = TRUE)]
dat[, prs_percentile := cut(prs, 
                            breaks = quantile(prs, probs = seq(0,1,0.01), na.rm = TRUE), 
                            labels = seq(1,100,1), 
                            include.lowest = TRUE)]
or_table <- data.table(prs_quantile = unique(dat$prs_quantile))
or_table[, `:=`(
  dcm_cases_in_group    = dat[prs_quantile == .BY[[1]] & gene_status != "No DCM", .N],
  dcm_controls_in_group = dat[prs_quantile == .BY[[1]] & gene_status == "No DCM", .N],
  dcm_cases_rest        = dat[prs_quantile != .BY[[1]] & gene_status != "No DCM", .N],
  dcm_controls_rest     = dat[prs_quantile != .BY[[1]] & gene_status == "No DCM", .N]
), by = prs_quantile]
or_table[, `:=`(
  odds_in_group = dcm_cases_in_group / dcm_controls_in_group,
  odds_in_rest  = dcm_cases_rest / dcm_controls_rest
)]
or_table[, `:=`(odds_ratio    = odds_in_group / odds_in_rest, 
                se_log_or     = sqrt(1 / dcm_cases_in_group + 
                                       1 / dcm_controls_in_group + 
                                       1 / dcm_cases_rest + 
                                       1 / dcm_controls_rest))]
or_table[, `:=`(log_or   = log(odds_ratio))]
or_table[, `:=`(ci_lower = exp(log_or - 1.96 * se_log_or),
                ci_upper = exp(log_or + 1.96 * se_log_or))]

# model the odds ratio across PGS quantile and extract the thresholds of 5 and 1.5
or_table[, pgs_midpoint := as.numeric(sub("-.*", "", as.character(prs_quantile))) + (as.numeric(sub(".*-", "", as.character(prs_quantile))) - as.numeric(sub("-.*", "", as.character(prs_quantile)))) / 2]
or_table[as.character(prs_quantile) == "<1", pgs_midpoint := 0.5]
or_table[as.character(prs_quantile) == ">99", pgs_midpoint := 99.5]

spline_model <- lm(log(odds_ratio) ~ splines::ns(pgs_midpoint, df = 6), data = or_table)
pred_prs_values <- seq(0,100,0.01)
predicted_or <- data.table(prs_quantile_val = pred_prs_values, 
                           pred_or          = exp(predict(spline_model, 
                                                          newdata = data.frame(pgs_midpoint = pred_pgs_values))))

# predict the thresholds from the fit 
prs_thresholds[, prs_cutoff := sapply(or_thresh, function(x) predicted_or[which.min(abs(pred_or-x)), prs_quantile_val])]
prs_thresholds[, prs := sapply(prs_cutoff, function(x) dat[ceiling(x)==as.numeric(prs_percentile), mean(prs, na.rm=T)])]

dat <- dat[
  prs_thresholds, 
  on = .(prs > prs), 
  `:=`(thresh_lab = paste0(i.thresh_lab, " PGS"))
]
dat[is.na(thresh_lab), `:=`(thresh_lab = "low PGS")]



ggplot(or_table, aes(x = pgs_midpoint, y = odds_ratio, ymin = ci_lower, ymax = ci_upper)) +
  geom_hline(yintercept = 1.0, color="darkgray", linetype = "dotted") +
  geom_pointrange() +
  geom_line(data = predicted_or, 
            aes(x = prs_quantile_val, y = pred_or, group=1), color = "darkred", inherit.aes = FALSE) +
  geom_segment(data = prs_thresholds, aes(x=0, xend = prs_cutoff, y = or_thresh, yend = or_thresh, color = thresh_lab), inherit.aes = F, linetype = "dashed") +
  geom_segment(data = prs_thresholds, aes(y=0, yend = or_thresh, x =prs_cutoff, xend=prs_cutoff, color = thresh_lab), inherit.aes = F, linetype = "dashed") +
  geom_text(data = prs_thresholds, aes(x=0, y = or_thresh, label = paste0(thresh_lab, " OR=",or_thresh, " (", sprintf("%.0f%%", prs_cutoff), ")"),  color = thresh_lab), vjust=-1, hjust=0, inherit.aes = F, show.legend = FALSE) +
  theme_classic(base_size = 12) +
  labs(x = "PGS quantile", y = "OR for DCM") +
  guides(color = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


ggplot(dat, aes(x = prs, fill = gene_status)) +
  geom_density(alpha=0.3) +
  geom_vline(data = prs_thresholds, aes(xintercept = prs), color="red") +
  geom_text(data = prs_thresholds, aes(x=prs, y = 0, label = paste0(thresh_lab, " OR=",or_thresh)), size=6, color = "red", vjust=-0.5, hjust=0, angle=90, inherit.aes = F, show.legend = FALSE) +
  scale_fill_viridis_d() +
  labs(x = "PGS", fill = "DCM status") +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y  = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank())


sunburst_data <- rbind(
  dat[dcm_status == "dcm", .(gene_status = "DCM", thresh_lab = "", N = .N)],
  dat[dcm_status == "dcm", .(thresh_lab = "", N = .N), by = "gene_status"],
  dat[dcm_status == "dcm", .N, by = .(gene_status, thresh_lab)]
)[, `:=`(id     = trimws(paste(gene_status, thresh_lab)),
         parent = fcase(gene_status=="DCM", "", 
                        thresh_lab==""    , "DCM", 
                        thresh_lab!=""    , as.character(gene_status)))]
sunburst_data[, label := fcase(gene_status=="DCM", "DCM", 
                               thresh_lab==""    , as.character(gene_status),
                               thresh_lab!=""    , thresh_lab)]
sunburst_data <- sunburst_data[, .(id, parent, N, label)]

# Create the sunburst plot
plot_ly(
  sunburst_data,
  ids     = ~id,
  labels  = ~label,
  parents = ~parent,
  values  = ~N,
  type    = 'sunburst',
  text = ~N,
  branchvalues = 'total'
) |> layout(colorway = c("#345799", "#e6550d", "#9ecae1"))




# Likelihood ratios and attributable risk faction
# calculate the odds ratio per Sdev of the PGS
model             <- glm(dcm_status=="dcm" ~ prs, data = dat, family = "binomial")
summary_model     <- summary(model)
odds_ratio_sd     <- exp(coef(summary_model)["prs", "Estimate"] * sd(dat$prs))
odds_ratio_sd_lci <- exp((coef(summary_model)["prs", "Estimate"] - 1.96*coef(summary_model)["prs", "Std. Error"]) * sd(dat$prs))
odds_ratio_sd_uci <- exp((coef(summary_model)["prs", "Estimate"] + 1.96*coef(summary_model)["prs", "Std. Error"]) * sd(dat$prs))

attributable_risk_frac <- function(prev, or_per_sd, centile) {
  odds   <- prev / (1-prev)
  u.A    <- log(or_per_sd)
  u.U    <- 0
  u.diff <- u.A - u.U
  z.U    <- qnorm(centile)
  z.A    <- z.U - u.A
  pgs_centile_aff <- pnorm(z.A)
  L.U    <- dnorm(z.U)
  L.A    <- dnorm(z.A)
  LR     <- L.A / L.U
  posttest.odds <- odds * LR
  posttest.P    <- posttest.odds / (posttest.odds + 1)
  afe           <- 1 - (prev / posttest.P)
  list(LR = LR, afe = afe)
}


dcm_prev <- 0.004 
likelihood_ratio <- data.table(
  dcm.P             = dcm_prev, 
  prs_or_per_sd     = odds_ratio_sd,
  prs_or_per_sd_lci = odds_ratio_sd_lci,
  prs_or_per_sd_uci = odds_ratio_sd_uci,
  prs_centile_unaff = seq(0.01,0.99,0.01)
)
likelihood_ratio[, c("LR", "afe")         := attributable_risk_frac(dcm.P, prs_or_per_sd, prs_centile_unaff)]
likelihood_ratio[, c("LR_lci", "afe_lci") := attributable_risk_frac(dcm.P, prs_or_per_sd_lci, prs_centile_unaff)]
likelihood_ratio[, c("LR_uci", "afe_uci") := attributable_risk_frac(dcm.P, prs_or_per_sd_uci, prs_centile_unaff)]



# Create the PGS centile vs proportion of DCM cases plot 
sd_lines <- data.table(y = c(ecdf(dat$prs)(mean(dat$prs) - 2*sd(dat$prs)),
                             ecdf(dat$prs)(mean(dat$prs) - sd(dat$prs)),
                             ecdf(dat$prs)(mean(dat$prs)),
                             ecdf(dat$prs)(mean(dat$prs) + sd(dat$prs)),
                             ecdf(dat$prs)(mean(dat$prs) + 2*sd(dat$prs))), 
                       x = 1,
                       text  = c("-2SD", "-1SD", "Mean", "+1SD", "+2SD"))

ggplot((dat[, .(dcm_count = sum(dcm_status == "dcm"), total_count = .N), by = prs_percentile]
              [order(as.numeric(prs_percentile)), ]
              [, `:=`(cumulative_dcm_count   = cumsum(dcm_count))]
              [, `:=`(cumulative_proportion  = cumulative_dcm_count / dat[dcm_status == "dcm", .N])]
              [, `:=`(se                     = sqrt((cumulative_proportion * (1 - cumulative_proportion)) / dat[dcm_status == "dcm", .N]))]), 
             aes(y = as.numeric(prs_percentile)/100, x = cumulative_proportion, group=1)) +
  geom_hline(data = sd_lines, aes(yintercept=y), linetype = "dotted", linewidth=1) +
  geom_text(data = sd_lines, aes(y=y, x=x, label=text), vjust=-0.3, hjust=0, size=5) +
  geom_line() + 
  geom_ribbon(aes(xmin = cumulative_proportion - 1.96*se, xmax = cumulative_proportion + 1.96*se), alpha = 0.2) + 
  labs(y = "PGS Percentile", x = "Proportion of DCM Cases") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent, expand = c(0,0.1)) + 
  theme_minimal(base_size = 18) +
  theme(panel.grid.minor = element_blank())



# plot the likelihood ratio by PGS centile
ggplot(likelihood_ratio, aes(x = prs_centile_unaff, y = LR)) +
  geom_line(color = "darkgray") +
  geom_point(size=3) +
  scale_x_continuous(breaks = seq(0,1,0.1), labels = scales::percent) +
  labs(x = "PGS centile in the unaffected", y = "Likelihood Ratio") +
  theme_classic(base_size = 18)



# plot the attributable risk factor for any given PGS centile
ggplot(likelihood_ratio[afe>=0], aes(x = prs_centile_unaff, y = afe)) +
  geom_ribbon(aes(ymin=afe_lci, ymax=afe_uci), alpha = 0.3) +
  geom_point(size=3) +
  scale_y_continuous(breaks = seq(0,1,0.1), labels = scales::percent) +
  scale_x_continuous(breaks = seq(0,1,0.05), labels = scales::percent) +
  labs(x = "PGS centile", y = "Attributable fraction among DCM") +
  theme_classic(base_size = 18)

