#======================#
#     HELPER FUNCTIONS
#======================#

library(marginaleffects)
library(broom.mixed)
library(dplyr)
library(tidyr)

# 1. Get absolute risk predictions
get_abs_risk <- function(model, outcome_label, group_var = "allocation_f", type = "response") {
  avg_predictions(model, variables = group_var, type = type) %>%
    as.data.frame() %>%
    rename(Group = !!sym(group_var), Risk = estimate, CI_Lower = conf.low, CI_Upper = conf.high) %>%
    mutate(
      Risk = round(Risk, 3),
      CI_Lower = round(CI_Lower, 3),
      CI_Upper = round(CI_Upper, 3),
      Risk_CI = sprintf("%.3f (%.3f, %.3f)", Risk, CI_Lower, CI_Upper),
      Outcome = outcome_label
    ) %>%
    select(Outcome, Group, Risk_CI)
}

# 2. Get risk difference
get_risk_difference <- function(model, outcome_label, group_var = "allocation_f", type = "response") {
  #avg_comparisons(model, variables = group_var, type = type) %>%
  comparisons(model, variables = group_var, type = type) %>%  
  as.data.frame() %>%
    slice(1) %>%
    mutate(
      estimate = round(estimate, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3),
      Group = paste0("RD (", contrast, ")"),
      Risk_CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
      Outcome = outcome_label
    ) %>%
    select(Outcome, Group, Risk_CI)
}

# 3. Get OR or IRR
get_or_irr <- function(model, outcome_label, term_name = "allocation_fIntervention", add_star = FALSE) {
  tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == term_name) %>%
    mutate(
      estimate = round(estimate, 2),
      conf.low = round(conf.low, 2),
      conf.high = round(conf.high, 2),
      OR_CI = sprintf("%.3f (%.3f, %.3f)%s", estimate, conf.low, conf.high, ifelse(add_star, "*", "")),
      p_value = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
      Outcome = outcome_label
    ) %>%
    select(Outcome, OR_CI, p_value)
}

# 4. Process binary outcome model (e.g., logistic)
process_binary_model <- function(model, outcome_label, term_name = "allocation_fIntervention") {
  abs_risk <- get_abs_risk(model, outcome_label)
  risk_diff <- get_risk_difference(model, outcome_label)
  or <- get_or_irr(model, outcome_label, term_name = term_name)
  
  bind_rows(abs_risk, risk_diff) %>%
    pivot_wider(names_from = Group, values_from = Risk_CI) %>%
    rename(
      Control = Control,
      Treatment = Intervention,
      `Risk Difference` = starts_with("RD")
    ) %>%
    select(Outcome, Control, Treatment, `Risk Difference`) %>%
    left_join(or, by = "Outcome")
}

# 5. Process count outcome model (e.g., Poisson)
process_count_model <- function(model, outcome_label, newdata, term_name = "allocation_fIntervention") {
  pred <- predictions(model, newdata = newdata, type = "response") %>%
    as.data.frame() %>%
    mutate(
      estimate = round(estimate, 2),
      conf.low = round(conf.low, 2),
      conf.high = round(conf.high, 2),
      Group = allocation_f,
      Count_CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
      Outcome = outcome_label
    ) %>%
    select(Outcome, Group, Count_CI)
  
  diff <- comparisons(model, newdata = newdata, type = "response") %>%
    as.data.frame() %>%
    slice(1) %>%
    mutate(
      estimate = round(estimate, 2),
      conf.low = round(conf.low, 2),
      conf.high = round(conf.high, 2),
      Group = paste0("RD (", contrast, ")"),
      Count_CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
      Outcome = outcome_label
    ) %>%
    select(Outcome, Group, Count_CI)
  
  summary_table <- bind_rows(pred, diff) %>%
    pivot_wider(names_from = Group, values_from = Count_CI) %>%
    rename(
      Control = Control,
      Treatment = Intervention,
      `Risk Difference` = starts_with("RD")
    ) %>%
    select(Outcome, Control, Treatment, `Risk Difference`)
  
  irr <- get_or_irr(model, outcome_label, term_name = term_name, add_star = TRUE)
  
  left_join(summary_table, irr, by = "Outcome")
}

# 6. Process outbreak size model (bootstrap, home-weighted)
process_outbreak_size_model <- function(model, data, outcome_label,
                                        term_name = "allocation_fIntervention",
                                        n_boot = 2000, seed = 123) {
  set.seed(seed)
  
  # Predict counts per row (ignoring random effects)
  pred_probs <- predict(model, newdata = data, type = "response", re.form = NA)
  total <- data$outbreak_size_secondary + data$noncases
  predicted_counts <- pred_probs * total
  
  df <- data.frame(
    allocation_f = data$allocation_f,
    predicted_count = predicted_counts
  )
  
  # Bootstrap to get 95% CI for each group and risk difference
  boot_results <- replicate(n_boot, {
    df_boot <- df %>%
      group_by(allocation_f) %>%
      sample_frac(replace = TRUE) %>%
      summarise(mean_pred = mean(predicted_count), .groups = "drop")
    
    c(
      Control = df_boot$mean_pred[df_boot$allocation_f == "Control"],
      Treatment = df_boot$mean_pred[df_boot$allocation_f == "Intervention"],
      RD = df_boot$mean_pred[df_boot$allocation_f == "Intervention"] -
        df_boot$mean_pred[df_boot$allocation_f == "Control"]
    )
  }, simplify = "matrix")
  
  # Compute means and CIs
  boot_summary <- apply(boot_results, 1, function(x) {
    mean_val <- round(mean(x), 2)
    ci <- round(quantile(x, c(0.025, 0.975)), 2)
    sprintf("%.3f (%.3f, %.3f)", mean_val, ci[1], ci[2])
  })
  
  # Compute odds ratio for allocation effect (ignoring random effects)
  coef_info <- summary(model)$coefficients
  if(!term_name %in% rownames(coef_info)) stop("term_name not found in model coefficients")
  est <- coef_info[term_name, "Estimate"]
  se  <- coef_info[term_name, "Std. Error"]
  
  or <- round(exp(est), 2)
  or_low <- round(exp(est - 1.96*se), 2)
  or_high <- round(exp(est + 1.96*se), 2)
  OR_CI <- sprintf("%.3f (%.3f, %.3f)", or, or_low, or_high)
  
  # Assemble final table
  final <- data.frame(
    Outcome = outcome_label,
    Control = boot_summary["Control"],
    Treatment = boot_summary["Treatment"],
    `Risk Difference` = boot_summary["RD"],
    OR_CI = OR_CI,
    p_value = as.character(signif(coef_info[term_name, "Pr(>|z|)"], 3)),
    check.names = FALSE
  )
  
  return(final)
}
