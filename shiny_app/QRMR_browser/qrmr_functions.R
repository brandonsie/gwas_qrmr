mr_fun <- function(input_row, rma_method = "EB", tau_index = 5:9, 
                   beta_index = c(11,19,27,35,43),
                   se_index = c(13,21,29,37,45)){
  snp_df <- data.frame(
    tau = input_row[tau_index] %>% as.numeric,
    snp_beta = input_row[beta_index] %>% as.numeric,
    snp_se = input_row[se_index] %>% as.numeric
  )
  snp_df$tau2 <- (snp_df$tau-0.5)^2
  
  snp_rma <- try(
    metafor::rma(
      yi = snp_beta, sei = snp_se,
      mods = cbind(tau, tau2),
      data = snp_df,
      method = rma_method #defualt REML. crashes 49% fisher does not converge
    ) # yi outcome, # sei standard error, # mods moderators, # data input df
  )
  if(class(snp_rma)[1] == "try-error"){
    mr_output <- rep(-9, 12)
  } else{
    mr_output <- c(snp_rma$beta, snp_rma$ci.lb, snp_rma$ci.ub, snp_rma$pval)
  }
  if(length(mr_output) != 12){
    print(paste0("wrong mr output. length ", length(mr_output)))
    print(mr_output)
    mr_output <- rep(-9, 12)
  }
  
  names(mr_output) <- c(
    "MR_b0", "MR_b1", "MR_b2",
    "MR_b0_L95", "MR_b1_L95", "MR_b2_L95",
    "MR_b0_U95", "MR_b1_U95", "MR_b2_U95",
    "MR_b0_pval", "MR_b1_pval", "MR_b2_pval")

  return(mr_output)
  
  
}
  

gg_mr_beta_scatter <- function(qrmr_filter){
  qrmr_filter$phenotype = factor(qrmr_filter$phenotype, levels = c('bmi', 'ht', 'a1c'))
  g <- ggplot(qrmr_filter, aes(x = MR_b1, y = MR_b2 , alpha = log10_min_pval, rsid = rsid, chr = chr)) +
    xlab("QRMR Linear Coefficient") + ylab("QRMR Quadratic Coefficient") +
    geom_point(aes(color = match_filter)) +
    theme_bw() +
    scale_color_discrete(name = "Match Filter") +
    scale_alpha_continuous(guide = FALSE)
    # scale_alpha(range = c(0.1, 0.5), breaks = c(50,100,150, 200), labels = c(50,100,150,200), name = "min(-log10(p))") 
}

gg_manhattan <- function(qrmr_sub, n_sample = NULL){
  if(!is.null(n_sample)){
   qrmr_sub <- qrmr_sub[sample(nrow(qrmr_sub), n_sample),] 
  }
  g <- ggplot(qrmr_sub, aes(x = pos, y = log10_min_pval, rsid = rsid, chr = chr)) +
    geom_point(aes(color = chr)) +
    facet_grid(phenotype~chr, scales = "free_x")
}

# gg_plotdf <- function(qrmrf, linear, this_snp){
#   this_qrmr_12 <- qrmrf[qrmrf$rsid == this_snp,] 
#   plot_df <- data.frame(
#     tau = this_qrmr_12[,grepl("^tau[0-9]$", colnames(this_qrmr_12))] %>% as.numeric %>% sort,
#     qrmr_12_beta = this_qrmr_12[grepl("beta_snp", colnames(this_qrmr_12))]  %>% as.numeric(),
#     qrmr_12_se = this_qrmr_12[grepl("se_snp", colnames(this_qrmr_12))]  %>% as.numeric()
#   )
#   plot_df$qrmr_12_estimate <- this_qrmr_12$MR_b0 + this_qrmr_12$MR_b1*plot_df$tau +
#   this_qrmr_12$MR_b2*(plot_df$tau - 0.5)^2
#   plot_df$qrmr_12_lower <- this_qrmr_12$MR_b0_L95 + this_qrmr_12$MR_b1_L95*plot_df$tau +
#     this_qrmr_12$MR_b2_L95*(plot_df$tau - 0.5)^2
#   plot_df$qrmr_12_upper <- this_qrmr_12$MR_b0_U95 + this_qrmr_12$MR_b1_U95*plot_df$tau +
#     this_qrmr_12$MR_b2_U95*(plot_df$tau - 0.5)^2
#   
#   plot_df
# }

gg_qrmr <- function(qrmrf, linear, this_snp){
  this_qrmr_12 <- qrmrf[qrmrf$rsid == this_snp,] 
  plot_df <- data.frame(
    tau = this_qrmr_12[,grepl("^tau[0-9]$", colnames(this_qrmr_12))] %>% as.numeric %>% sort,
    qrmr_12_beta = this_qrmr_12[grepl("beta_snp", colnames(this_qrmr_12))]  %>% as.numeric(),
    qrmr_12_se = this_qrmr_12[grepl("se_snp", colnames(this_qrmr_12))]  %>% as.numeric()
  )
  plot_df$qrmr_12_estimate <- this_qrmr_12$MR_b0 + this_qrmr_12$MR_b1*plot_df$tau +
    this_qrmr_12$MR_b2*(plot_df$tau - 0.5)^2
  plot_df$qrmr_12_lower <- this_qrmr_12$MR_b0_L95 + this_qrmr_12$MR_b1_L95*plot_df$tau +
    this_qrmr_12$MR_b2_L95*(plot_df$tau - 0.5)^2
  plot_df$qrmr_12_upper <- this_qrmr_12$MR_b0_U95 + this_qrmr_12$MR_b1_U95*plot_df$tau +
    this_qrmr_12$MR_b2_U95*(plot_df$tau - 0.5)^2
  

  
  # ggplot(plot_df) + geom_point(aes(x = tau %>% as.numeric, y = qrmr_12_upper))
  
  # 
  taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  this_g <- ggplot(plot_df, aes(x = tau)) +
    theme_bw() +
    ggtitle(paste0(this_qrmr_12$phenotype, ": ", this_snp),
            subtitle = paste0("b0: ", round(this_qrmr_12$MR_b0, 3),
                              ", b1: ", round(this_qrmr_12$MR_b1, 3),
                              ", b2: ", round(this_qrmr_12$MR_b2, 3),
                              ", log10p0: ", round(this_qrmr_12$MR_b0_log10p, 3),
                              ", log10p1: ", round(this_qrmr_12$MR_b1_log10p, 3),
                              ", log10p2: ", round(this_qrmr_12$MR_b2_log10p, 3))) +
    #12 covariates
    geom_point(aes(x = tau + 0.01, y = qrmr_12_beta)) +
    geom_errorbar(aes(x = tau + 0.01, ymin = qrmr_12_beta - 1.96*qrmr_12_se, ymax = qrmr_12_beta + 1.96*qrmr_12_se)) +
    geom_ribbon(aes(ymin = qrmr_12_lower, ymax = qrmr_12_upper), alpha = 0.2, fill = 'blue') +
    ylab("Effect Size") +
    theme(text = element_text(size = 15)) +
    scale_x_continuous(breaks = taus, labels = taus)

  if(!is.null(linear)){

    this_linear <- linear[linear$ID == this_snp,]
    this_g <- this_g +
      #linear
      geom_hline(data = this_linear, aes(yintercept = BETA), color = 'red') +
      geom_hline(data = this_linear, aes(yintercept = U95), color = 'red', linetype = 'dashed') +
      geom_hline(data = this_linear, aes(yintercept = L95), color = 'red', linetype = 'dashed')

  }
  
  # shift simulator
  num_obs <- 500000
  phenotype_ref <- data.frame(p = c('bmi', 'ht', 'a1c'),
                              mean = c(27.2, 167.5, 32.9),
                              sd = c(5.57, 15.6, 13.4))
  this_ref <- phenotype_ref %>% dplyr::filter(p == this_linear$phenotype)
  p0 = rnorm(num_obs, mean = this_ref$mean, sd = this_ref$sd)
  q0 = ecdf(p0)(p0)
  n1 <- data.frame(
    genotype = 0,
    phenotype = p0
  )
  
  b0 = this_qrmr_12$MR_b0 #-2.221
  b1 = this_qrmr_12$MR_b1 #3.3
  b2 = this_qrmr_12$MR_b2 #10.5
  
  
  effect <- function(q){b0 + b1*q + b2*(0.5-q)^2}
  p1 = p0 + effect(q0)
  q1 = ecdf(p1)(p1)
  p2 = p1 + effect(q1)
  
  
  n1 <- dplyr::bind_rows(
    n1,
    data.frame(genotype = 1, phenotype = p1),
    data.frame(genotype = 2, phenotype = p2)
  )
  
  n1$genotype <- n1$genotype %>% as.factor
  
  
  g_sim <- ggplot(n1) + 
    theme_bw() +
    geom_histogram(aes(x = phenotype, fill = genotype), 
                   bins =  50, position = 'dodge') +
    facet_grid(genotype~.) +
    ggtitle("Simulated Quantile-Dependent Distribution Shift") +
  ylab('density') +
    xlab(this_linear$phenotype) +
    theme(
      text = element_text(size = 15),
      # axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = c(0.9, 0.85)
    )
  
  
  # merge
  gridExtra::grid.arrange(this_g, g_sim, ncol = 2)
}


