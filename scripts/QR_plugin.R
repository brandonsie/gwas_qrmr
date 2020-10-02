#' R plugin to perform quantile regression and standard error computation with 
#' plink 1.9 / 1.07 with quantreg::rq and quantreg::summary.rq
#' 
#' @param PHENO vector of length n with quantitative trait scores
#' @param GENO matrix of n observation by x variants. 
#' @param CLUSTER vector of length n with cluster identities
#' @param COVAR matrix of n observation by p covariates. If no covar, == FALSE
#' 
#' @return Matrix 46 row by x col of 46 quantreg outputs per variant.

Rplink <- function(PHENO, GENO, CLUSTER, COVAR){

  # Initial setup 
  `%>%` <- magrittr::`%>%` 
  taus <- c(0.1, 0.25, 0.5, 0.75, 0.9) # quantiles to report
  
  # identify whether covariates were passed by plink
  has_covar <- ifelse(class(COVAR) == "logical", FALSE, TRUE) 
  
  qr_fun <- function(x){ # qr_fun is evaluated for each variant
    debug <- FALSE # if(debug) print extra information about intermediate steps
    base_q2 <- rep(-9, 45) # If quantreg fails, output vector of -9  
    
    # exclude observations with missing genotype
    sub_gene <- PHENO[!is.na(x)]
    sub_snp <- x[!is.na(x)]
    if(has_covar){
      sub_cov <- COVAR[!is.na(x),] %>% as.data.frame 
      colnames(sub_cov) <- paste0("c", (1:ncol(sub_cov)))
    } else{
      sub_cov <- NULL
    }
    
    # check unique genotype values. cannot run quantreg::rq if 100% homogenous
    unique_x <- x[!is.na(x)] %>% unique
    if(length(unique_x) == 1){
      if(debug) print("only one unique x")
      return(c(length(base_q2), base_q2))
    } else{
      if(debug) print(table(x))
    }
    
    # initialize data frame for input into quantreg::rq
    qr_df <- data.frame(x = x, y = PHENO)[!is.na(x),]
    if(has_covar){qr_df <- dplyr::bind_cols(qr_df, sub_cov)}
    
    # Prepare regression eqn. `y ~ x + c1 + c2 + ... cp` for p covariates
    if(has_covar){
      this_formula <- as.formula(
        paste0("y ~ x+", paste0(colnames(sub_cov), collapse = "+"))
      )
    } else{
      this_formula <- as.formula("y ~ x")
    } 
      
    # quantreg::rq
    rq_method = "pfn" #quantile regression method. ?quantreg::rq
    q1 <- try(
      quantreg::rq(this_formula, tau = taus, data = qr_df, method = rq_method)
    )
    
    # quantreg::summary.rq to estimate standard error
    if(class(q1)[1] == "try-error"){
      if(debug) print("returning from error")
      q2 <- base_q2
    } else{
      if(debug) print("no error")
        
      make_summary <- TRUE
      se_method <- "nid" # standard error method. ?quantreg::summary.rq  #sfn=boot? BLB?
      if(make_summary){
        if(debug) print(paste0("Start summary:", Sys.time()))
        sum1 <- summary(q1, se = se_method)
        q2 <- c(taus, sum1 %>% lapply(function(x) x$coefficients[1:2,] %>% as.numeric) %>% unlist)
        # above subsetting keeps only intercept and snp estimates. discards covariate estimates.
        # output fields: 
        # taus(4). beta_intercept. beta_x. lower_intercept. lower_beta. 
        # upper_intercept. upper_beta. next tau beta and x.
          
      } else{
        q2 <- c(taus, q1$coefficients[1:2,] %>% as.numeric) 
        # taus, then intercept and x for each tau. no bounds
      }
    }
    
    # plink expects output per variant as a vector of the form c(length(r), r)
    q3 <- c(length(q2), q2)
    if(debug) print(paste("q3 length:", length(q3)))
    if(debug) print(head(q3))
    q3
  } # end qr_fun()
 
  use_future = TRUE
  if(use_future){
    # call qr fun for each variant
    slurm_cores <- min(future::availableCores(""), 20)
    print(paste("used cores: ", slurm_cores, ", ", Sys.time()))
    future::plan(future::multisession, workers = slurm_cores)
    
    q_out <- furrr::future_map(GENO %>% as.data.frame, qr_fun)
  } else{
    q_out <- purrr::map(GENO %>% as.data.frame, qr_fun)  
    
  }
 
  q_out <- q_out %>% unlist %>% matrix(nrow = 46) #%>% as.numeric
  q_out
  
} # end Rplink()

      
