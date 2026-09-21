# Project : PRISMA early pregnancy HbA1c and adverse maternal and newborn outcomes
# Description: This project is to assess the relationships between early HbA1c and adverse maternal and newborn outcomes
# Author: Wen-Chien Yang 
# updated: 2026-09-11
      
  ####**********************************************************************
  ####                              Content                             ####
  ####********************************************************************** 
  # I.  Descriptive statistics                                     line 24
  # II. Descriptive statistics by site                             line 103
  # III. HbA1c histogram                                           line 208
  # IV.  Primary analysis (regression and MA)                      line 246
  # V.   Sensitivity analysis, among A1c < 5.7%                    line 496     
  # VI.  Sensitivity analysis, among no anemia                     line 671
  # VII. Sensitivity analysis, among mild/moderate anemia          line 841
  # VIII.Sensitivity analysis, excluding late A1c measurement      line 1011
  # IX.  Sensitivity analysis, using rule in and rule out cutoffs  line 1184 
  # X.   GAM : site specific curve                                 line 1304  
  # XI.  ROC analysis among overall population                     line 1529
  # XII. ROC analysis by site                                      line 1631
  # XIII.AUC figure                                                line 1800
  # XIV. Site specific estimates                                   line 1852
      
#*******************************************************************************
#                        I. Descriptive statistics ----   
#*******************************************************************************
  ## Create variable list for loop 
  vars_con <- c("MAT_AGE",   "SCHOOL_yrs", "WEIGHT_ENROLL", "HEIGHT_ENROLL", "BMI", "A1C", "BOE_GA_WKS_ENROLL", "HBA1C_GA_WKS", "OGTT_GA_WKS")
  vars_cat <- c("AGE_GROUP", "SCHOOL_GROUP", "MULTIPARITY", "EMERGENT_CS", "SMOKE", "PREVPREG_GDM", "PREVPREG_MACRO", "BMI4CAT")

    ## 1) Statistics for continuous variables ----
        for (var in vars_con) {
          
          cat("\nSummary statistics for:", var, "\n")
          
          n_val    <- sum(!is.na(d[[var]]))
          mean_val <- mean(d[[var]], na.rm = TRUE)
          sd_val   <- sd(d[[var]], na.rm = TRUE)
          
          cat("N (non-missing):", n_val, "\n")
          cat("Mean:", round(mean_val, 1), "\n")
          cat("SD:", round(sd_val, 1), "\n")
          
          cat(paste0(round(mean_val,1), " (", round(sd_val,1), "), n=", n_val), "\n")}
      
    ## 2) Statistics for categorical variables ----
        for (var in vars_cat) {
          
          cat("\n", var, "\n")
          x <- d[[var]][!is.na(d[[var]])]
          denom <- length(x)
          tab <- table(x)
    
          for (cate in names(tab)) {
            
            n_cate <- tab[[cate]] 
            pct    <- round(n_cate / denom * 100, 1)
            
            cat(cate, ": ", paste0(n_cate, "/", denom, " (", pct, "%)"), "\n")}} 

    ## 3) Statistics for continuous variables by GDM status ----   
     con_table <- d%>%  filter(!is.na(GDM)) %>%   
        pivot_longer(
          cols = all_of(vars_con),
          names_to = "Variable",
          values_to = "Value"
        ) %>%
        group_by(Variable, GDM) %>%
        summarise(
          N = sum(!is.na(Value)),
          Mean = round(mean(Value, na.rm = TRUE), 1),
          SD = round(sd(Value, na.rm = TRUE), 1),
          .groups = "drop"
        ) %>%
        mutate(`Mean (SD)` = sprintf("%.1f (%.1f), n=%d", Mean, SD, N)) %>%
        select(Variable, GDM, `Mean (SD)`) %>%
        pivot_wider(
          names_from = GDM,
          values_from = `Mean (SD)`) 
      
    ## 4) Statistics for categorical variables by GDM status ----
      cat_table <- d %>% filter(!is.na(GDM)) %>%
        
        mutate(across(all_of(vars_cat), as.character)) %>%
        pivot_longer(cols = all_of(vars_cat),
                     names_to = "Variable",
                     values_to = "Category") %>%
    
        filter(!is.na(Category)) %>%
        group_by(Variable, GDM) %>%
        mutate(Denominator = n()) %>%
        group_by(Variable, GDM, Category) %>%
        
        summarise(
          N = n(),
          Denominator = first(Denominator),
          Percent = round(N / first(Denominator) * 100, 1),
          cell = paste0(N, "/", Denominator, " (", Percent, "%)"),
          .groups = "drop") %>% 
        select(Variable, Category, GDM, cell) %>%
        pivot_wider(names_from = GDM, values_from = cell)

#*******************************************************************************
#                    II. Descriptive statistics by site ----  
#*******************************************************************************
    ## 1) Frequency and proportions of GDM by site ---- 
    gdm_site <- d %>%
      filter(!is.na(GDM), !is.na(SITE)) %>%
      count(SITE, GDM)%>%
      group_by(SITE) %>%
      mutate(
        Percent = round(n / sum(n) * 100, 1),
        Statistics = paste0(n, "/", sum(n), " (", Percent, "%)")) %>%
      ungroup()
    
    ## 2) Statistics for continuous variables ----
    cont_table_site <- d %>%
      filter(!is.na(SITE)) %>%
      pivot_longer(
        cols = all_of(vars_con),
        names_to = "Variable",
        values_to = "Value"
      ) %>%
      mutate(Variable = factor(Variable, levels = vars_con)) %>%
      group_by(SITE, Variable) %>%
      summarise(
        N    = sum(!is.na(Value)),
        Mean = round(mean(Value, na.rm = TRUE), 1),
        SD   = round(sd(Value, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>%
      mutate(Summary = sprintf("%.1f (%.1f), n=%d", Mean, SD, N)) %>%
      select(Variable, SITE, Summary) %>%
      pivot_wider(
        names_from = SITE,
        values_from = Summary) %>%
      arrange(Variable)
  
    ## 3) Statistics for categorical variables by site ----
    cate_table_site <- d %>%
      filter(!is.na(SITE)) %>%
      mutate(across(all_of(vars_cat), as.character)) %>%
      pivot_longer(
        cols = all_of(vars_cat),
        names_to = "Variable",
        values_to = "Category"
      ) %>%
      filter(!is.na(Category)) %>%
      group_by(Variable, SITE) %>%
      mutate(Denominator = n()) %>%
      group_by(Variable, SITE, Category) %>%
      summarise(
        N = n(),
        Denominator = first(Denominator),
        Percent = round(N/ Denominator * 100, 1),
        Summary = paste0(N, "/", Denominator, " (", Percent, "%)"),
        .groups = "drop"
      ) %>%
      select(Variable, Category, SITE, Summary) %>%
      pivot_wider(
        names_from  = SITE,
        values_from = Summary)
  
    ## 4) Statistics for continuous variables by site and GDM status ----   
    cont_table_site_gdm <- d %>%
      filter(!is.na(GDM), !is.na(SITE)) %>%
      pivot_longer(
        cols = all_of(vars_con),
        names_to  = "Variable",
        values_to = "Value"
      ) %>%
      group_by(SITE, Variable, GDM) %>%
      summarise(
        N = sum(!is.na(Value)),
        Mean = round(mean(Value, na.rm = TRUE), 1),
        SD = round(sd(Value, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>%
      mutate(`Mean (SD)` = sprintf("%.1f (%.1f), n=%d", Mean, SD, N)) %>%
      select(SITE, Variable, GDM, `Mean (SD)`) %>%
      pivot_wider(
        names_from = GDM,
        values_from = `Mean (SD)`) 
    
    ## 5) Statistics for categorical variables by SITE and by GDM status ----
    cate_table_site_gdm <- d %>%
      filter(!is.na(GDM), !is.na(SITE)) %>%
      mutate(across(all_of(vars_cat), as.character)) %>%
      pivot_longer(
        cols = all_of(vars_cat),
        names_to = "Variable",
        values_to = "Category"
      ) %>%
      filter(!is.na(Category)) %>%
      group_by(SITE, Variable, GDM) %>%
      mutate(Denominator = n()) %>%
      group_by(SITE, Variable, GDM, Category) %>%
      summarise(N = n(),
                Denominator = first(Denominator),
                Percent = round(N / first(Denominator) * 100, 1),
                cell = paste0(N, "/", Denominator, " (", Percent, "%)"),
                .groups = "drop") %>% 
      select(SITE, Variable, Category, GDM, cell) %>%
      pivot_wider(
        names_from = GDM,
        values_from = cell) 
    
#*******************************************************************************
#                              III. HbA1c histogram ----
#*******************************************************************************
    ## 1) Histogram by site and by GDM status ----
    
    # Create two separate datasets for GDM cases and non-GDM cases 
    dt.nogdm <- d%>% filter(GDM == 0)
    dt.gdm <- d%>% filter(GDM == 1)
    
    # Create histogram counts for GDM == 1 to estimate scale ratio
    hist0 <- hist(dt.nogdm$A1C, breaks = seq(0, 7, 0.5), plot = FALSE)
    hist1 <- hist(dt.gdm$A1C, breaks = seq(0, 7, 0.5), plot = FALSE)
    
    # Calculate a ratio to match the scales
    scale_ratio <-  max(hist0$counts) / max(hist1$counts) - 3  
  
    # Create histogram
    fig  <- ggplot() +
      geom_histogram(data = dt.nogdm, aes(x = A1C, fill = "No GDM"),
                     binwidth = 0.5, alpha = 0.45, colour = "black") +
      
      geom_histogram(data = dt.gdm, aes(x = A1C, y = after_stat(count) * scale_ratio, fill = "GDM"),
                     binwidth = 0.5, alpha = 0.45, colour = "black") +
      
      facet_wrap2(~ Site, scales = "free_y",
                    strip.position = "top", axes = "all") + 
      
      scale_x_continuous(breaks = 0:7, limits = c(0, 7)) +
      scale_y_continuous(name = "Count (No GDM)", sec.axis = sec_axis(~ . / scale_ratio, name = "Count (GDM)")) +
      scale_fill_manual(
        name = "GDM status",
        values = c("No GDM" = "#1B9E77", "GDM" = "#D95F02"))+
      
      theme_minimal() +
      labs(x = "HbA1c (%)") + 
      theme(axis.title.y.right = element_text(color = "black"),
            axis.title.y.left  = element_text(color = "black"))

#*******************************************************************************
#                    IV.  Primary analysis (regression and MA) ---- 
#*******************************************************************************

  # Create list for running loop
  SITES <-c("India-CMC", "India-SAS", "Pakistan", "Kenya", "Zambia") 
  OUTCOMES <- c("GDM", "HDP", "EMERGENT_CS", "PRETERM", "LGA", "STILLBIRTH")

  ## 1) model 1: HbA1c ----
    Result_model1 <- list ()      # create list to store meta-analysis result for each outcome
    
    for (outcome in OUTCOMES){    # A loop for each outcome in OUTCOMES list
      
      res_1 <- list()             # To store the estimate for each site for a specific outcome  
      
      for (site in SITES) {       # Within a specific outcome, a loop for each site 
        
        f <- as.formula(paste(outcome, "~ a1c"))
        
        model <- glm(f, data = subset(d, SITE == site), family = poisson(link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) # get robust variances for modified Poisson
        
        N_obs <- nobs(model)                                              # get the N of observation for each regression
        
        res_1[[site]] <- data.frame(
                         site  = site,
                         N_obs = N_obs,
                         logRR = logRR,
                         SE    = SE)}
                      
      Res_1 <- do.call(rbind, res_1)          # rowbind site-specific results, Res_1 is a table of 5 rows (sites) with logRR and SE 
      
      meta_model1 <- rma(yi   = Res_1$logRR,  # run meta-analysis from the metafor pacakge, yi : the effect size estimate  
                         sei  = Res_1$SE,     # sei : standard error   
                         data = Res_1,
                         method = "REML")     # use restricted maximum likelihood
      
      print(meta_model1)
      
      pooled_RR_model1  <- exp(coef(meta_model1))
      pooled_LCI_model1 <- exp(meta_model1$ci.lb)
      pooled_UCI_model1 <- exp(meta_model1$ci.ub)
      k_sites           <- meta_model1$k
      tau2              <- meta_model1$tau2
      I2                <- meta_model1$I2
      N_total_obs       <- sum(Res_1$N_obs, na.rm = TRUE)  # Total N of observations across sites for this outcome regression

      Result_model1[[outcome]] <- data.frame(
        Outcome = outcome,
        RR    = pooled_RR_model1, 
        CI    =  paste0(round(pooled_LCI_model1, 2), "-", round(pooled_UCI_model1, 2)),
        p     =  ifelse(meta_model1$pval < 0.001, "<0.001", sprintf("%.3f", meta_model1$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
     Result_Model1 <- do.call(rbind, Result_model1)
        
    ## create MA results for all outcomes for model 1 - only HbA1c 
    MA_result1<- cbind(Result_Model1[1,], Result_Model1[2, ], Result_Model1[3,], Result_Model1[4,] , Result_Model1[5, ], Result_Model1[6,] )

  ## 2) model 2: HbA1c + age + BMI + EDUCATION + MULTIPARITY + HBA1C_GA ----
    # The setup for model 2 is the same as model 1, but model 2 includes more covariates in regression 
    
    Result_model2 <- list ()
    
    for (outcome in OUTCOMES){
      res_2 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS"))
        model <- glm(f, data = subset(d, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) # get robust variances for modified Poisson
        N_obs <- nobs(model) 
        
        res_2[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_2 <- do.call(rbind, res_2)
      
      meta_model2 <- rma(yi = Res_2$logRR,
                         sei = Res_2$SE,
                         data = Res_2,
                         method = "REML")
      
      print(meta_model2)
      
      pooled_RR_model2  <- exp(coef(meta_model2))
      pooled_LCI_model2 <- exp(meta_model2$ci.lb)
      pooled_UCI_model2 <- exp(meta_model2$ci.ub)
      k_sites           <- meta_model2$k
      tau2              <- meta_model2$tau2
      I2                <- meta_model2$I2 
      N_total_obs       <- sum(Res_2$N_obs, na.rm = TRUE) 
      
      Result_model2[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model2,
        CI      =  paste0(round(pooled_LCI_model2, 2), "-", round(pooled_UCI_model2, 2)),
        p       = ifelse(meta_model2$pval < 0.001, "<0.001", sprintf("%.3f", meta_model2$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2      = paste0((round(I2)), "%"))}
    
    Result_Model2 <- do.call(rbind, Result_model2)
    MA_result2<- cbind(Result_Model2[1,] , Result_Model2[2, ], Result_Model2[3,], Result_Model2[4,] , Result_Model2[5, ], Result_Model2[6,] )
  
  ## 3) model 3: HbA1c + age + BMI + EDUCATION + HBA1C_GA + previous GDM (among those with prior preg) ----
    # Create a dataset for those with prior pregnancies
    d_1 <- subset(d, MULTIPARITY == 1) 
    
    Result_model3 <- list ()
    
    for (outcome in OUTCOMES){
      res_3 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM"))
        model <- glm(f, data = subset(d_1, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_3[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_3 <- do.call(rbind, res_3)
      
      meta_model3 <- rma(yi = Res_3$logRR,
                         sei = Res_3$SE,
                         data = Res_3,
                         method = "REML")
      
      print(meta_model3)
      
      pooled_RR_model3  <- exp(coef(meta_model3))
      pooled_LCI_model3 <- exp(meta_model3$ci.lb)
      pooled_UCI_model3 <- exp(meta_model3$ci.ub)
      k_sites           <- meta_model3$k
      tau2              <- meta_model3$tau2
      I2                <- meta_model3$I2
      N_total_obs       <- sum(Res_3$N_obs, na.rm = TRUE) # Total N of observations across sites for this outcome
      
      Result_model3[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model3, 
        CI      =  paste0(round(pooled_LCI_model3, 2), "-", round(pooled_UCI_model3, 2)),
        p       = ifelse(meta_model3$pval < 0.001, "<0.001", sprintf("%.3f", meta_model3$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
    Result_Model3 <- do.call(rbind, Result_model3)
    MA_result3<- cbind(Result_Model3[1,] , Result_Model3[2, ], Result_Model3[3,], Result_Model3[4,] , Result_Model3[5, ], Result_Model3[6,] )
    
    ## combine results from model 1, 2, 3 
    HBA1C_RR_MA<-rbind (MA_result1, MA_result2, MA_result3) 

  ## ** Forestplot ----
    # Create dataset for forestplot 
    Result_Model1$Model <- "Model 1"
    Result_Model2$Model <- "Model 2"
    Result_Model3$Model <- "Model 3"
    
    DF <- rbind(as.data.frame(Result_Model1), as.data.frame(Result_Model2), as.data.frame(Result_Model3))
    
    # Rename varaibles
    DF$Outcome[DF$Outcome == "EMERGENT_CS"] <- "emergent CS"
    DF$Outcome[DF$Outcome == "PRETERM"]     <- "preterm birth"
    DF$Outcome[DF$Outcome == "STILLBIRTH"]  <- "stillbirth"
    
    # Reorder outcome order 
    outcome_order <- c("GDM", "HDP", "emergent CS", "preterm birth", "LGA", "stillbirth")
    DF$Outcome <- factor(DF$Outcome, levels = outcome_order)
    
    # Reorder model order
    DF$Model <- factor(DF$Model, levels = c("Model 3", "Model 2", "Model 1"))
    
    # Create a column for results (RR + CI) called label 
    DF <- DF %>% mutate(label = paste0 (sprintf("%.2f (%.2f-%.2f)", RR, LCI, UCI)))
      
    # Prevent figure from overlapping  
    pd <- position_dodge(width = 0.7) 
    
    # Create forestplot using DF dataset 
    forestplot <- ggplot(DF, aes(x = RR, y = fct_rev(Outcome), color = Model, group = Model)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      
      geom_errorbarh(
        aes(xmin = LCI, xmax = UCI),
        height = 0.18,
        size = 0.7,
        position = pd) +
      
      geom_point(size = 4, position = pd) +
      
      geom_text(
        data = DF,
        aes(x = 1.19, y = fct_rev(Outcome), label = label, group = Model),
        inherit.aes = FALSE,
        color = "black",
        hjust = 0,
        size = 4.5,
        position = pd) +
      
      scale_color_manual(
        values = c(
          "Model 1" = "#1B9E77",
          "Model 2" = "#D95F02",
          "Model 3" = "#7570B3" ),
        
        breaks = c("Model 1", "Model 2", "Model 3"),
        labels = c(
          "Model 1 included HbA1c",
          "Model 2 included HbA1c, age groups, BMI categories, education years, GA at HbA1c measurement (weeks), and multiparity status",
          "Model 3 included HbA1c, age groups, BMI categories, education years, GA at HbA1c measurement (weeks), and previous GDM status (among multiparous women)" )) +
      
      scale_x_continuous(limits = c(0.96, 1.20),
                         breaks = seq(0.95, 1.20, by = 0.05)) +

      labs(x = "Risk ratio of adverse outcomes (per 0.1% increase in HbA1c)",
           y = NULL,
           color = NULL) +
      
      coord_cartesian(clip = "off") +
      guides(color = guide_legend(ncol = 1)) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(color = "gray85"),
            axis.text.y = element_text(size = 14, color = "black"),
            axis.text.x = element_text(color = "gray30"),
            legend.position = "bottom",
            legend.text = element_text(size = 11),
            plot.margin = margin(5.5, 150, 5.5, 5.5)) 
    
#*******************************************************************************
#                 V.  Sensitivity analysis, among A1c < 5.7% ---- 
#*******************************************************************************
  # create d_2 for model 6-10 for those with HbA1c < 5.7%
  d_2 <- d %>% subset (HBA1C_PRCNT < 5.7)
  nrow(d_2) # 10409

  ## 1) model 4: HbA1c ----
    Result_model4 <- list ()      # to store meta-analysis result for each outcome
    
    for (outcome in OUTCOMES){    # A loop for each outcome in OUTCOMES
      
      res_4 <- list()
      for (site in SITES) {       # Within each outcome, a loop for each site 
        
        f <- as.formula(paste(outcome, "~ a1c"))
        
        model <- glm(f, data = subset(d_2, SITE == site),
                        family = poisson(link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_4[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_4 <- do.call(rbind, res_4)  
      
      meta_model4 <- rma(yi = Res_4$logRR,  
                         sei = Res_4$SE,
                         data = Res_4,
                         method = "REML")
      
      print(meta_model4)
      
      pooled_RR_model4  <- exp(coef(meta_model4))
      pooled_LCI_model4 <- exp(meta_model4$ci.lb)
      pooled_UCI_model4 <- exp(meta_model4$ci.ub)
      k_sites           <- meta_model4$k
      tau2              <- meta_model4$tau2
      I2                <- meta_model4$I2
      N_total_obs       <- sum(Res_4$N_obs, na.rm = TRUE) # Total N of observations across sites for this outcome
    
      Result_model4[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model4, 
        CI      =  paste0(round(pooled_LCI_model4, 3), "-", round(pooled_UCI_model4, 3)),
        p       = ifelse(meta_model4$pval < 0.001, "<0.001", sprintf("%.3f", meta_model4$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2      = paste0((round(I2)), "%"))}
     
      Result_Model4 <- do.call(rbind, Result_model4)
    
          ## create MA results for all outcomes for model 4 - HbA1c (among < 5.7)
          MA_result4<- cbind(Result_Model4[1,] , Result_Model4[2, ], Result_Model4[3,], Result_Model4[4,] , Result_Model4[5, ], Result_Model4[6,] )
          
  ## 2) model 5: HbA1c + age + BMI + EDUCATION + MULTIPARITY + HBA1C_GA ----
    
    Result_model5 <- list ()
    
    for (outcome in OUTCOMES){
      res_5 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS"))
        model <- glm(f, data = subset(d_2, SITE == site),  # use d_2, the dataset for those with prior pregnancies
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"])  
        N_obs <- nobs(model)
        
        res_5[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_5 <- do.call(rbind, res_5)
      
      meta_model5 <- rma(yi = Res_5$logRR,
                         sei = Res_5$SE,
                         data = Res_5,
                         method = "REML")
      
      summary(meta_model5)
      
      pooled_RR_model5 <- exp(coef(meta_model5))
      pooled_LCI_model5 <- exp(meta_model5$ci.lb)
      pooled_UCI_model5 <- exp(meta_model5$ci.ub)
      k_sites           <- meta_model5$k
      tau2              <- meta_model5$tau2
      I2                <- meta_model5$I2
      N_total_obs       <- sum(Res_5$N_obs, na.rm = TRUE)  
      
      Result_model5[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model5, 
        CI      =  paste0(round(pooled_LCI_model5, 2), "-", round(pooled_UCI_model5, 2)),
        p       = ifelse(meta_model5$pval < 0.001, "<0.001", sprintf("%.3f", meta_model5$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2      = paste0((round(I2)), "%"))}
  
        Result_Model5 <- do.call(rbind, Result_model5)
        MA_result5<- cbind(Result_Model5[1,] , Result_Model5[2, ], Result_Model5[3,], Result_Model5[4,] , Result_Model5[5, ], Result_Model5[6,] )
        
  ## 3) model 6: HbA1c + age + BMI + EDUCATION + HBA1C_GA + previous GDM (among those with prior preg)----
    # Create a dataset of those <5.7 and with prior pregnancies
    d_3 <- subset(d_2, MULTIPARITY == 1) 
    nrow(d_2)
    nrow(d_3)
    
    Result_model6 <- list ()
    
    for (outcome in OUTCOMES){
      res_6 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM"))
        model <- glm(f, data = subset(d_3, SITE == site), # use d_3, <.5.7 and with prior pregnancies
                        family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
        N_obs <- nobs(model)
        
        res_6[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_6 <- do.call(rbind, res_6)
      
      meta_model6 <- rma(yi = Res_6$logRR,
                         sei = Res_6$SE,
                         data = Res_6,
                         method = "REML")
      
      print(meta_model6)
      
      pooled_RR_model6 <- exp(coef(meta_model6))
      pooled_LCI_model6 <- exp(meta_model6$ci.lb)
      pooled_UCI_model6 <- exp(meta_model6$ci.ub)
      k_sites           <- meta_model6$k
      tau2              <- meta_model6$tau2
      I2                <- meta_model6$I2
      N_total_obs       <- sum(Res_6$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model6[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model6, 
        CI      = paste0(round(pooled_LCI_model6, 2), "-", round(pooled_UCI_model6, 2)),
        p       = ifelse(meta_model6$pval < 0.001, "<0.001", sprintf("%.3f", meta_model6$pval)),
        n_studies = k_sites, 
        N_obs   = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
    Result_Model6 <- do.call(rbind, Result_model6)
    MA_result6<- cbind(Result_Model6[1,] , Result_Model6[2, ], Result_Model6[3,], Result_Model6[4,] , Result_Model6[5, ], Result_Model6[6,] )
    
    # combine models from model 4 5 6 
    HBA1C_5.7_RR_MA<-rbind (MA_result4, MA_result5, MA_result6) 
    
    # setwd ("D:/Users/wenchienyang/Documents/Wen-Chien/GDM/GDM_early HbA1c/figures and tables")
    # write_xlsx (HBA1C_5.7_RR_MA, "HBA1C_RR_sensitivity_5.7.xlsx") 
 
#*******************************************************************************
#               VI.  Sensitivity analysis, among no anemia ---- 
#*******************************************************************************
  # create d_noanemia for those without anemia 
  table(d$ANE)
  d_noanemia <- d %>% subset (ANE == "No anemia")
  nrow(d_noanemia) # 7187
    
  ## 1) model 7: HbA1c ----
    Result_model7 <- list ()     
    
    for (outcome in OUTCOMES){    
      
      res_7 <- list()
      
      for (site in SITES) {       
        
        f <- as.formula(paste(outcome, "~ a1c"))
        
        model <- glm(f, data = subset(d_noanemia, SITE == site), # use no anemia dataset
                     family = poisson(link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_7[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_7 <- do.call(rbind, res_7)  
      
      meta_model7 <- rma(yi = Res_7$logRR,   
                         sei = Res_7$SE,
                         data = Res_7,
                         method = "REML")
      
      print(meta_model7)
      
      pooled_RR_model7  <- exp(coef(meta_model7))
      pooled_LCI_model7 <- exp(meta_model7$ci.lb)
      pooled_UCI_model7 <- exp(meta_model7$ci.ub)
      k_sites           <- meta_model7$k
      tau2              <- meta_model7$tau2
      I2                <- meta_model7$I2
      N_total_obs       <- sum(Res_7$N_obs, na.rm = TRUE) 
      
      Result_model7[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model7, 
        CI      = paste0(round(pooled_LCI_model7, 2), "-", round(pooled_UCI_model7, 2)),
        p       = ifelse(meta_model7$pval < 0.001, "<0.001", sprintf("%.3f", meta_model7$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2        = paste0((round(I2)), "%"))}
    
        Result_Model7 <- do.call(rbind, Result_model7)
        
        ## create MA results for all outcomes for model 11 
        MA_result7<- cbind(Result_Model7[1,] , Result_Model7[2, ], Result_Model7[3,], Result_Model7[4,] , Result_Model7[5, ], Result_Model7[6,] )
        
  ## 2) model 8: HbA1c + age + BMI + EDUCATION + MULTIPARITY + HBA1c_GA ----
    Result_model8 <- list ()
    
    for (outcome in OUTCOMES){
      res_8 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS"))
        model <- glm(f, data = subset(d_noanemia, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
        N_obs <- nobs(model)
        
        res_8[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_8 <- do.call(rbind, res_8)
      
      meta_model8 <- rma(yi = Res_8$logRR,
                          sei = Res_8$SE,
                          data = Res_8,
                          method = "REML")
      
      summary(meta_model8)
      
      pooled_RR_model8 <- exp(coef(meta_model8))
      pooled_LCI_model8 <- exp(meta_model8$ci.lb)
      pooled_UCI_model8 <- exp(meta_model8$ci.ub)
      k_sites           <- meta_model8$k
      tau2              <- meta_model8$tau2
      I2                <- meta_model8$I2
      N_total_obs       <- sum(Res_8$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model8[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model8, 
        CI      =  paste0(round(pooled_LCI_model8, 2), "-", round(pooled_UCI_model8, 2)),
        p       = ifelse(meta_model8$pval < 0.001, "<0.001", sprintf("%.3f", meta_model8$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
    Result_Model8 <- do.call(rbind, Result_model8)
    MA_result8<- cbind(Result_Model8[1,] , Result_Model8[2, ], Result_Model8[3,], Result_Model8[4,] , Result_Model8[5, ], Result_Model8[6,] )
    
  ## 3) model 9: HbA1c + age + BMI + EDUCATION + HBA1C_GA + previous GDM (among those with prior preg) ----
    # Create a dataset of those without anemia and with prior pregnancies
    d_noanemia_multiparous <- subset(d_noanemia, MULTIPARITY == 1) 
    Result_model9 <- list ()
    
    for (outcome in OUTCOMES){
      res_9 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM"))
        model <- glm(f, data = subset(d_noanemia_multiparous, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_9[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_9 <- do.call(rbind, res_9)
      
      meta_model9 <- rma(yi = Res_9$logRR,
                         sei = Res_9$SE,
                         data = Res_9,
                         method = "REML")
      
      print(meta_model9)
      
      pooled_RR_model9  <- exp(coef(meta_model9))
      pooled_LCI_model9 <- exp(meta_model9$ci.lb)
      pooled_UCI_model9 <- exp(meta_model9$ci.ub)
      k_sites           <- meta_model9$k
      tau2              <- meta_model9$tau2
      I2                <- meta_model9$I2
      N_total_obs       <- sum(Res_9$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model9[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model9, 
        CI      = paste0(round(pooled_LCI_model9, 2), "-", round(pooled_UCI_model9, 2)),
        p       = ifelse(meta_model9$pval < 0.001, "<0.001", sprintf("%.3f", meta_model9$pval)),
        n_studies = k_sites, 
        N_obs   = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
    Result_Model9 <- do.call(rbind, Result_model9)
    MA_result9<- cbind(Result_Model9[1,] , Result_Model9[2, ], Result_Model9[3,], Result_Model9[4,] , Result_Model9[5, ], Result_Model9[6,] )
    
    # combine models from model 7 8 9 
    HBA1C_NO_ANEMIA<-rbind (MA_result7, MA_result8, MA_result9) 

#*******************************************************************************
#             VII. Sensitivity analysis, among mild/moderate anemia ---- 
#*******************************************************************************
  # create d_anemia for those with mild and moderate anemia 
  table(d$ANE)
  d_anemia <- d %>% subset (ANE == "Mild anemia"|ANE == "Moderate anemia")
  nrow(d_anemia) # 3731
    
  ## 1) model 10: HbA1c ----
    Result_model10 <- list ()     
    
    for (outcome in OUTCOMES){    
      
      res_10 <- list()
      
      for (site in SITES) {       
        
        f <- as.formula(paste(outcome, "~ a1c"))
        
        model <- glm(f, data = subset(d_anemia, SITE == site), # use no anemia dataset
                     family = poisson(link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_10[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_10 <- do.call(rbind, res_10)  
      
      meta_model10 <- rma(yi = Res_10$logRR,   
                         sei = Res_10$SE,
                         data = Res_10,
                         method = "REML")
      
      print(meta_model10)
      
      pooled_RR_model10  <- exp(coef(meta_model10))
      pooled_LCI_model10 <- exp(meta_model10$ci.lb)
      pooled_UCI_model10 <- exp(meta_model10$ci.ub)
      k_sites           <- meta_model10$k
      tau2              <- meta_model10$tau2
      I2                <- meta_model10$I2
      N_total_obs       <- sum(Res_10$N_obs, na.rm = TRUE) 
      
      Result_model10[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model10, 
        CI      = paste0(round(pooled_LCI_model10, 3), "-", round(pooled_UCI_model10, 3)),
        p       = ifelse(meta_model10$pval < 0.001, "<0.001", sprintf("%.3f", meta_model10$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
    Result_Model10 <- do.call(rbind, Result_model10)
    
        # create MA results for all outcomes for model 11 
        MA_result10<- cbind(Result_Model10[1,] , Result_Model10[2, ], Result_Model10[3,], Result_Model10[4,] , Result_Model10[5, ], Result_Model10[6,] )
        
  ## 2) model 11: HbA1c + age + BMI + EDUCATION + MULTIPARITY + HBA1C_GA ----
    Result_model11 <- list ()
    
    for (outcome in OUTCOMES){
      res_11 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS"))
        model <- glm(f, data = subset(d_anemia, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
        N_obs <- nobs(model)
        
        res_11[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_11 <- do.call(rbind, res_11)
      
      meta_model11 <- rma(yi = Res_11$logRR,
                         sei = Res_11$SE,
                         data = Res_11,
                         method = "REML")
      
      summary(meta_model11)
      
      pooled_RR_model11 <- exp(coef(meta_model11))
      pooled_LCI_model11 <- exp(meta_model11$ci.lb)
      pooled_UCI_model11 <- exp(meta_model11$ci.ub)
      k_sites           <- meta_model11$k
      tau2              <- meta_model11$tau2
      I2                <- meta_model11$I2
      N_total_obs       <- sum(Res_11$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model11[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model11, 
        CI      =  paste0(round(pooled_LCI_model11, 2), "-", round(pooled_UCI_model11, 2)),
        p       = ifelse(meta_model11$pval < 0.001, "<0.001", sprintf("%.3f", meta_model11$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2      = paste0((round(I2)), "%"))}
    
    Result_Model11 <- do.call(rbind, Result_model11)
    MA_result11<- cbind(Result_Model11[1,] , Result_Model11[2, ], Result_Model11[3,], Result_Model11[4,] , Result_Model11[5, ], Result_Model11[6,] )
    
  ## 3) model 12: HbA1c + age + BMI + EDUCATION + HBA1C_GA + previous GDM (among those with prior preg) ----
    # Create a dataset of those with anemia and with prior pregnancies
    d_anemia_multiparous <- subset(d_anemia, MULTIPARITY == 1)
    Result_model12 <- list ()
    
    for (outcome in OUTCOMES){
      res_12 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM"))
        model <- glm(f, data = subset(d_anemia_multiparous, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_12[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_12 <- do.call(rbind, res_12)
      
      meta_model12 <- rma(yi = Res_12$logRR,
                         sei = Res_12$SE,
                         data = Res_12,
                         method = "REML")
      
      print(meta_model12)
      
      pooled_RR_model12  <- exp(coef(meta_model12))
      pooled_LCI_model12 <- exp(meta_model12$ci.lb)
      pooled_UCI_model12 <- exp(meta_model12$ci.ub)
      k_sites           <- meta_model12$k
      tau2              <- meta_model12$tau2
      I2                <- meta_model12$I2
      N_total_obs       <- sum(Res_12$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model12[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model12, 
        CI      = paste0(round(pooled_LCI_model12, 2), "-", round(pooled_UCI_model12, 2)),
        p       = ifelse(meta_model12$pval < 0.001, "<0.001", sprintf("%.3f", meta_model12$pval)),
        n_studies = k_sites, 
        N_obs   = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
    Result_Model12 <- do.call(rbind, Result_model12)
    MA_result12<- cbind(Result_Model12[1,] , Result_Model12[2, ], Result_Model12[3,], Result_Model12[4,] , Result_Model12[5, ], Result_Model12[6,] )
    
    # combine models from model 10 11 12 
    HBA1C_ANEMIA<-rbind (MA_result10, MA_result11, MA_result12) 

#*******************************************************************************
#         VIII. Sensitivity analysis, excluding late A1c measurement ----
#*******************************************************************************
  # create d_nolate among HBA1C before GA18, excluding those HBA1C measured >= GA 18
  d_nolate <- d %>% subset (HBA1C_GA_WKS <18)
  nrow(d_nolate) # 9392
    
  ## 1) model 13: HbA1c ----
    Result_model13 <- list ()     
    
    for (outcome in OUTCOMES){    
      
      res_13 <- list()
      
      for (site in SITES) {       
        
        f <- as.formula(paste(outcome, "~ a1c"))
        
        model <- glm(f, data = subset(d_nolate, SITE == site),
                     family = poisson(link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"])  
        N_obs <- nobs(model)
        
        res_13[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_13 <- do.call(rbind, res_13)  
      
      meta_model13 <- rma(yi = Res_13$logRR,  
                          sei = Res_13$SE,
                          data = Res_13,
                          method = "REML")
      
      print(meta_model13)
      
      pooled_RR_model13  <- exp(coef(meta_model13))
      pooled_LCI_model13 <- exp(meta_model13$ci.lb)
      pooled_UCI_model13 <- exp(meta_model13$ci.ub)
      k_sites            <- meta_model13$k
      tau2               <- meta_model13$tau2
      I2                 <- meta_model13$I2
      N_total_obs        <- sum(Res_13$N_obs, na.rm = TRUE) 
      
      Result_model13[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = round(pooled_RR_model13, 2),
        CI      =  paste0(round(pooled_LCI_model13, 2), "-", round(pooled_UCI_model13, 2)),
        p       = ifelse(meta_model13$pval < 0.001, "<0.001", sprintf("%.3f", meta_model13$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2      = paste0((round(I2)), "%")) }
      
    
    Result_Model13 <- do.call(rbind, Result_model13)
    
        # create MA results for all outcomes for model 13 
        MA_result13<- cbind(Result_Model13[1,] , Result_Model13[2, ], Result_Model13[3,], Result_Model13[4,] , Result_Model13[5, ], Result_Model13[6,] )
        
  ## 2) model 14: HbA1c + age + BMI + EDUCATION + MULTIPARITY + HBA1C_GA ----  
    Result_model14 <- list ()
    
    for (outcome in OUTCOMES){
      res_14 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS"))
        model <- glm(f, data = subset(d_nolate, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
        N_obs <- nobs(model)
        
        res_14[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_14 <- do.call(rbind, res_14)
      
      meta_model14 <- rma(yi = Res_14$logRR,
                          sei = Res_14$SE,
                          data = Res_14,
                          method = "REML")
      
      summary(meta_model14)
      
      pooled_RR_model14 <- exp(coef(meta_model14))
      pooled_LCI_model14 <- exp(meta_model14$ci.lb)
      pooled_UCI_model14 <- exp(meta_model14$ci.ub)
      k_sites           <- meta_model14$k
      tau2              <- meta_model14$tau2
      I2                <- meta_model14$I2
      N_total_obs       <- sum(Res_14$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model14[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = round(pooled_RR_model14,2), 
        CI      = paste0(round(pooled_LCI_model14, 2), "-", round(pooled_UCI_model14, 2)),
        p       = ifelse(meta_model14$pval < 0.001, "<0.001", sprintf("%.3f", meta_model14$pval)),
        n_studies = k_sites, 
        N_obs     = N_total_obs,
        I2     = paste0((round(I2)), "%"))}
    
        Result_Model14 <- do.call(rbind, Result_model14)
        MA_result14<- cbind(Result_Model14[1,] , Result_Model14[2, ], Result_Model14[3,], Result_Model14[4,] , Result_Model14[5, ], Result_Model14[6,] )
        
  ## 3) model 15: HbA1c + age + BMI + EDUCATION + HBA1C_GA + previous GDM (among those with prior preg)----
      # Create a dataset of those with no late HbA1c measurement and those with prior pregnancies
      d_nolate_multiparous <- subset(d_nolate, MULTIPARITY == 1) 
      Result_model15 <- list ()
    
    for (outcome in OUTCOMES){
      res_15 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM"))
        model <- glm(f, data = subset(d_nolate_multiparous, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["a1c"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["a1c","a1c"]) 
        N_obs <- nobs(model)
        
        res_15[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_15 <- do.call(rbind, res_15)
      
      meta_model15 <- rma(yi = Res_15$logRR,
                          sei = Res_15$SE,
                          data = Res_15,
                          method = "REML")
      
      print(meta_model15)
      
      pooled_RR_model15 <- exp(coef(meta_model15))
      pooled_LCI_model15 <- exp(meta_model15$ci.lb)
      pooled_UCI_model15 <- exp(meta_model15$ci.ub)
      k_sites           <- meta_model15$k
      tau2              <- meta_model15$tau2
      I2                <- meta_model15$I2
      N_total_obs       <- sum(Res_15$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model15[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = round(pooled_RR_model15, 2),
        CI      = paste0(round(pooled_LCI_model15, 2), "-", round(pooled_UCI_model15, 2)),
        p       = ifelse(meta_model15$pval < 0.001, "<0.001", sprintf("%.3f", meta_model15$pval)),
        n_studies = k_sites, 
        N_obs   = N_total_obs,
        I2     = paste0((round(I2)), "%"))} 
    
    Result_Model15 <- do.call(rbind, Result_model15)
    MA_result15<- cbind(Result_Model15[1,] , Result_Model15[2, ], Result_Model15[3,], Result_Model15[4,] , Result_Model15[5, ], Result_Model15[6,] )
    
    # combine models from model 13 14 15 
    HBA1C_NO_LATE<-rbind (MA_result13, MA_result14, MA_result15) 
    
    setwd ("D:/Users/wenchienyang/Documents/Wen-Chien/GDM/GDM_early HbA1c/figures and tables")
    write_xlsx (HBA1C_NO_LATE, "HBA1C_RR_sensitivity_nolateA1C.xlsx") 
    
#*******************************************************************************
#       IX.  Sensitivity analysis, using rule in and rule out cutoffs ---- 
#*******************************************************************************
      
  # These two cutoffs were derived from the AUC/ROC analyses 
  # Create HbA1c cutoffs using the two threshold approach  
  d$ruleout <- ifelse (d$A1C >= 5.0, ">=5.0", "<5.0" )
  d$rulein  <- ifelse (d$A1C >= 5.5, ">=5.5", "<5.5" )

  ## 1) model 16: rule out cutoff ----
    d$ruleout <- relevel(factor(d$ruleout), ref = ">=5.0")
    
    Result_model16 <- list ()
    
    for (outcome in OUTCOMES){
      res_16 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ ruleout + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS")) # use binary cutoff as exposure in adjusted models
        model <- glm(f, data = subset(d, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["ruleout<5.0"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["ruleout<5.0","ruleout<5.0"]) # robust variances for modified Poisson
        N_obs <- nobs(model)
        
        res_16[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_16 <- do.call(rbind, res_16)
      
      meta_model16 <- rma(yi = Res_16$logRR,
                          sei = Res_16$SE,
                          data = Res_16,
                          method = "REML")
      
      summary(meta_model16)
      
      pooled_RR_model16  <- exp(coef(meta_model16))
      pooled_LCI_model16 <- exp(meta_model16$ci.lb)
      pooled_UCI_model16 <- exp(meta_model16$ci.ub)
      k_sites           <- meta_model16$k
      tau2              <- meta_model16$tau2
      I2                <- meta_model16$I2
      N_total_obs       <- sum(Res_16$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model16[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      = pooled_RR_model16,
        CI      =  paste0(round(pooled_LCI_model16, 2), "-", round(pooled_UCI_model16, 2)),
        p       = ifelse(meta_model16$pval < 0.001, "<0.001", sprintf("%.3f", meta_model16$pval)),
        n_studies = k_sites, 
        N_obs   = N_total_obs,
        I2      = paste0((round(I2)), "%"))}
    
      Result_Model16 <- do.call(rbind, Result_model16)
      MA_result16<- cbind(Result_Model16[1,] , Result_Model16[2, ], Result_Model16[3,], Result_Model16[4,] , Result_Model16[5, ], Result_Model16[6,] )
      
  ## 2) model 17: rule in cutoff ----
    
    Result_model17 <- list ()
    
    for (outcome in OUTCOMES){
      res_17 <- list()
      
      for (site in SITES) {
        
        f <- as.formula(paste(outcome, "~ rulein + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS"))
        model <- glm(f, data = subset(d, SITE == site),
                     family = poisson (link = "log"))
        
        logRR <- coef(model)["rulein>=5.5"]
        SE    <- sqrt(sandwich::vcovHC(model, type = "HC0")["rulein>=5.5","rulein>=5.5"]) 
        N_obs <- nobs(model)
        
        res_17[[site]] <- data.frame(
          site  = site,
          N_obs = N_obs,
          logRR = logRR,
          SE    = SE)}
      
      Res_17 <- do.call(rbind, res_17)
      
      meta_model17 <- rma(yi = Res_17$logRR,
                          sei = Res_17$SE,
                          data = Res_17,
                          method = "REML")
      
      summary(meta_model17)
      
      pooled_RR_model17  <- exp(coef(meta_model17))
      pooled_LCI_model17 <- exp(meta_model17$ci.lb)
      pooled_UCI_model17 <- exp(meta_model17$ci.ub)
      k_sites           <- meta_model17$k
      tau2              <- meta_model17$tau2
      I2                <- meta_model17$I2
      N_total_obs       <- sum(Res_17$N_obs, na.rm = TRUE) ## Total N of observations across sites for this outcome
      
      Result_model17[[outcome]] <- data.frame(
        Outcome = outcome,
        RR      =  pooled_RR_model17, 
        CI      =  paste0(round(pooled_LCI_model17, 3), "-", round(pooled_UCI_model17, 3)),
        p       = ifelse(meta_model17$pval < 0.001, "<0.001", sprintf("%.3f", meta_model17$pval)),
        n_studies = k_sites, 
        N_obs   = N_total_obs,
        I2      = paste0((round(I2)), "%"))}
    
    Result_Model17 <- do.call(rbind, Result_model17)
    MA_result17<- cbind(Result_Model17[1,] , Result_Model17[2, ], Result_Model17[3,], Result_Model17[4,] , Result_Model17[5, ], Result_Model17[6,] )
    
    # combine models from model 16 17 
    HBA1C_cutoff_RR_MA<-rbind (MA_result16, MA_result17) 
    
    setwd ("D:/Users/wenchienyang/Documents/Wen-Chien/GDM/GDM_early HbA1c/figures and tables")
    write_xlsx (HBA1C_cutoff_RR_MA, "HBA1C_RR_sensitivity_cutoff.xlsx")  
  
#*******************************************************************************
#                     X. GAM : site specific curve ---- 
#*******************************************************************************
    # Note that models fit 1-6 have covariates
    # For GAM to generate, need to provide a value for covariates  
      
      # GDM 
      fit1 <- gam(GDM ~ s(A1C, by = Site) + Site + MAT_AGE + BMI + SCHOOL_YRS + HBA1C_GA_WKS + MULTIPARITY, data = d, family = binomial, method = "REML" )
      summary(fit1)

      # HDP
      fit2 <- gam(HDP ~ s(A1C, by = Site)+ Site + MAT_AGE + BMI + SCHOOL_YRS + HBA1C_GA_WKS + MULTIPARITY, data = d, family = binomial, method = "REML")
      summary(fit2)
      
      # emergent CS 
      fit3 <- gam(EMERGENT_CS ~ s(A1C, by = Site) + Site + MAT_AGE + BMI + SCHOOL_YRS + HBA1C_GA_WKS + MULTIPARITY, data = d, family = binomial, method = "REML")
      summary(fit3)
      
      # preterm birth 
      fit4 <- gam(PRETERM ~ s(A1C, by = Site) + Site + MAT_AGE + BMI + SCHOOL_YRS + HBA1C_GA_WKS + MULTIPARITY, data = d, family = binomial, method = "REML")
      summary(fit4) 
      
      # LGA
      fit5 <- gam(LGA ~ s(A1C, by = Site) + Site + MAT_AGE + BMI + SCHOOL_YRS + HBA1C_GA_WKS + MULTIPARITY, data = d, family = binomial, method = "REML")
      summary(fit5)

      # stillbirth
      fit6 <- gam(STILLBIRTH ~ s(A1C, by = Site) + Site + MAT_AGE + BMI + SCHOOL_YRS + HBA1C_GA_WKS + MULTIPARITY, data = d, family = binomial, method = "REML")
      summary(fit6)
      
   # 1) create a prediction dataset (create 200 evenly spaced A1C values between the smallest and largest from our data)
      site_levels <- levels(factor(d$Site))       # all site levels
      
      site_a1c_seq <- d %>%
        group_by(Site) %>%
        reframe(A1C = seq(min(A1C, na.rm = TRUE),
                          max(A1C, na.rm = TRUE),
                          length.out = 200)) 
      
      # NOTE : create dataset new.data that has information that can be used to feed the fit1 to 6 to generate predicted probabilities
        # A1C       = a1c_seq,
        # Site      = site_levels,
        # MAT_AGE   = mean(d$MAT_AGE, na.rm=T),               # mean age, used to represent an average woman
        # BMI       = mean(d$BMI, na.rm=T),                   # mean BMI, used to represent an average woman
        # SCHOOL_YRS  = mean(d$SCHOOL_YRS, na.rm=T),          # mean school years, used to represent an average woman
        # MULTIPARITY = 1,                                    # multiparous = 1 since most women are multiparous
        # HBA1C_GA_WKS = mean(d$HBA1C_GA_WKS, na.rm = TRUE))  # mean HBA1c GA weeks, used to represent an average woman

      site_covariates <- d %>%
        group_by(Site) %>%
        summarise(
          MAT_AGE      = mean(MAT_AGE, na.rm = TRUE),
          BMI          = mean(BMI, na.rm = TRUE),
          SCHOOL_YRS   = mean(SCHOOL_YRS, na.rm = TRUE),
          MULTIPARITY  = 1,
          HBA1C_GA_WKS = mean(HBA1C_GA_WKS, na.rm = TRUE),
          .groups = "drop")
      
      new.data <- site_a1c_seq %>% left_join(site_covariates, by = "Site")  
       
   # 2) predict the probability on the link scale (log-odds) 
      predict.gdm        <- predict(fit1, newdata = new.data, type = "link", se.fit = TRUE)
      predict.hdp        <- predict(fit2, newdata = new.data, type = "link", se.fit = TRUE)
      predict.cs         <- predict(fit3, newdata = new.data, type = "link", se.fit = TRUE)
      predict.preterm    <- predict(fit4, newdata = new.data, type = "link", se.fit = TRUE)
      predict.lga        <- predict(fit5, newdata = new.data, type = "link", se.fit = TRUE)
      predict.stillbirth <- predict(fit6, newdata = new.data, type = "link", se.fit = TRUE)
      
   # 3) Convert to probability using inverse logit
      # Create function inv_logit
      inv_logit <- function(x) exp(x) / (1 + exp(x)) 
      
      # Get THE predicted probability (plan not to run preterm and stillbirth)
      new.data$prob.gdm  <- inv_logit(predict.gdm$fit)
      new.data$upper.gdm <- inv_logit(predict.gdm$fit + 1.96 * predict.gdm$se.fit)
      new.data$lower.gdm <- inv_logit(predict.gdm$fit - 1.96 * predict.gdm$se.fit)
      
      new.data$prob.hdp  <- inv_logit(predict.hdp$fit)
      new.data$upper.hdp <- inv_logit(predict.hdp$fit + 1.96 * predict.hdp$se.fit)
      new.data$lower.hdp <- inv_logit(predict.hdp$fit - 1.96 * predict.hdp$se.fit)
      
      new.data$prob.cs  <- inv_logit(predict.cs$fit)
      new.data$upper.cs <- inv_logit(predict.cs$fit + 1.96 * predict.cs$se.fit)
      new.data$lower.cs <- inv_logit(predict.cs$fit - 1.96 * predict.cs$se.fit)
     
      new.data$prob.lga  <- inv_logit(predict.lga$fit)
      new.data$upper.lga <- inv_logit(predict.lga$fit + 1.96 * predict.lga$se.fit)
      new.data$lower.lga <- inv_logit(predict.lga$fit - 1.96 * predict.lga$se.fit)
      # Remove preterm and stillbirth from figure
      
    # 4) Use ggplot to create site specific curves   
      
      site_cols <- c(
        "South India" = "#173f5f",
        "North India" = "#20639b",
        "Pakistan"    = "#3CAEA3",
        "Kenya"       = "#F6D55C",
        "Zambia"      = "#ed553b")
      
      ## 1) GDM ---- 
      gam_gdm <- ggplot(new.data, aes(x = A1C, y = prob.gdm, color =Site, group = Site)) +
        facet_wrap(~ Site, ncol = 1) +
        geom_line(linewidth = 0.7) +
        
        geom_ribbon(aes(ymin = lower.gdm, ymax = upper.gdm), alpha = 0.15, color = NA) +
        scale_color_manual(values = site_cols) +
        scale_fill_manual(values = site_cols) +
        scale_x_continuous(
          breaks = seq(3.4, 6.4, by = 0.2),    
          limits = c(3.4, 6.4)) + 
        
        scale_y_continuous(
          breaks = seq(0, 0.5, by = 0.1),
          labels = percent_format(accuracy = 1),
          limits = c(0, 0.45)) +
        
        labs(title = "HbA1c and GDM",
             x = "Early HbA1c (%)",
             y = "Predicted probability") + 
        
        geom_vline(xintercept = 5.7, linetype = "dashed", color = "gray40") +
        geom_vline(xintercept = 5.4, linetype = "dashed", color = "gray40") +
        
        annotate("text", x = 5.7, y = 0.02, label = "5.7%", hjust = -0.05,  size = 3.5) + 
        annotate("text", x = 5.4, y = 0.02, label = "5.4%", hjust = -0.05,  size = 3.5) + 
        
        theme_classic() +
        theme(strip.background = element_blank(),
              legend.position = "none") 
   
      ## 2) HDP ---- 
      gam_hdp <- ggplot(new.data, aes(x = A1C, y = prob.hdp, color =Site, group = Site)) +
        facet_wrap(~ Site, ncol = 1) +
        geom_line(linewidth = 0.7) +
        
        geom_ribbon(aes(ymin = lower.hdp, ymax = upper.hdp), alpha = 0.15, color = NA) +
        scale_color_manual(values = site_cols) +
        scale_fill_manual(values = site_cols) +
        scale_x_continuous(
          breaks = seq(3.4, 6.4, by = 0.2),    
          limits = c(3.4, 6.4)) +
        
        scale_y_continuous(
          breaks = seq(0, 0.5, by = 0.1),
          labels = percent_format(accuracy = 1),
          limits = c(0, 0.45)) +
        
        labs(title = "HbA1c and HDP",
             x = "Early HbA1c (%)") + 
         
        geom_vline(xintercept = 5.7, linetype = "dashed", color = "gray40") +
        geom_vline(xintercept = 5.4, linetype = "dashed", color = "gray40") +
        
        annotate("text", x = 5.7, y = 0.02, label = "5.7%", hjust = -0.05,  size = 3.5) + 
        annotate("text", x = 5.4, y = 0.02, label = "5.4%", hjust = -0.05,  size = 3.5) + 
        
        theme_classic() +
        theme(strip.background = element_blank(),
              legend.position = "none",
              axis.title.y = element_blank())
      
      ## 3) Unplanned CS ---- 
      gam_cs <- ggplot(new.data, aes(x = A1C, y = prob.cs, color =Site, group = Site)) +
        facet_wrap(~ Site, ncol = 1) +
        geom_line(linewidth = 0.7) +
        
        geom_ribbon(aes(ymin = lower.cs, ymax = upper.cs), alpha = 0.15, color = NA) +
        scale_color_manual(values = site_cols) +
        scale_fill_manual(values = site_cols) +
        scale_x_continuous(
          breaks = seq(3.4, 6.4, by = 0.2),    
          limits = c(3.4, 6.4)) +
        
        scale_y_continuous(
          breaks = seq(0, 0.5, by = 0.1),
          labels = percent_format(accuracy = 1),
          limits = c(0, 0.45)) +
        labs(title = "HbA1c and emergent C/S",
             x = "Early HbA1c (%)") +
        
        geom_vline(xintercept = 5.7, linetype = "dashed", color = "gray40") +
        geom_vline(xintercept = 5.4, linetype = "dashed", color = "gray40") +
        
        annotate("text", x = 5.7, y = 0.02, label = "5.7%", hjust = -0.05,  size = 3.5) + 
        annotate("text", x = 5.4, y = 0.02, label = "5.4%", hjust = -0.05,  size = 3.5) + 
        
        theme_classic() +
        theme(strip.background = element_blank(),
              legend.position = "none",
              axis.title.y = element_blank()) 

      ## 4) LGA ---- 
      gam_lga <- ggplot(new.data, aes(x = A1C, y = prob.lga, color =Site, group = Site)) +
        facet_wrap(~ Site, ncol = 1) +
        geom_line(linewidth = 0.7) +
        
        geom_ribbon(aes(ymin = lower.lga, ymax = upper.lga), alpha = 0.15, color = NA) +
        scale_color_manual(values = site_cols) +
        scale_fill_manual(values = site_cols) +
        scale_x_continuous(
          breaks = seq(3.4, 6.4, by = 0.2),    
          limits = c(3.4, 6.4)) +
        
        scale_y_continuous(
          breaks = seq(0, 0.2, by = 0.1),
          labels = percent_format(accuracy = 1),
          limits = c(0, 0.2)) +
        
        labs(title = "HbA1c and LGA",
             x = "Early HbA1c (%)") +
          
        geom_vline(xintercept = 5.7, linetype = "dashed", color = "gray40") +
        geom_vline(xintercept = 5.4, linetype = "dashed", color = "gray40") +
        
        annotate("text", x = 5.7, y = 0.02, label = "5.7%", hjust = -0.05,  size = 3.0) + 
        annotate("text", x = 5.4, y = 0.02, label = "5.4%", hjust = -0.05,  size = 3.0) + 
        
        theme_classic() +
        theme(strip.background = element_blank(),
              legend.position = "none",
              axis.title.y = element_blank())
      
      # Create final GAM figure  
      gam_fig <- wrap_plots(gam_gdm, gam_hdp, gam_cs, gam_lga, ncol = 4)
      
#*******************************************************************************
#                XI. ROC analysis among overall population ----
#*******************************************************************************
    # GLM model that includes only A1C in the model
    dt_m1 <- d %>% drop_na(A1C)  # remove missingness 
    model1 <- glm(GDM ~ A1C, data = dt_m1, family = binomial (link = "logit")) 
  
    # Create predicted probabilities and ROC  
    dt_m1$prob_m1 <- predict(model1, type = "response") # predicted probability of GDM   
    range(dt_m1$prob_m1)                                # Nothing greater than 1 as the predicted probabilities can only be 0-1 (0% to 100%)
      
    roc_m1 <- roc(dt_m1$GDM, dt_m1$prob_m1)             # ROC analysis: For each observation, dt_m1$GDM = actual GDM status, dt_m1$prob_m1 = predicted probabilitie of GDM 
    
    length(unique(dt_m1$prob_m1))                       # the N of unique predicted probability (or threshold)
    length(unique(dt_m1$A1C)) 
    round(auc(roc_m1), 3)                               # get AUC
  
      ## Use Delong's method to calculate AUC's CI  
       ci.auc(roc_m1, conf.level = 0.95, method = "delong")
  
    # create sen, spec and threshold table (NOTE : unique predicted probabilities = unique threshold)
    roc_df_all <- data.frame(
       threshold = roc_m1$thresholds,
       sensitivity = roc_m1$sensitivities,
       specificity = roc_m1$specificities)
    
    # create a variable youden (Youden's index) for each threshold 
    roc_df_all$youden <- roc_df_all$sensitivity + roc_df_all$specificity - 1
    
    # create a variable distance to indicate the distance between each threshold to leftupper (0, 1) in ROC (Euclidean distance)
    roc_df_all$distance <- sqrt((1 - roc_df_all$sensitivity)^2 + (roc_df_all$specificity - 1)^2)
    
    # add corresponding A1C values that generate those predicted probabilities
          # Use the inverse of the regression model - provide predicted probabilities, b0 and b1 to calcualte back to HbA1c  
          intercept <- coef(model1)[1] # intercept 
          beta <- coef(model1)[2]      # coefficient for A1C
          summary(model1)
  
          roc_df_all$corres_a1c<- (log(roc_df_all$threshold / (1 - roc_df_all$threshold)) - intercept) / beta
  
          # extract the unique A1C
          unique_a1c_prob <- dt_m1 %>%
            distinct(A1C, prob_m1) %>%
            arrange(A1C)
          
          unique_a1c_prob <- dplyr::rename(
            unique_a1c_prob, threshold = prob_m1) 
      
    # add true positive 
    roc_df_all$TP <- sapply(roc_df_all$threshold, function(thresh) {
      sum(dt_m1$prob_m1 >= thresh & dt_m1$GDM == 1, na.rm = TRUE)})
    
    # add false positive 
    roc_df_all$FP <- sapply(roc_df_all$threshold, function(thresh) {
      sum(dt_m1$prob_m1 >= thresh & dt_m1$GDM == 0, na.rm = TRUE)})
    
    # add true negative 
    roc_df_all$TN <- sapply(roc_df_all$threshold, function(thresh) {
      sum(dt_m1$prob_m1 < thresh & dt_m1$GDM == 0, na.rm = TRUE)})
    
    # add false negative 
    roc_df_all$FN <- sapply(roc_df_all$threshold, function(thresh) {
      sum(dt_m1$prob_m1 < thresh & dt_m1$GDM == 1, na.rm = TRUE)})
    
    # sum up -> should be the entire sample 
    roc_df_all$Total <- with(roc_df_all, TP + FP + TN + FN) 
       all(roc_df_all$Total == nrow(dt_m1)) # check - TRUE! 
    
    # false discovery rate = FP / (TP + FP)
    roc_df_all$FDR <- with(roc_df_all, ifelse(TP + FP > 0, FP / (TP + FP), NA))
    
    # positive predictive value = TP / (TP + FP) = positive predictive value 
    roc_df_all$PPV <-  roc_df_all$TP / (roc_df_all$TP + roc_df_all$FP)
    
    # negative predictive value = TN / (TN + FN) = positive predictive value 
    roc_df_all$NPV <-  roc_df_all$TN / (roc_df_all$TN + roc_df_all$FN)
    
    # accuracy 
    roc_df_all$accuracy <- (roc_df_all$TP + roc_df_all$TN )/ roc_df_all$Total 
    
    # retrieve HbA1c, sen, spe for a specific Youden and smallest distance 
    roc_df_all[which.max(roc_df_all$youden), ]
    roc_df_all[which.min(roc_df_all$distance), ]
    
    print(roc_df_all) 

          # extract the a1c that has the largest youden 
          roc_df_all$corres_a1c[which.max(roc_df_all$youden)]
          roc_df_all[which.max(roc_df_all$youden), ]
          
          # extract the a1c that has the smallest Euclidean 
          roc_df_all$corres_a1c[which.min(roc_df_all$distance)]
          roc_df_all[which.min(roc_df_all$distance), ]
          
          # extract the a1c that has the sen that is closest to 80% - rule out threshold 
          roc_df_all$corres_a1c[which.min(abs(roc_df_all$sensitivity - 0.8))]
          roc_df_all[which.min(abs(roc_df_all$sensitivity - 0.8)), ]
          
          # extract the a1c that has the spec that is closest to 80% - rule in threshold 
          roc_df_all$corres_a1c[which.min(abs(roc_df_all$specificity - 0.8))]
          roc_df_all[which.min(abs(roc_df_all$specificity - 0.8)), ] 

#*******************************************************************************
#                         XII. ROC analysis by site ---- 
#*******************************************************************************
  # create site list 
  site_list <- c("North India", "South India", "Pakistan", "Kenya", "Zambia")
  
  # Create an empty list to store all ROC data for each site
  roc_results_list <- list()
  auc_list <- c()    
  roc_obj_list <-list ()
  
  # create a for loop
  for (i in seq_along(site_list)) { 
    site_name <- site_list[[i]]
    
    # Create site specific dataset 
    dt_site <- d %>% filter(Site == site_name) %>% drop_na(A1C)
    
    # Run log binomial model
    model <- glm(GDM ~ A1C, data = dt_site, family = binomial(link = "logit"))
    
    # create predicted probability 
    dt_site$prob <- predict(model, type = "response")
    
    # Use ROC and create ROC object
    roc_obj <- roc(dt_site$GDM, dt_site$prob)
    
    # Store ROC result for each site
    roc_obj_list[[site_name]] <- roc_obj  
    
    # Create ROC data frame  
    roc_df_site <- data.frame(threshold   = roc_obj$thresholds,
                              sensitivity = roc_obj$sensitivities,
                              specificity = roc_obj$specificities) 
    
    # Add Youden's index and Euclidean distance 
    roc_df_site$youden <- roc_df_site$sensitivity + roc_df_site$specificity - 1
    roc_df_site$distance <- sqrt((1 - roc_df_site$sensitivity)^2 + (roc_df_site$specificity - 1)^2)
    
    # Obtain corresponding HbA1C 
         # Use inverse of the regression model - provide predicted probabilities, b0 and b1 to calculate back to HbA1c  
         intercept <- coef(model)[1] # intercept from the model 
         beta <- coef(model)[2]      # coefficient for A1C from the model 
         
    roc_df_site$corres_a1c <- (log(roc_df_site$threshold / (1 - roc_df_site$threshold)) - intercept) / beta 
    
    # Obtain TP FP TN FN 
    roc_df_site$TP <- sapply(roc_df_site$threshold, function(thresh) {
      sum(dt_site$prob >= thresh & dt_site$GDM == 1)})
    
    roc_df_site$FP <- sapply(roc_df_site$threshold, function(thresh) {
      sum(dt_site$prob >= thresh & dt_site$GDM == 0)})
    
    roc_df_site$TN <- sapply(roc_df_site$threshold, function(thresh) {
      sum(dt_site$prob < thresh & dt_site$GDM == 0)})
    
    roc_df_site$FN <- sapply(roc_df_site$threshold, function(thresh) {
      sum(dt_site$prob < thresh & dt_site$GDM == 1)})
    
    # Obtain total, fdr, PPV, NPV, accuracy 
    roc_df_site$Total <- with(roc_df_site, TP + FP + TN + FN)
    
    roc_df_site$FDR <- with(roc_df_site, ifelse(TP + FP > 0, FP / (TP + FP), NA))
    
    roc_df_site$PPV <- with(roc_df_site, ifelse(TP + FP > 0, TP / (TP + FP), NA))
     
    roc_df_site$NPV <- with(roc_df_site, ifelse(TN + FN > 0, TN / (TN + FN), NA))
    
    roc_df_site$accuracy <- (roc_df_site$TP + roc_df_site$TN )/ roc_df_site$Total 
    
    roc_results_list[[site_name]] <- roc_df_site }  # Store results in list 
    
   ## retrieve HbA1c based on youden and Euclidean distance for each site
  
    # India SAS (North India)
      # extract the a1c that has the largest youden 
      roc_results_list[[1]]$corres_a1c[which.max(roc_results_list[[1]]$youden)]
      roc_results_list[[1]][which.max(roc_results_list[[1]]$youden), ]
      
      # extract the a1c that has the smallest Euclidean 
      roc_results_list[[1]]$corres_a1c[which.min(roc_results_list[[1]]$distance)]
      roc_results_list[[1]][which.min(roc_results_list[[1]]$distance), ]
      
      # extract the a1c that has the sen that is closest to 80% 
      roc_results_list[[1]]$corres_a1c[which.min(abs(roc_results_list[[1]]$sensitivity - 0.8))]
      roc_results_list[[1]][which.min(abs(roc_results_list[[1]]$sensitivity - 0.8)), ]
      
      # extract the a1c that has the spec that is closest to 80% 
      roc_results_list[[1]]$corres_a1c[which.min(abs(roc_results_list[[1]]$specificity - 0.8))]
      roc_results_list[[1]][which.min(abs(roc_results_list[[1]]$specificity - 0.8)), ]
  
    # India CMC (South India)
      # extract the a1c that has the largest youden 
      roc_results_list[[2]]$corres_a1c[which.max(roc_results_list[[2]]$youden)]
      roc_results_list[[2]][which.max(roc_results_list[[2]]$youden), ]
      
      # extract the a1c that has the smallest Euclidean 
      roc_results_list[[2]]$corres_a1c[which.min(roc_results_list[[2]]$distance)]
      roc_results_list[[2]][which.min(roc_results_list[[2]]$distance), ]
      
      # extract the a1c that has the sen that is closest to 80% 
      roc_results_list[[2]]$corres_a1c[which.min(abs(roc_results_list[[2]]$sensitivity - 0.8))]
      roc_results_list[[2]][which.min(abs(roc_results_list[[2]]$sensitivity - 0.8)), ]
      
      # extract the a1c that has the spec that is closest to 80% 
      roc_results_list[[2]]$corres_a1c[which.min(abs(roc_results_list[[2]]$specificity - 0.8))]
      roc_results_list[[2]][which.min(abs(roc_results_list[[2]]$specificity - 0.8)), ]
    
    # Pakistan  
      # extract the a1c that has the largest youden 
      roc_results_list[[3]]$corres_a1c[which.max(roc_results_list[[3]]$youden)]
      roc_results_list[[3]][which.max(roc_results_list[[3]]$youden), ]
      
      # extract the a1c that has the smallest Euclidean 
      roc_results_list[[3]]$corres_a1c[which.min(roc_results_list[[3]]$distance)]
      roc_results_list[[3]][which.min(roc_results_list[[3]]$distance), ]
      
      # extract the a1c that has the sen that is closest to 80% 
      roc_results_list[[3]]$corres_a1c[which.min(abs(roc_results_list[[3]]$sensitivity - 0.8))]
      roc_results_list[[3]][which.min(abs(roc_results_list[[3]]$sensitivity - 0.8)), ]
      
      # extract the a1c that has the spec that is closest to 80% 
      roc_results_list[[3]]$corres_a1c[which.min(abs(roc_results_list[[3]]$specificity - 0.8))]
      roc_results_list[[3]][which.min(abs(roc_results_list[[3]]$specificity - 0.8)), ]
      
    # Kenya 
      # extract the a1c that has the largest youden 
      roc_results_list[[4]]$corres_a1c[which.max(roc_results_list[[4]]$youden)]
      roc_results_list[[4]][which.max(roc_results_list[[4]]$youden), ]
      
      # extract the a1c that has the smallest Euclidean 
      roc_results_list[[4]]$corres_a1c[which.min(roc_results_list[[4]]$distance)]
      roc_results_list[[4]][which.min(roc_results_list[[4]]$distance), ]
      
      # extract the a1c that has the sen that is closest to 80% 
      roc_results_list[[4]]$corres_a1c[which.min(abs(roc_results_list[[4]]$sensitivity - 0.8))]
      roc_results_list[[4]][which.min(abs(roc_results_list[[4]]$sensitivity - 0.8)), ]
      
      # extract the a1c that has the spec that is closest to 80% 
      roc_results_list[[4]]$corres_a1c[which.min(abs(roc_results_list[[4]]$specificity - 0.8))]
      roc_results_list[[4]][which.min(abs(roc_results_list[[4]]$specificity - 0.8)), ]
    
    # Zambia  
      # extract the a1c that has the largest youden 
      roc_results_list[[5]]$corres_a1c[which.max(roc_results_list[[5]]$youden)]
      roc_results_list[[5]][which.max(roc_results_list[[5]]$youden), ]
      
      # extract the a1c that has the smallest Euclidean 
      roc_results_list[[5]]$corres_a1c[which.min(roc_results_list[[5]]$distance)]
      roc_results_list[[5]][which.min(roc_results_list[[5]]$distance), ]
      
      # extract the a1c that has the sen that is closest to 80% 
      roc_results_list[[5]]$corres_a1c[which.min(abs(roc_results_list[[5]]$sensitivity - 0.8))]
      roc_results_list[[5]][which.min(abs(roc_results_list[[5]]$sensitivity - 0.8)), ]
      
      # extract the a1c that has the spec that is closest to 80% 
      roc_results_list[[5]]$corres_a1c[which.min(abs(roc_results_list[[5]]$specificity - 0.8))]
      roc_results_list[[5]][which.min(abs(roc_results_list[[5]]$specificity - 0.8)), ]
      
    write_xlsx(roc_results_list, "roc by site table.xlsx") 
  
    ### 95% CI for AUC by site using DeLong's method
      ci.auc(roc_obj_list[["Kenya"]], conf.level = 0.95, method = "delong")
      ci.auc(roc_obj_list[["Zambia"]], conf.level = 0.95, method = "delong")
      ci.auc(roc_obj_list[["South India"]], conf.level = 0.95, method = "delong")
      ci.auc(roc_obj_list[["North India"]], conf.level = 0.95, method = "delong")
      ci.auc(roc_obj_list[["Pakistan"]], conf.level = 0.95, method = "delong")

    
#*******************************************************************************
#                             XIII. AUC figure ----
#*******************************************************************************

  # build a named list of ROC objects (names become the legend)
  roc_list <- list(
    "All sites"     = roc_m1,
    "North India"   = roc_obj_list[[1]],
    "South India"   = roc_obj_list[[2]],
    "Pakistan"      = roc_obj_list[[3]],
    "Kenya"         = roc_obj_list[[4]],
    "Zambia"        = roc_obj_list[[5]])
    
  auc_ci_list <- list(
      "All sites"     = ci.auc(roc_m1, conf.level = 0.95, method = "delong"),
      "North India"   = ci.auc( roc_obj_list[[1]], conf.level = 0.95, method = "delong"),
      "South India"   = ci.auc( roc_obj_list[[2]], conf.level = 0.95, method = "delong"),
      "Pakistan"      = ci.auc( roc_obj_list[[3]], conf.level = 0.95, method = "delong"),
      "Kenya"         = ci.auc( roc_obj_list[[4]], conf.level = 0.95, method = "delong"),
      "Zambia"        = ci.auc( roc_obj_list[[5]], conf.level = 0.95, method = "delong"))
    

  # make legend labels that include AUC
    auc_labels <- sapply(names(roc_list), function(x) {
                  paste0(x, " AUC=", 
                         round(as.numeric(auc(roc_list[[x]])), 3), " (95% CI ", 
                         round(ci.auc(roc_list[[x]],  conf.level = 0.95, method = "delong")[1], 3), "-",
                         round(ci.auc(roc_list[[x]],  conf.level = 0.95, method = "delong")[3], 3), ")")    
    })
    names(auc_labels) <- names(roc_list) 

  # create ROC figure 
    auc.figure <- ggroc(roc_list, legacy.axes = TRUE, size = 1.1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      labs(
        #title = "HbA1c in predicting GDM - ROC and AUC comparison",
        x = "1 - Specificity",
        y = "Sensitivity",
        color = NULL
      ) +
      scale_color_discrete(labels = auc_labels) +
      
      scale_x_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0), limits = c(0,1)) +
      scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.0), limits = c(0,1)) +
      
      theme_classic(base_size = 12.5) +
      theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
    
    ggsave("AUC figure.tif", plot = auc.figure, width = 10, height = 9, units = "in", dpi = 600)
    
#*******************************************************************************
#                         XIV. Site specific estimates ----
#*******************************************************************************

  # create site list
  Sites <- c("South India", "North India", "Pakistan", "Kenya", "Zambia") 
  
  ### 1) GDM ----
    # Model 1
    GDM_result_1 <- list()      # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(GDM ~ a1c, 
             data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) 
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      GDM_result_1[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    GDM_result_model1 <- do.call(rbind, GDM_result_1)  
    
    # Model 2 
    GDM_result_2 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {         # Within a specific outcome, a loop for each site 
      
      fit <- glm(GDM ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) 
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
      
      GDM_result_2[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    GDM_result_model2 <- do.call(rbind, GDM_result_2)  
    
    # Model 3 
    GDM_result_3 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(GDM~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM,
                 data = subset(d_1, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) 
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE)   
      
      GDM_result_3[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    GDM_result_model3 <- do.call(rbind, GDM_result_3)  
    GDM_result<-rbind(GDM_result_model1, GDM_result_model2, GDM_result_model3)
  
  ### 2) HDP ----
    # Model 1
    HDP_result_1 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {        # Within a specific outcome, a loop for each site 
      
      fit <- glm(HDP ~ a1c, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) 
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
      
      HDP_result_1[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
      HDP_result_model1 <- do.call(rbind, HDP_result_1)  
    
    # Model 2 
    HDP_result_2 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {        # Within a specific outcome, a loop for each site 
      
      fit <- glm(HDP ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) 
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE)  
      
      HDP_result_2[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    HDP_result_model2 <- do.call(rbind, HDP_result_2)  
    
    # Model 3 
    HDP_result_3 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {         # Within a specific outcome, a loop for each site 
      
      fit <- glm(HDP~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM,
                 data = subset(d_1, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      HDP_result_3[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    HDP_result_model3 <- do.call(rbind, HDP_result_3)  
    HDP_result<-rbind(HDP_result_model1, HDP_result_model2, HDP_result_model3)
    
    ### 3) emergent CS ----
    # Model 1 
    EMERGENT_CS_result_1 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {                 # Within a specific outcome, a loop for each site 
      
      fit <- glm(EMERGENT_CS ~ a1c, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) 
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE)
      
      EMERGENT_CS_result_1[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    EMERGENT_CS_result_model1 <- do.call(rbind, EMERGENT_CS_result_1)  
    
    # Model 2 
    EMERGENT_CS_result_2 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {                 # Within a specific outcome, a loop for each site 
      
      fit <- glm(EMERGENT_CS ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      EMERGENT_CS_result_2[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    EMERGENT_CS_result_model2 <- do.call(rbind, EMERGENT_CS_result_2)  
    
    # Model 3 
    EMERGENT_CS_result_3 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {                 # Within a specific outcome, a loop for each site 
      
      fit <- glm(EMERGENT_CS~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM,
                 data = subset(d_1, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      EMERGENT_CS_result_3[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    EMERGENT_CS_result_model3 <- do.call(rbind, EMERGENT_CS_result_3)  
    EMERGENT_CS_result<-rbind(EMERGENT_CS_result_model1, EMERGENT_CS_result_model2, EMERGENT_CS_result_model3)
    
    ### 4) Preterm birth ----
    # Model 1 
    PRETERM_result_1 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {            # Within a specific outcome, a loop for each site 
      
      fit <- glm(PRETERM ~ a1c, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      PRETERM_result_1[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    PRETERM_result_model1 <- do.call(rbind, PRETERM_result_1)  
    
    # Model 2 
    PRETERM_result_2 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(PRETERM ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      PRETERM_result_2[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    PRETERM_result_model2 <- do.call(rbind, PRETERM_result_2)  
    
    # Model 3 
    PRETERM_result_3 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(PRETERM ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM,
                 data = subset(d_1, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      PRETERM_result_3[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    PRETERM_result_model3 <- do.call(rbind, PRETERM_result_3)  
    PRETERM_result<-rbind(PRETERM_result_model1, PRETERM_result_model2, PRETERM_result_model3)
    
    ### 5) LGA ----
    # Model 1 
    LGA_result_1 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(LGA ~ a1c, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      LGA_result_1[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    LGA_result_model1 <- do.call(rbind, LGA_result_1)  
    
    # Model 2 
    LGA_result_2 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(LGA ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      LGA_result_2[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    LGA_result_model2 <- do.call(rbind, LGA_result_2)  
    
    # Model 3 
    LGA_result_3 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(LGA~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM,
                 data = subset(d_1, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      LGA_result_3[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    LGA_result_model3 <- do.call(rbind, LGA_result_3)  
    LGA_result<-rbind(LGA_result_model1, LGA_result_model2, LGA_result_model3)
    
    ### 6) STILLBIRTH ----
    # Model 1 
    STILLBIRTH_result_1 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(STILLBIRTH ~ a1c, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      STILLBIRTH_result_1[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    STILLBIRTH_result_model1 <- do.call(rbind, STILLBIRTH_result_1)  
    
    # Model 2 
    STILLBIRTH_result_2 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(STILLBIRTH ~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + MULTIPARITY + HBA1C_GA_WKS, 
                 data = subset(dt, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      STILLBIRTH_result_2[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 ( " (", LCI, "-", UCI, ")"),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    STILLBIRTH_result_model2 <- do.call(rbind, STILLBIRTH_result_2)  
    
    # Model 3 
    STILLBIRTH_result_3 <- list()        # To store site-specific estimate for a specific outcome  
    
    for (site in Sites) {       # Within a specific outcome, a loop for each site 
      
      fit <- glm(STILLBIRTH~ a1c + AGE_GROUP + BMI4CAT + SCHOOL_YRS + HBA1C_GA_WKS + PREVPREG_GDM,
                 data = subset(d_1, Site == site), family = poisson(link = "log"))
      
      logRR <- coef(fit)["a1c"]
      SE    <- sqrt(sandwich::vcovHC(fit, type = "HC0")["a1c","a1c"]) # robust variances for modified Poisson
      RR    <- round(exp(logRR), 2)
      LCI   <- round(exp(logRR - 1.96 * SE), 2)
      UCI   <- round(exp(logRR + 1.96 * SE), 2)
      N_obs <- nobs(fit)
      zval  <- logRR / SE
      pval  <- 2 * pnorm(abs(zval), lower.tail = FALSE) # need to calculate p value  
      
      STILLBIRTH_result_3[[site]] <- data.frame(
        site  = site,
        N_obs = N_obs,
        RR    = RR,
        CI    = paste0 (LCI, "-", UCI),
        p     = ifelse(
          is.na(pval),
          NA_character_,
          ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))))}
    
    STILLBIRTH_result_model3 <- do.call(rbind, STILLBIRTH_result_3)  
    STILLBIRTH_result<-rbind(STILLBIRTH_result_model1, STILLBIRTH_result_model2, STILLBIRTH_result_model3)
