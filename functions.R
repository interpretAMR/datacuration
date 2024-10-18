# Author: Kat Holt (kat.holt@lshtm.ac.uk) & Zoe Dyson (zoe.dyson@lshtm.ac.uk)
# Title: functions.R
# Date: 18/10/2024


# Package management
# List required packages
packages <- c("tidyverse", "patchwork", "logistf")

# Load packages & install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))



# solo_marker: Summarise resistance markers by AST phenotype data
#
#   Parameters: 
#             1. ast: A data frame containing AST data - must have columns 
#                'Antibiotic', `Resistance phenotype`, sample IDs in 'BioSample'
#             2. amr_binary: A data frame containing binary presence/absence data 
#                 for markers, sample IDs in 'sample_accession' from Enterobase
#             3. gene_class_list: A data frame containing AMR determinant data -
#                 must have columns 'gene', 'Class' and 'Subclass'
#             4. refgene_class OR refgene_subclass: A string for the class or 
#                 subclass for which to identify 'solo' markers
#
#   Returns: A list containing
#             1. solo_stats:
#.            2. solo_count: 
#             3. solo_count_plot: 
#             4. combined_plot:
#
solo_marker <- function(ast, antibiotic, gene_class_list, amr_binary, 
                        refgene_class=NULL, refgene_subclass=NULL, 
                        plot_cols=c("resistant"="IndianRed", 
                                    "intermediate"="orange", 
                                    "nonsusceptible"="orange", 
                                    "susceptible"="lightgrey", 
                                    "not defined"="white")) {
  
  # extract AST data for this drug
  ast <- ast %>%
    mutate(`Resistance phenotype`=if_else(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible"), `Resistance phenotype`, "not defined"))
  ast_drug <- ast %>% filter(Antibiotic==antibiotic)
  
  # extract list of relevant drug markers
  if (!is.null(refgene_class)) {
    genes <- afp_matched %>% 
      filter(Class %in% refgene_class) %>% 
      pull(`Gene symbol`)
    
  }
  else if (!is.null(refgene_subclass)) {
    genes <- afp_matched %>% 
      filter(grepl(refgene_subclass[1], Subclass)) %>% 
      pull(`Gene symbol`)
    
    if (length(refgene_subclass)>1) {
      for (i in 2:length(refgene_subclass)) {
        genes <- c(genes, afp_matched %>% 
                     filter(grepl(refgene_subclass[i], Subclass)) %>% 
                     pull(`Gene symbol`))
        
      }
    }
  }
  else {stop("must specify refgene_class or refgene_subclass")}
  
  
  # create longform binary of relevant drug markers from input binary matrix
  calls_long <- amr_binary %>% 
    select(sample_accession,any_of(genes)) %>% 
    pivot_longer(!sample_accession, names_to="gene")
  
  # marker count vs phenotype
  ast_drug_genes <- calls_long %>% 
    group_by(sample_accession) %>% 
    summarise(count=sum(value)) %>% 
    ungroup() %>%
    right_join(ast_drug, join_by("sample_accession"=="BioSample")) %>% 
    relocate(count, .after=sample_accession)
  
  # solo genes
  solo <- ast_drug_genes %>% 
    filter(count==1) %>% 
    left_join(calls_long) %>% 
    filter(value==1)
  
  solo_count <- solo %>% 
    group_by(gene) %>% 
    count()
  
  solo_count_plot <- solo_count %>% 
    ggplot(aes(x=gene, y=n)) + 
    geom_col() + 
    coord_flip() + 
    theme_light()
  
  # solo PPV
  solo_stats_R <- solo %>% 
    group_by(gene) %>% 
    summarise(total=sum(`Resistance phenotype` %in% c("resistant", 
                                                      "intermediate", 
                                                      "susceptible")), 
              n=sum(`Resistance phenotype`=="resistant"), 
              p=n/total, se=sqrt(p*(1-p)/total), 
              ci.lower=max(0,p-1.96*se), 
              ci.upper=min(1,p+1.96*se)) %>% 
    mutate(category="resistant")
  
  solo_stats_NWT <- solo %>% 
    group_by(gene) %>% 
    summarise(total=sum(`Resistance phenotype` %in% c("resistant", 
                                                      "intermediate", 
                                                      "susceptible")), 
              n=sum(`Resistance phenotype`%in% c("resistant", "intermediate")), 
              p=n/total, se=sqrt(p*(1-p)/total), 
              ci.lower=max(0,p-1.96*se), 
              ci.upper=min(1,p+1.96*se)) %>% 
    mutate(category="nonsusceptible")
  
  solo_stats <- bind_rows(solo_stats_R, solo_stats_NWT) %>% 
    relocate(category, .before=n)
  
  # solo plots
  
  pd <- position_dodge(width=0.8) # position dodge for coefficient plots
  
  ppv_plot <- solo_stats %>% 
    ggplot(aes(y=gene, group=category, col=category)) +
    geom_vline(xintercept=0.5, linetype=2) +
    geom_linerange(aes(xmin=ci.lower, xmax=ci.upper), position=pd) +
    geom_point(aes(x=p), position=pd) + 
    theme_bw() +
    scale_y_discrete(labels=paste0("(n=",solo_count$n,")"), position="right") + 
    labs(y="", x="PPV", col="Category") + 
    scale_colour_manual(values=plot_cols) + xlim(0,1)
  
  solo_pheno_plot <- solo %>% 
    mutate(`Resistance phenotype`=fct_relevel(`Resistance phenotype`, 
                                              "not defined", 
                                              "susceptible", 
                                              "intermediate", 
                                              "resistant")) %>%
    ggplot(aes(x=gene, fill=`Resistance phenotype`)) + 
    geom_bar(stat="count", position="fill") + 
    scale_fill_manual(values=plot_cols) + 
    coord_flip() + 
    geom_text(aes(label=..count..), stat="count", position=position_fill(vjust = .5), size=3) + 
    theme_light() + labs(x="", y="Proportion", fill="Phenotype")
  
  if (!is.null(refgene_class)) {
    header=paste0("Solo markers for class '", tolower(refgene_class), "'") 
  }else if (!is.null(refgene_subclass)) {
    header=paste0("Solo markers for subclass '", tolower(refgene_subclass), "'") 
  }
  
  combined_plot <- solo_pheno_plot + 
    ppv_plot + 
    plot_layout(axes="collect", guides="collect") + 
    plot_annotation(title=header, subtitle=paste("vs", antibiotic, "phenotype"))
  
  return(list(solo_stats=solo_stats, solo_count, solo_count_plot, combined_plot=combined_plot))
}


solo_marker_atb <- function(ast_matched, antibiotic, afp_matched, refgene_class=NULL, refgene_subclass=NULL, 
                            plot_cols=c("resistant"="IndianRed", 
                                        "intermediate"="orange", 
                                        "susceptible"="lightgrey", 
                                        "nonsusceptible"="orange")) {
  
  # extract list of relevant drug markers
  if (!is.na(refgene_class)) {
    
    genes <- afp_matched %>% 
      filter(Class %in% refgene_class) %>% 
      pull(`Gene symbol`)
  }
  else if (!is.na(refgene_subclass)) {
    genes <- afp_matched %>% 
      filter(grepl(refgene_subclass[1], Subclass)) %>% 
      pull(`Gene symbol`)
    
    if (length(refgene_subclass)>1) {
      for (i in 2:length(refgene_subclass)) {
        genes <- c(genes, afp_matched %>% 
                     filter(grepl(refgene_subclass[i], Subclass)) %>% 
                     pull(`Gene symbol`))
      }
    }
  }
  else {stop("must specify refgene_class or refgene_subclass")}
  
  afp_solo_strains <- afp_matched %>% filter(`Gene symbol` %in% genes) %>% 
    group_by(Name) %>% count() %>% filter(n==1) %>% pull(Name)
  
  # including all samples with AST for this drug and AFP results, even if no genes reported
  afp_solo_ast <- afp_matched %>% 
    filter(`Gene symbol` %in% genes) %>% 
    filter(Name %in% afp_solo_strains) %>% 
    right_join(ast_matched %>% 
                 filter(Antibiotic==antibiotic & !is.na(`Resistance phenotype`)), join_by("Name"=="BioSample"))
  
  # solo genes
  solo_count <- afp_solo_ast %>% 
    group_by(`Gene symbol`) %>% 
    count()
  
  solo_stats <- NULL
  solo_count_plot <- NULL
  combined_plot <- NULL
  
  if (nrow(solo_count)>0) {
    solo_count_plot <- solo_count %>% 
      filter(!is.na(`Gene symbol`)) %>%
      ggplot(aes(x=`Gene symbol`, y=n)) + 
      geom_col() + 
      coord_flip() + 
      theme_light()
    
    # solo PPV
    solo_stats_R <- afp_solo_ast %>% 
      group_by(`Gene symbol`) %>% 
      summarise(total=sum(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible")), 
                n=sum(`Resistance phenotype`=="resistant"), 
                p=n/total, 
                se=sqrt(p*(1-p)/total), 
                ci.lower=max(0,p-1.96*se), 
                ci.upper=min(1,p+1.96*se)) %>% 
      mutate(category="resistant")
    solo_stats_NWT <- afp_solo_ast %>% 
      group_by(`Gene symbol`) %>% 
      summarise(total=sum(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible")), 
                n=sum(`Resistance phenotype`%in% c("resistant", "intermediate")), 
                p=n/total, 
                se=sqrt(p*(1-p)/total), 
                ci.lower=max(0,p-1.96*se), 
                ci.upper=min(1,p+1.96*se)) %>% 
      mutate(category="nonsusceptible")
    solo_stats <- bind_rows(solo_stats_R, solo_stats_NWT) %>% 
      relocate(category, .before=n) %>%
      mutate(Antibiotic=antibiotic) %>%
      filter(!is.na(`Gene symbol`))
    
    # solo plots
    
    pd <- position_dodge(width=0.8) # position dodge for coefficient plots
    
    ppv_plot <- solo_stats %>% 
      filter(!is.na(`Gene symbol`)) %>%
      ggplot(aes(y=`Gene symbol`, group=category, col=category)) +
      geom_vline(xintercept=0.5, linetype=2) +
      geom_linerange(aes(xmin=ci.lower, xmax=ci.upper), position=pd) +
      geom_point(aes(x=p), position=pd) + 
      theme_bw() +
      scale_y_discrete(labels=paste0("(n=",solo_count$n,")"), position="right") + 
      labs(y="", x="PPV", col="Category") + 
      scale_colour_manual(values=plot_cols) + xlim(0,1)
    
    solo_pheno_plot <- afp_solo_ast %>% 
      filter(!is.na(`Gene symbol`)) %>%
      mutate(`Resistance phenotype`=fct_relevel(`Resistance phenotype`, "susceptible", "intermediate", "resistant")) %>%
      ggplot(aes(x=`Gene symbol`, fill=`Resistance phenotype`)) + 
      geom_bar(stat="count", position="fill") + 
      scale_fill_manual(values=plot_cols) + 
      coord_flip() + 
      geom_text(aes(label=..count..), stat="count", position=position_fill(vjust = .5), size=3) + 
      theme_light() + labs(x="", y="Proportion", fill="Phenotype")
    
    if (!is.na(refgene_class)) {
      header=paste0("Solo markers for class '", tolower(refgene_class), "'") 
    }else if (!is.na(refgene_subclass)) {
      header=paste0("Solo markers for subclass '", tolower(refgene_subclass), "'") 
    }
    
    combined_plot <- solo_pheno_plot + 
      ppv_plot + 
      plot_layout(axes="collect", guides="collect") + 
      plot_annotation(title=header, subtitle=paste("vs", antibiotic, "phenotype"))
    
  }
  
  return(list(solo_stats=solo_stats, solo_count=solo_count, solo_count_plot=solo_count_plot, combined_plot=combined_plot))
}

# getBinMatrix: 
#
#   Parameters: 
#             1. ast_matched:
#             2. antibiotic:
#             3. afp_matched:
#             4. refgene_class or refgene_subclass:
#             5. maf:
#             6. plot_cols: A list for a colour palette for use in plotting -
#                default is Resistant=Red, Intermediate/Non-susceptible = Orange,
#                Susceptible=Gray.
#
#   Returns: 
#             afp_binary_markers:  
#
getBinMatrix <- function(ast_matched, antibiotic, afp_matched, refgene_class=NULL, 
                         refgene_subclass=NULL, maf=10, 
                         plot_cols=c("resistant"="IndianRed", 
                                     "intermediate"="orange", 
                                     "nonsusceptible"="orange", 
                                     "susceptible"="lightgrey")) {
  
  # extract list of relevant drug markers
  if (!is.null(refgene_class)) {
    genes <- afp_matched %>% 
      filter(Class %in% refgene_class) %>% 
      pull(`Gene symbol`)
  }
  else if (!is.null(refgene_subclass)) {
    genes <- afp_matched %>% 
      filter(grepl(refgene_subclass[1], Subclass)) %>% 
      pull(`Gene symbol`)
    if (length(refgene_subclass)>1) {
      for (i in 2:length(refgene_subclass)) {
        genes <- c(genes, afp_matched %>% filter(grepl(refgene_subclass[i], Subclass)) %>% 
                     pull(`Gene symbol`))
      }
    }
  }
  else {stop("must specify refgene_class or refgene_subclass")}
  
  afp_binary_markers <- afp_matched %>% 
    filter(`Gene symbol` %in% genes) %>% 
    group_by(Name, `Gene symbol`) %>% count() %>%
    ungroup() %>%
    right_join(ast_matched %>% filter(Antibiotic==antibiotic) %>% select(BioSample), 
               join_by("Name"=="BioSample")) %>%
    distinct() %>%
    pivot_wider(id_cols=Name, names_from=`Gene symbol`, values_from=n, values_fill=0) 
  
  if ("NA" %in% colnames(afp_binary_markers)) {
    afp_binary_markers <- afp_binary_markers %>% 
      select(-c("NA"))
  }
  
  afp_binary_markers <- afp_binary_markers %>% 
    select(Name, sort(names(.))) %>%
    inner_join(ast_matched %>% filter(Antibiotic==antibiotic & !is.na(`Resistance phenotype`)) %>% 
                 select(BioSample, `Resistance phenotype`), 
               join_by("Name"=="BioSample")) %>%
    mutate(resistant=case_when(`Resistance phenotype`=="resistant" ~ 1,
                               `Resistance phenotype` %in% c("intermediate", "susceptible") ~ 0,
                               TRUE ~ NA)) %>%
    select(-`Resistance phenotype`, -Name) %>%
    filter(!is.na(resistant)) %>%
    select(resistant, where(~sum(.x)>=maf))
  
  return(afp_binary_markers)
  
}


# logreg_details: Tabulate estimate, 95% CI, p-value from logistf regression model
#
#   Parameters: 
#             model: A lostif regression model object
#
#   Returns: 
#             model_summary: a tibble containing 95% CI and p-values from model 
#
logreg_details <- function(model) {
  
  # Create summary tibble
  model_summary <- cbind(est=model$coefficients, 
                         ci.lower=model$ci.lower, 
                         ci.upper=model$ci.upper, 
                         pval=model$prob) %>%
    as_tibble(rownames="marker")
  
  # Return summary tibble
  return(model_summary)
}


# glm_details: Tabulate estimate, 95% CI, p-value from glm regression model
#
#   Parameters: 
#             model: A glm regression model object
#
#   Returns: 
#             model_summary: a tibble containing 95% CI and p-values from model 
#
glm_details <- function(model) {
  
  # get CI data
  ci <- confint(model) %>% 
    as_tibble(rownames="marker") %>%
    rename(ci.lower=`2.5 %`, ci.upper=`97.5 %`)
  
  # Create summary tibble
  model_summary <- summary(model)$coef %>% 
    as_tibble(rownames="marker") %>% 
    rename(est=Estimate, pval=`Pr(>|z|)`) %>% 
    select(marker, est, pval) %>% 
    left_join(ci)
  
  # Return summary tibble
  return(model_summary)
}


# getPPV: Calculate Positive Predictive Value (PPV) for a list of drug/class pairs
#
#   Parameters: 
#             1. ast_matched: A data frame containing AST data matched to (2)
#             2. afp_matched: A data frame containing AMRFP data matched to (1)
#             3. test_list: A list of drugs/drug classes to test
#             4. outdir: An output directory for storing results  
#
#   Returns: 
#             1. solo_stats: A data frame of PPV data
#             2. creates pdf formatted figures for each test carried out in user
#                supplied output directory (no object returned for this)
#
getPPV <- function(ast_matched = ast_matched, afp_matched=afp_matched, test_list=NULL, outdir=NULL) {
  
  # Create empty data frame to store PPV data
  solo_stats <- NULL
  
  # Iterate thorugh drug list testing each drug
  for (i in 1:nrow(test_list)) {
    
    # Call solo_marker_atb()
    solo <- solo_marker_atb(ast_matched = ast_matched, 
                            afp_matched=afp_matched, 
                            antibiotic = test_list[i,1],
                            refgene_class = test_list[i,2], 
                            refgene_subclass = test_list[i,3])
    
    # extract subclass or class data (depending on which is present)
    class=test_list[i,2]
    
    if(is.na(class)) {
      class=test_list[i,3]
    }
    
    # Save figure to user supplied output directory
    ggsave(height=9,
           width=8, 
           solo$combined_plot, 
           file=paste0(outdir,"SoloPPV_",test_list[i,1],"_",class,".pdf"))
    
    # Add data for drug to data frame
    solo_stats <- bind_rows(solo_stats, solo$solo_stats)
  }
  
  # Return data frame of PPV data
  return(solo_stats)
}



getBinMatMIC <- function(ast_matched, antibiotic, afp_matched, refgene_class=NULL, refgene_subclass=NULL, plot_cols=c("resistant"="IndianRed", "intermediate"="orange", "nonsusceptible"="orange", "susceptible"="lightgrey")) {
  
  # extract list of relevant drug markers
  if (!is.null(refgene_class)) {
    genes <- afp_matched %>% filter(Class %in% refgene_class) %>% pull(`Gene symbol`)
  }
  else if (!is.null(refgene_subclass)) {
    genes <- afp_matched %>% filter(grepl(refgene_subclass[1], Subclass)) %>% pull(`Gene symbol`)
    if (length(refgene_subclass)>1) {
      for (i in 2:length(refgene_subclass)) {
        genes <- c(genes, afp_matched %>% filter(grepl(refgene_subclass[i], Subclass)) %>% pull(`Gene symbol`)) %>% unique()
      }
    }
  }
  else {stop("must specify refgene_class or refgene_subclass")}
  
  afp_binary_markers <- afp_matched %>% 
    filter(`Gene symbol` %in% genes) %>% 
    group_by(Name, `Gene symbol`) %>% count() %>%
    ungroup() %>%
    right_join(ast_matched %>% filter(Antibiotic==antibiotic) %>% select(BioSample), 
               join_by("Name"=="BioSample")) %>%
    distinct() %>%
    pivot_wider(id_cols=Name, names_from=`Gene symbol`, values_from=n, values_fill=0) 
  
  if ("NA" %in% colnames(afp_binary_markers)) {
    afp_binary_markers <- afp_binary_markers %>%select(-c("NA"))
  }
  
  afp_binary_markers <- afp_binary_markers %>% 
    select(Name, sort(names(.))) %>%
    inner_join(ast_matched %>% filter(Antibiotic==antibiotic) %>% select(BioSample, `Resistance phenotype`, `MIC (mg/L)`, `Disk diffusion (mm)`), 
               join_by("Name"=="BioSample")) %>%
    mutate(resistant=case_when(`Resistance phenotype`=="resistant" ~ 1,
                               `Resistance phenotype` %in% c("intermediate", "susceptible") ~ 0,
                               TRUE ~ NA)) %>%
    select(-`Resistance phenotype`, -Name)
  
  return(afp_binary_markers)
  
}
