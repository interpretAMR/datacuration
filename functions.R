# inputs
# ec_amr_meta: binary AMRfinder plus results, with metadata fields:
#         - `Clermont Type (ClermonTyping)`
#         - `Pathovar`
# marker: name of marker to summarise

summarise_marker_dist_ecoli <- function(ec_amr_meta=ec_amr_meta, marker) {
  clermont <- ec_amr_meta %>% group_by(`Clermont Type (ClermonTyping)`) %>% summarise(freq=mean(eval(parse(text=marker))))
  clermont_plot <- clermont %>% ggplot(aes(y=`Clermont Type (ClermonTyping)`, x=freq)) + geom_col() + theme_light() + ggtitle(marker)
  pathovar <- ec_amr_meta %>% group_by(Pathovar) %>% summarise(freq=mean(eval(parse(text=marker))))
  pathovar_plot <- pathovar %>% ggplot(aes(y=Pathovar, x=freq)) + geom_col() + theme_light() + ggtitle(marker)
  return(list(clermont=clermont, clermont_plot=clermont_plot, pathovar=pathovar, pathovar_plot=pathovar_plot))
}

# ast: must have columns 'Antibiotic', `Resistance phenotype`, sample IDs in 'BioSample'
# amr_binary: binary presence/absence data for markers, sample IDs in 'sample_accession'
# gene_class_list: must have columns 'gene', 'Class' and 'Subclass'
# regene_class OR regene_subclass: class or subclass for which to identify 'solo' markers
solo_marker <- function(ast, antibiotic, gene_class_list, amr_binary, refgene_class=NULL, refgene_subclass=NULL, plot_cols=c("resistant"="IndianRed", "intermediate"="orange", "nonsusceptible"="orange", "susceptible"="lightgrey", "not defined"="white")) {

  # extract AST data for this drug
  ast <- ast %>%
    mutate(`Resistance phenotype`=if_else(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible"), `Resistance phenotype`, "not defined"))
  ast_drug <- ast %>% filter(Antibiotic==antibiotic)
  
  # extract list of relevant drug markers
  if (!is.null(refgene_class)) {
    genes <- gene_class_list %>% filter(grepl(refgene_class, Class)) %>% pull(gene)
  }
  else if (!is.null(refgene_subclass)) {
    genes <- gene_class_list %>% filter(grepl(refgene_subclass, Subclass)) %>% pull(gene)
  }
  else {stop("must specify refgene_class or refgene_subclass")}
  
  # create longform binary of relevant drug markers from input binary matrix
  calls_long <- amr_binary %>% select(sample_accession,any_of(genes)) %>% pivot_longer(!sample_accession, names_to="gene")
  
  # marker count vs phenotype
  ast_drug_genes <- calls_long %>% group_by(sample_accession) %>% 
    summarise(count=sum(value)) %>% ungroup() %>%
    right_join(ast_drug, join_by("sample_accession"=="BioSample")) %>% 
    relocate(count, .after=sample_accession)
  
  # solo genes
  solo <- ast_drug_genes %>% filter(count==1) %>% left_join(calls_long) %>% filter(value==1)
  solo_count <- solo %>% group_by(gene) %>% count()
  solo_count_plot <- solo_count %>% ggplot(aes(x=gene, y=n)) + geom_col() + coord_flip() + theme_light()
  
  # solo PPV
  
  solo_stats_R <- solo %>% group_by(gene) %>% summarise(total=sum(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible")), n=sum(`Resistance phenotype`=="resistant"), p=n/total, se=sqrt(p*(1-p)/total), ci.lower=max(0,p-1.96*se), ci.upper=min(1,p+1.96*se)) %>% mutate(category="resistant")
  solo_stats_NWT <- solo %>% group_by(gene) %>% summarise(total=sum(`Resistance phenotype` %in% c("resistant", "intermediate", "susceptible")), n=sum(`Resistance phenotype`%in% c("resistant", "intermediate")), p=n/total, se=sqrt(p*(1-p)/total), ci.lower=max(0,p-1.96*se), ci.upper=min(1,p+1.96*se)) %>% mutate(category="nonsusceptible")
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
    mutate(`Resistance phenotype`=fct_relevel(`Resistance phenotype`, "not defined", "susceptible", "intermediate", "resistant")) %>%
    ggplot(aes(x=gene, fill=`Resistance phenotype`)) + 
    geom_bar(stat="count", position="fill") + 
    scale_fill_manual(values=plot_cols) + 
    coord_flip() + 
    geom_text(aes(label=..count..), stat="count", position=position_fill(vjust = .5), size=3) + 
    theme_light() + labs(x="", y="Proportion", fill="Phenotype")
  
  if (!is.null(refgene_class)) {header=paste0("Solo markers for class '", tolower(refgene_class), "'") }
  else if (!is.null(refgene_subclass)) {header=paste0("Solo markers for subclass '", tolower(refgene_subclass), "'") }
  combined_plot <- solo_pheno_plot + ppv_plot + plot_layout(axes="collect", guides="collect") + plot_annotation(title=header, subtitle=paste("vs", antibiotic, "phenotype"))
  
  return(list(solo_stats=solo_stats, solo_count, solo_count_plot, combined_plot=combined_plot))
}