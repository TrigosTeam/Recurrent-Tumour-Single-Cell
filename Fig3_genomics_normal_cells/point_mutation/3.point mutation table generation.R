library(R.utils)
library(vcfR)
library(ggplot2)
library(data.table)

point_mutation <-list()

SO_term <- list(
HIGH = c("transcript_ablation","splice_acceptor_variant","frameshift_variant","splice_donor_variant","stop_gained","transcript_amplification","start_lost","stop_lost"),
MODERATE = c("inframe_insertion","regulatory_region_ablation"   , "inframe_deletion", "missense_variant","protein_altering_variant"),
LOW = c("incomplete_terminal_codon_variant", "splice_region_variant","splice_donor_5th_base_variant","splice_donor_region_variant","splice_polypyrimidine_tract_variant","start_retained_variant","stop_retained_variant","synonymous_variant" ),
MODIFIER = c("3_prime_UTR_variant"         ,       "5_prime_UTR_variant"          ,      "coding_sequence_variant" ,          
            "downstream_gene_variant"    ,        "feature_elongation"      ,           "feature_truncation",                
            "intergenic_variant"      ,           "intron_variant"          ,           "mature_miRNA_variant",              
            "NMD_transcript_variant"   ,          "non_coding_transcript_exon_variant", "non_coding_transcript_variant",     
            "regulatory_region_amplification"  ,  "regulatory_region_variant",         
            "TF_binding_site_variant"     ,       "TFBS_ablation"          ,            "TFBS_amplification",                
            "upstream_gene_variant"     ))
format <- c("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|CADD_PHRED|CADD_RAW|Condel|FATHMM_converted_rankscore|FATHMM_pred|FATHMM_score|Interpro_domain|LRT_Omega|LRT_converted_rankscore|LRT_pred|LRT_score|MutationAssessor_pred|MutationAssessor_rankscore|MutationAssessor_score|MutationTaster_AAE|MutationTaster_converted_rankscore|MutationTaster_model|MutationTaster_pred|MutationTaster_score|PROVEAN_converted_rankscore|PROVEAN_pred|PROVEAN_score|Polyphen2_HDIV_pred|Polyphen2_HDIV_rankscore|Polyphen2_HDIV_score|Polyphen2_HVAR_pred|Polyphen2_HVAR_rankscore|Polyphen2_HVAR_score|SIFT4G_converted_rankscore|SIFT4G_pred|SIFT4G_score|SIFT_converted_rankscore|SIFT_pred|SIFT_score|clinvar_MedGen_id|clinvar_OMIM_id|clinvar_Orphanet_id|clinvar_clnsig|clinvar_hgvs|clinvar_id|clinvar_review|clinvar_trait|clinvar_var_source|gnomADg|gnomADg_AF|gnomADg_AF_raw|gnomADg_AF_fin_female|gnomADg_AF_afr_male|gnomADg_AF_afr|gnomADg_AF_eas_female|gnomADg_AF_afr_female|gnomADg_AF_fin_male|gnomADg_AF_nfe_female|gnomADg_AF_amr|gnomADg_AF_eas|gnomADg_AF_asj_male|gnomADg_AF_oth_female|gnomADg_AF_female|gnomADg_AF_eas_male|gnomADg_AF_nfe|gnomADg_AF_fin|gnomADg_AF_nfe_male|gnomADg_AF_asj_female|gnomADg_AF_asj|gnomADg_AF_oth|gnomADg_AF_male|gnomADg_AF_amr_male|gnomADg_AF_amr_female|gnomADg_AF_oth_male|gnomADg_AF_ami_female|gnomADg_AF_ami_male|gnomADg_AF_ami|gnomADg_AF_sas_female|gnomADg_AF_sas_male|gnomADg_AF_sas")
f <- strsplit(format, "[|]")[[1]]
paths <-system("realpath ~/Sequencing_updated_reference/*/*/vep/*intersect.vep.vcf.gz", intern = T)
paths_split <- sapply(strsplit(paths, split = "/"), function(x) paste(x[5], x[6], sep = "_"))
paths <- setNames(paths, paths_split)

point_mutation <-list()
for (i in paths_split){
  vcf <- read.vcfR(paths[i])
  #queryMETA(vcf, element = 'AD')
  fix_df <- as.data.frame(vcf@fix)
  info <- extract.info(vcf, "CSQ")
  mylist <- gsub("[|]", ",", info)
  mylist <- strsplit(mylist, ",")
  CSQlist <- lapply(seq_along(mylist), function(idx) {
    x <- mylist[[idx]]
    temp <- matrix(x, ncol = 142, byrow = T)
    colnames(temp) <- f
    return(temp)
  })
  
  SGT = extract.info(vcf,"SGT")
  DP = as.numeric(extract.info(vcf, "DP"))
  gt_dp <- extract.gt(vcf, element = "DP", as.numeric = T)
  gt_ad <- extract.gt(vcf, element = "AD", as.numeric = F)
  gt_af <- extract.gt(vcf, element = "AF", as.numeric = T)
  gt_df <- as.data.table(cbind(gt_dp, gt_ad, gt_af), keep.rownames = F)
  colnames(gt_df) <- c( "normal_DP", "tumor_DP", "normal_AD", "tumor_AD", "normal_AF", "tumor_AF")
  
  #print(paste(c(length(SGT), length(DP),nrow(gt_df), length(CSQlist))))
  SO_table <- rbindlist(lapply(seq_along(CSQlist), function(idx){
    temp <- cbind(fix_df[idx, 1:6], SGT = SGT[idx],DP = DP[idx], gt_df[idx,], CSQlist[[idx]])
    return(temp)
  }))
  #"SOMATIC;QSS=0;TQSS=1;NT=ref;QSS_NT=0;TQSS_NT=1;SGT=AA->AA;DP=263;MQ=35.64;MQ0=19;ReadPosRankSum=1.74;SNVSB=1.87;SomaticEVS=0.08;CombinedSomEVS=12.14;
  point_mutation[[i]] <- SO_table
  print(i) 
}
saveRDS(point_mutation, "point_mutation_table_intersection.Rds")
