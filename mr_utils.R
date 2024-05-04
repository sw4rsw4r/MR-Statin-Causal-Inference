library(dplyr)
library(vroom)
library(rtracklayer)
library(reticulate)
library(Matrix) # install.packages("Matrix", repos = "http://R-Forge.R-project.org")
library(prop.coloc)
library(ggplot2)
library(reshape2)
library(easyGgplot2) # https://github.com/kassambara/easyGgplot2
library(MendelianRandomization)
source("PCA-LIML-function.R")
np <- reticulate::import("numpy")

get_gene_region <- function(gene, window = 100000) {
  if (gene == "HMGCR") {
    ### target: HMGCR gene with REF: GRCh37
    # https://www.ncbi.nlm.nih.gov/gene/3156
    CHR <- "5"
    gene_start <- 74632993
    gene_end <- 74657941
    gene_name <- "ENSG00000113161"
  }
  if (gene == "LPL") {
    CHR <- "8"
    gene_start <- 19796764
    gene_end <- 19824770
    gene_name <- "ENSG00000175445"
  }
  if (gene == "ANGPTL4") {
    CHR <- "19"
    gene_start <- 8429039
    gene_end <- 8439254
    gene_name <- "ENSG00000167772"
  }
  if (gene == "PCSK9") {
    # https://www.ncbi.nlm.nih.gov/gene/255738
    CHR <- "1"
    gene_start <- 55505221
    gene_end <- 55530525
    gene_name <- "ENSG00000169174"
  }
  return(list(CHR = CHR, locus_lower = max(gene_start - window, 0), locus_upper = gene_end + window, gene_name = gene_name))
}

get_gene_region_hg38 <- function(gene, window = 100000) {
  if (gene == "HMGCR") {
    ### target: HMGCR gene with REF: GRCh37
    # https://www.ncbi.nlm.nih.gov/gene/3156
    CHR <- "5"
    gene_start <- 75336529
    gene_end <- 75362116
    gene_name <- "ENSG00000113161"
  }
  if (gene == "LPL") {
    CHR <- "8"
    gene_start <- 19939253
    gene_end <- 19967259
    gene_name <- "ENSG00000175445"
  }
  if (gene == "ANGPTL4") {
    CHR <- "19"
    gene_start <- 8364155
    gene_end <- 8374370
    gene_name <- "ENSG00000167772"
  }
  if (gene == "PCSK9") {
    # https://www.ncbi.nlm.nih.gov/gene/255738
    CHR <- "1"
    gene_start <- 55039548
    gene_end <- 55064852
    gene_name <- "ENSG00000169174"
  }
  return(list(CHR = CHR, locus_lower = max(gene_start - window, 0), locus_upper = gene_end + window, gene_name = gene_name))
}


load_CAD <- function(gene, window) {
  region <- get_gene_region(gene, window)

  # http://www.cardiogramplusc4d.org/data-downloads/
  # wget http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz
  # https://static-content.springer.com/esm/art%3A10.1038%2Fng.3913/MediaObjects/41588_2017_BFng3913_MOESM1_ESM.pdf
  fname <- "data/CAD/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz" # hg19
  tb <- vroom::vroom(fname)
  case_prob1 <- 10801 / (10801 + 137371)
  case_prob2 <- 60801 / (60801 + 123504)
  case_prob3 <- 42335 / (42335 + 78240)

  res <- with(tb, dplyr::tibble(
    beta = ifelse(effect_allele_freq < .5, logOR, -logOR),
    se = se_gc,
    varbeta = se_gc^2,
    snp = snptestid,
    effect = ifelse(effect_allele_freq < .5, effect_allele, noneffect_allele),
    other = ifelse(effect_allele_freq < .5, noneffect_allele, effect_allele),
    chrom = chr,
    position = bp_hg19,
    MAF = ifelse(effect_allele_freq < .5, effect_allele_freq, 1 - effect_allele_freq),
    type = "cc",
    s = mean(c(case_prob1, case_prob2, case_prob3)),
    nsample = n_samples
  )) %>% dplyr::filter(chrom == region$CHR, position >= region$locus_lower, position <= region$locus_upper)
  return(res)
}


load_T2D <- function(type = "EUR", gene, window) {
  region <- get_gene_region(gene, window)

  # https://diagram-consortium.org/downloads.html
  # Download Ancestry specific GWAS meta-analysis summary statistics: European
  # Published in Mahajan et al (2022)
  fname <- paste0("data/T2D/DIAMANTE-", type, ".sumstat.txt")
  if (type == "EUR") {
    N <- 492191 * 0.511
  }
  if (type == "TA") {
    N <- 492191
  }
  case_prop <- 180834 / (180834 + 1159055)
  df <- vroom::vroom(fname)
  res <- with(df, dplyr::tibble(
    beta = ifelse(effect_allele_frequency < .5, `Fixed-effects_beta`, -`Fixed-effects_beta`),
    se = `Fixed-effects_SE`,
    varbeta = `Fixed-effects_SE`^2,
    snp = rsID,
    effect = ifelse(effect_allele_frequency < .5, toupper(effect_allele), toupper(other_allele)),
    other = ifelse(effect_allele_frequency < .5, toupper(other_allele), toupper(effect_allele)),
    chrom = `chromosome(b37)`,
    position = `position(b37)`,
    MAF = ifelse(effect_allele_frequency < .5, effect_allele_frequency, 1 - effect_allele_frequency),
    type = "cc",
    s = case_prop,
    nsample = N
  )) %>% dplyr::filter(!is.na(snp), chrom == region$CHR, position >= region$locus_lower, position <= region$locus_upper)
  return(res)
}

load_GWAS <- function(pheno, gene, data_type = "quant", case_prop = NA, window) {
  # https://www.ebi.ac.uk/gwas/
  if (pheno == "Pancreatitis") {
    fname <- "data/GWAS/34662886-GCST90077706-EFO_0000278.h.tsv.gz"
    nsample <- 329052
    data_type <- "cc"
    case_prop <- 868 / (868 + 328184)
  }
  if (pheno == "Hepatitis") {
    fname <- "data/GWAS/34737426-GCST90041716-HP_0012115-Build37.f.tsv.gz"
    nsample <- 456348
    data_type <- "cc"
    case_prop <- 236 / (236 + 456112)
  }
  if (pheno == "Myositis") {
    fname <- "data/GWAS/34737426-GCST90043770-EFO_0000783-Build37.f.tsv.gz"
    nsample <- 456348
    data_type <- "cc"
    case_prop <- 226 / (226 + 456122)
  }
  if (pheno == "Myalgia") {
    fname <- "data/GWAS/34662886-GCST90080544-EFO_0000408-Build38.f.tsv.gz"
    nsample <- 384034
    data_type <- "cc"
    case_prop <- 1350 / (1350 + 382684)
  }
  if (pheno == "Alzheimer") {
    fname <- "data/GWAS/35379992-GCST90027158-MONDO_0004975-Build38.f.tsv.gz"
    nsample <- 487511
    data_type <- "cc"
    case_prop <- 85934 / (85934 + 401577)
  }
  if (pheno == "Parkinson") {
    fname <- "data/GWAS/34737426-GCST90043734-EFO_0002508-Build37.f.tsv.gz"
    nsample <- 456348
    data_type <- "cc"
    case_prop <- 294 / (294 + 456054)
  }
  if (pheno == "Atherosclerosis") {
    fname <- "data/GWAS/34737426-GCST90044006-EFO_0003914-Build37.f.tsv.gz"
    nsample <- 456348
    data_type <- "cc"
    case_prop <- 106 / (106 + 456242)
  }
  if (pheno == "Amyotrophic lateral sclerosis") {
    fname <- "data/GWAS/34873335-GCST90027164-MONDO_0004976-Build37.f.tsv.gz"
    nsample <- 138086
    data_type <- "cc"
    case_prop <- 27205 / (27205 + 110881)
  }
  if (pheno == "Acute Insulin response") {
    fname <- "data/GWAS/28490609-GCST004575-EFO_0006831.h.tsv.gz" # GRCh38
    nsample <- 4765
  }
  if (pheno == "Fasting Insulin") {
    fname <- "data/GWAS/22581228-GCST005185-EFO_0004466.h.tsv.gz" # GRCh38
    nsample <- 51750
  }
  if (pheno == "Fasting Glucose") {
    fname <- "data/GWAS/22581228-GCST005186-EFO_0004465.h.tsv.gz"
    nsample <- 58074
  }
  if (pheno == "BMI") {
    fname <- "data/GWAS/GCST90179150_buildGRCh37.tsv" # GRCh37
    nsample <- 694649
  }
  if (pheno == "LDL-C") {
    fname <- "data/GWAS/GCST90239658_buildGRCh37.tsv" # GRCh37
    nsample <- 1320016
  }
  if (pheno == "HDL-C") {
    # correct rsID using LDL-C data
    # Rscript data/GWAS/correct_HDLC_rsID.R
    fname <- "data/GWAS/GCST90239652_buildGRCh37_corrected_rsID.tsv" # GRCh37
    nsample <- 1320016
  }
  if (pheno == "Triglyceride") {
    fname <- "data/GWAS/GCST90239664_buildGRCh37.tsv" # GRCh37
    nsample <- 1320016
  }
  if (pheno == "WH ratio") {
    fname <- "data/GWAS/GCST90020025.h.tsv.gz" # GRCh38
    nsample <- 219872
  }
  if (pheno == "Leptin") {
    fname <- "data/GWAS/33067605-GCST90012076-EFO_0005000.h.tsv.gz" # GRCh38
    nsample <- 21758
  }
  if (pheno == "Sterol") {
    fname <- "data/GWAS/34503513-GCST90060133-EFO_0010231-Build37.f.tsv.gz"
    nsample <- 13814
  }
  if (pheno == "Cortisol") {
    fname <- "data/GWAS/GCST90200378_buildGRCh38.tsv.gz"
    nsample <- 8193
  }
  if (pheno == "Testosterone") {
    fname <- "data/GWAS/GCST90239819_buildGRCh37.tsv.gz"
    nsample <- 176212
  }
  if (pheno == "Estradiol") {
    fname <- "data/GWAS/34255042-GCST90020092-EFO_0004697-Build38.f.tsv.gz"
    nsample <- 163985
    data_type <- "cc"
    case_prop <- 37461 / (37461 + 126524)
  }
  if (pheno == "Vitamin D") {
    fname <- "data/GWAS/GCST90162562_buildGRCh37.tsv"
    nsample <- 64988
  }
  if (pheno == "Bile acid") {
    fname <- "data/GWAS/34503513-GCST90060135-EFO_0010231-Build37.f.tsv.gz"
    nsample <- 13814
  }
  if (pheno == "Aldosterone") {
    fname <- "data/GWAS/GCST90012609_buildGRCh37.tsv.gz"
    nsample <- 1128
  }
  if (pheno == "CRP") {
    fname <- "data/GWAS/33328453-GCST90019446-EFO_0004458-Build37.f.tsv.gz"
    nsample <- 10708
  }
  if (pheno == "HbA1c") {
    fname <- "data/GWAS/34059833-GCST90002244-EFO_0004541-Build37.f.tsv.gz" # GRCh38
    nsample <- 146806
  }
  if (pheno == "Ubiquinone") {
    fname <- "data/GWAS/35668104-GCST90024608-EFO_0021486-Build37.f.tsv.gz"
    nsample <- 4492
  }
  if (pheno == "Erectile dysfunction") {
    fname <- "data/GWAS/34737426-GCST90044277-EFO_0004234-Build37.f.tsv.gz"
    nsample <- 208808
    data_type <- "cc"
    case_prop <- 357 / (357 + 208451)
  }
  if (pheno == "Nonalcoholic fatty liver disease") {
    fname <- "data/GWAS/GCST90275041.h.tsv.gz"
    nsample <- 32941
    data_type <- "cc"
    case_prop <- 6623 / (6623 + 26318)
  }

  df <- vroom(fname)
  genomic_version <- ifelse(grepl("GRCh37", fname) || grepl("Build37", fname), "GRCh37", "GRCh38")

  if (fname == "data/GWAS/34255042-GCST90020092-EFO_0004697-Build38.f.tsv.gz") genomic_version <- "GRCh37"

  if (genomic_version == "GRCh38") {
    region <- get_gene_region_hg38(gene, window)
    df_filt <- df %>% dplyr::filter(
      chromosome == region$CHR,
      base_pair_location >= region$locus_lower,
      base_pair_location <= region$locus_upper
    )
    if (fname == "data/GWAS/GCST90020025.h.tsv.gz") df_filt <- df_filt %>% rename(rs_id = "variant_id")
    if (fname == "data/GWAS/GCST90275041.h.tsv.gz") df_filt <- df_filt %>% mutate(variant_id = rs_id)
    if (fname %in% c("data/GWAS/34662886-GCST90080544-EFO_0000408-Build38.f.tsv.gz", "data/GWAS/34662886-GCST90077706-EFO_0000278.h.tsv.gz")) {
      df_filt <- df_filt %>% rename(odds_ratio = "beta")
    }

    gr <- with(df_filt, GenomicRanges::GRanges(
      seqnames = paste0("chr", chromosome),
      beta = ifelse(effect_allele_frequency < .5, beta, -beta),
      se = standard_error,
      snp = variant_id,
      effect = ifelse(effect_allele_frequency < .5, effect_allele, other_allele),
      other = ifelse(effect_allele_frequency < .5, other_allele, effect_allele),
      MAF = ifelse(effect_allele_frequency < .5, effect_allele_frequency, 1 - effect_allele_frequency),
      ranges = IRanges(start = base_pair_location - 1, end = base_pair_location)
    ))

    chain_hg38Tohg19 <- rtracklayer::import.chain("data/liftOver/hg38ToHg19.over.chain")
    lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)
    res <- with(as.data.frame(lifted_over), dplyr::tibble(beta,
      se,
      varbeta = se^2,
      snp,
      effect = toupper(effect),
      other = toupper(other),
      chrom = sub("chr", "", seqnames),
      position = end,
      MAF,
      type = data_type,
      s = case_prop,
      nsample = nsample
    ))
  } else if (genomic_version == "GRCh37") {
    region <- get_gene_region(gene, window)
    df_filt <- df %>% dplyr::filter(
      chromosome == region$CHR,
      base_pair_location >= region$locus_lower,
      base_pair_location <= region$locus_upper
    )
    if (fname == "data/GWAS/GCST90012609_buildGRCh37.tsv.gz") {
      df_filt <- df_filt %>% rename(BETA = "beta", SE = "standard_error", A1 = "effect_allele", A2 = "other_allele", AF1 = "effect_allele_frequency")
    }
    if (fname == "data/GWAS/GCST90162562_buildGRCh37.tsv") {
      df_filt <- df_filt %>% rename(SNP = "variant_id")
    }
    res <- with(df_filt, dplyr::tibble(
      beta = ifelse(effect_allele_frequency < 0.5, beta, -beta),
      se = standard_error,
      varbeta = standard_error^2,
      snp = variant_id,
      effect = toupper(ifelse(effect_allele_frequency < 0.5, effect_allele, other_allele)),
      other = toupper(ifelse(effect_allele_frequency < 0.5, other_allele, effect_allele)),
      chrom = chromosome,
      position = base_pair_location,
      MAF = ifelse(effect_allele_frequency < 0.5, effect_allele_frequency, 1 - effect_allele_frequency),
      type = data_type,
      s = case_prop,
      nsample = nsample
    ))
  }
  res <- res %>%
    dplyr::filter(!is.na(beta) & !is.na(effect) & !is.na(MAF)) %>%
    distinct()
  return(res)
}

load_pQTL <- function(pheno, gene, data_type = "quant", case_prop = NA, window) {
  # https://www.decode.com/summarydata/
  # Grímur Hjörleifsson Eldjarn, Egil Ferkingstad et al. Large-scale plasma proteomics comparisons through genetics and disease associations
  if (pheno == "HMGCR") {
    fname <- "data/pQTL/deCODE/2023_Large_scale_plasma/Proteomics_SMP_PC0_5230_99_HMGCR_HMGR.txt.gz" # GRCh38
  }
  if (pheno == "Phosphomevalonate kinase") {
    fname <- "data/pQTL/deCODE/2023_Large_scale_plasma/Proteomics_SMP_PC0_12450_42_PMVK_PMVK.txt.gz"
  }

  df <- vroom::vroom(fname)
  df <- df %>%
    rename(
      Chrom = "chromosome", Pos = "base_pair_location", Beta = "beta", SE = "standard_error", rsids = "variant_id",
      A1 = "effect_allele", A0 = "other_allele", ImpFreqA1 = "effect_allele_frequency"
    ) %>%
    mutate(chromosome = sub("chr", "", chromosome))

  nsample <- mean(df$N)

  region <- get_gene_region_hg38(gene, window)
  df_filt <- df %>% dplyr::filter(
    chromosome == region$CHR,
    base_pair_location >= region$locus_lower,
    base_pair_location <= region$locus_upper
  )

  gr <- with(df_filt, GenomicRanges::GRanges(
    seqnames = paste0("chr", chromosome),
    beta = ifelse(effect_allele_frequency < .5, beta, -beta),
    se = standard_error,
    snp = variant_id,
    effect = ifelse(effect_allele_frequency < .5, effect_allele, other_allele),
    other = ifelse(effect_allele_frequency < .5, other_allele, effect_allele),
    MAF = ifelse(effect_allele_frequency < .5, effect_allele_frequency, 1 - effect_allele_frequency),
    ranges = IRanges(start = base_pair_location - 1, end = base_pair_location)
  ))
  chain_hg38Tohg19 <- rtracklayer::import.chain("data/liftOver/hg38ToHg19.over.chain")
  lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)
  res <- with(as.data.frame(lifted_over), dplyr::tibble(beta,
    se,
    varbeta = se^2,
    snp,
    effect = toupper(effect),
    other = toupper(other),
    chrom = sub("chr", "", seqnames),
    position = end,
    MAF,
    type = data_type,
    s = case_prop,
    nsample = nsample
  ))
  res <- res %>%
    dplyr::filter(!is.na(beta) & !is.na(effect) & !is.na(MAF)) %>%
    distinct()
  return(res)
}


load_eQTL <- function(cell_type, gene, window = 100000) {
  # https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv

  if (cell_type == "Brain caudate") {
    fname <- "data/eQTL/QTD000156.all.tsv.gz"
  } else if (cell_type == "Small intestine") {
    fname <- "data/eQTL/QTD000321.all.tsv.gz"
  } else if (cell_type == "Liver") {
    fname <- "data/eQTL/QTD000266.all.tsv.gz"
  } else if (cell_type == "Pancreas") {
    fname <- "data/eQTL/QTD000296.all.tsv.gz"
  } else if (cell_type == "Muscle") {
    fname <- "data/eQTL/QTD000281.all.tsv.gz"
  } else if (cell_type == "Blood") {
    fname <- "data/eQTL/QTD000356.all.tsv.gz"
  } else if (cell_type == "Fat") {
    fname <- "data/eQTL/QTD000534.all.tsv.gz"
  } else if (cell_type == "Adrenal gland") {
    fname <- "data/eQTL/QTD000126.all.tsv.gz"
  } else if (cell_type == "Testis") {
    fname <- "data/eQTL/QTD000336.all.tsv.gz"
  } else if (cell_type == "Ovary") {
    fname <- "data/eQTL/QTD000291.all.tsv.gz"
  } else if (cell_type == "Skin sun exposed") {
    fname <- "data/eQTL/QTD000316.all.tsv.gz"
  } else if (cell_type == "Skin not sun exposed") {
    fname <- "data/eQTL/QTD000311.all.tsv.gz"
  } else if (cell_type == "Adipose subcutaneous") {
    fname <- "data/eQTL/QTD000116.all.tsv.gz"
  } else if (cell_type == "Adipose visceral") {
    fname <- "data/eQTL/QTD000121.all.tsv.gz"
  } else if (cell_type == "Brain cerebellum") {
    fname <- "data/eQTL/QTD000166.all.tsv.gz"
  } else if (cell_type == "Brain cortex") {
    fname <- "data/eQTL/QTD000171.all.tsv.gz"
  } else if (cell_type == "Brain hippocampus") {
    fname <- "data/eQTL/QTD000181.all.tsv.gz"
  } else if (cell_type == "Brain putamen") {
    fname <- "data/eQTL/QTD000196.all.tsv.gz"
  } else if (cell_type == "Brain spinal cord") {
    fname <- "data/eQTL/QTD000201.all.tsv.gz"
  } else if (cell_type == "Brain substantia nigra") {
    fname <- "data/eQTL/QTD000206.all.tsv.gz"
  } else {
    stop("unexpected cell_type")
  }

  region <- get_gene_region_hg38(gene, window)
  cmd <- with(region, paste0("tabix ", fname, " ", CHR, ":", locus_lower, "-", locus_upper))
  df <- read.delim(pipe(cmd), header = F)
  colnames(df) <- c("ensg_id", "chrom", "position", "other", "effect", "chrpos", "variant_ma_samples", "MAF", "pval", "beta", "se", "is_snp", "ac", "an", "r2", "geneid1", "geneid2", "median_tpm", "snp")
  df_filt <- df %>% filter(ensg_id == region$gene_name)
  if (nrow(df_filt) == 0) {
    return(NULL)
  }

  gr <- with(df_filt, GenomicRanges::GRanges(
    seqnames = paste0("chr", chrom),
    beta = beta,
    se = se,
    snp = snp,
    effect = effect,
    other = other,
    MAF = MAF,
    ranges = IRanges(start = position - 1, end = position)
  ))

  chain_hg38Tohg19 <- rtracklayer::import.chain("data/liftOver/hg38ToHg19.over.chain")
  lifted_over <- rtracklayer::liftOver(gr, chain_hg38Tohg19)
  res <- with(as.data.frame(lifted_over), dplyr::tibble(beta,
    se,
    varbeta = se^2,
    snp,
    effect,
    other,
    chrom = sub("chr", "", seqnames),
    position = end,
    MAF,
    type = "quant",
    s = NA,
    nsample = 7526
  ))
  return(res)
}


lookup_ncbi <- function(query_in, ID, max_search = 10) {
  cnt <- 0
  query_out <- NULL
  while (is.null(query_out)) {
    query_out <- tryCatch(
      {
        readLines(query_in)
      },
      error = function(cond) NULL
    )
    if (cnt > max_search) break
    message("Search ", ID, " .. ", cnt)
    cnt <- cnt + 1
  }
  return(query_out)
}

fill_rsID <- function(lst_data_combined) {
  if (sum(sapply(lst_data_combined, function(x) sum(is.na(x$snp)))) == 0) {
    return(lst_data_combined)
  }

  ID_map <- list()
  for (idx1 in 1:length(lst_data_combined)) {
    print(names(lst_data_combined)[idx1])
    lst_data_combined[[idx1]] <- lst_data_combined[[idx1]] %>% dplyr::mutate(posID = paste0(chrom, ":", position))
    for (idx2 in 1:nrow(lst_data_combined[[idx1]])) {
      posID <- lst_data_combined[[idx1]][idx2, ]$posID
      rsID <- lst_data_combined[[idx1]][idx2, ]$snp
      if (is.na(rsID)) next
      ID_map[[posID]] <- unique(c(ID_map[[posID]], rsID))
    }
  }
  if (length(ID_map) == 0) {
    ID_map <- list()
    for (idx1 in 1:length(lst_data_combined)) {
      print(names(lst_data_combined)[idx1])
      lst_data_combined[[idx1]] <- lst_data_combined[[idx1]] %>% dplyr::mutate(posID = paste0(chrom, ":", position))
      for (idx2 in 1:nrow(lst_data_combined[[idx1]])) {
        posID <- lst_data_combined[[idx1]][idx2, ]$posID
        chrom <- lst_data_combined[[idx1]][idx2, ]$chrom
        pos <- lst_data_combined[[idx1]][idx2, ]$position
        query <- paste0("https://www.ncbi.nlm.nih.gov/snp/?term=", chrom, "%3A", pos)
        query_out <- lookup_ncbi(query, posID)
        line_selected <- grep("snp_info=", query_out, value = T)
        if (length(line_selected) == 0) {
          next
        } else {
          line_selected <- line_selected[[1]]
        }
        rsID <- strsplit(line_selected, 'snp_info=\"|:')[[1]][11]
        ID_map[[posID]] <- rsID
      }
    }
  }
  ID_map <- ID_map[sapply(ID_map, length) == 1]
  for (idx1 in 1:length(lst_data_combined)) {
    print(names(lst_data_combined)[idx1])
    for (idx2 in 1:nrow(lst_data_combined[[idx1]])) {
      posID <- lst_data_combined[[idx1]][idx2, ]$posID
      rsID <- lst_data_combined[[idx1]][idx2, ]$snp
      if (is.na(rsID) & posID %in% names(ID_map)) {
        lst_data_combined[[idx1]][idx2, "snp"] <- ID_map[[posID]]
      }
    }
    lst_data_combined[[idx1]] <- lst_data_combined[[idx1]] %>% dplyr::filter(!is.na(snp))
  }
  return(lst_data_combined)
}

count_snps <- function(lst_data_combined, dir_output) {
  lst <- names(lst_data_combined)
  do_append <- F
  for (N in 1:length(lst)) {
    set <- as.matrix(combn(lst, N))
    for (idx in 1:ncol(set)) {
      n_snps <- length(Reduce("intersect", lapply(lst_data_combined[set[, idx]], function(x) x$snp)))
      line_to_write <- paste(paste(set[, idx], collapse = "_"), n_snps)
      write.table(line_to_write, file.path(dir_output, "00_n_snps.txt"), quote = F, row.names = F, col.names = F, append = do_append)
      do_append <- T
    }
  }
}

load_ld_mat_UKBB <- function(path_ld_mat, ld_merged) {
  fname_save <- sub(".npz$", ".RDS", path_ld_mat)
  if (file.exists(fname_save)) {
    mat <- readRDS(fname_save)
  } else {
    npz_file <- np$load(path_ld_mat)
    i <- as.numeric(npz_file$f[["row"]])
    j <- as.numeric(npz_file$f[["col"]])
    v <- as.numeric(npz_file$f[["data"]])
    dims <- as.numeric(npz_file$f[["shape"]])
    mat <- Matrix::sparseMatrix(i, j, x = v, index1 = FALSE, dims = dims)
    saveRDS(mat, file = fname_save)
  }

  filt <- sort(ld_merged$idx)
  mat_filt <- as.matrix(mat[filt, filt])
  mat_t <- t(mat_filt)
  lower_tri <- lower.tri(mat_filt, diag = T)
  mat_t[lower_tri] <- mat_t[lower_tri] + mat_filt[lower_tri]
  rownames(mat_t) <- colnames(mat_t) <- ld_merged$snp
  return(mat_t)
}

harmonize_ld <- function(ld_mat, ld_merged) {
  flip <- with(ld_merged, ifelse(effect_1 == effect_2, 1, -1))
  ld_mat_corrected <- ld_mat * flip %o% flip
  return(ld_mat_corrected)
}

harmonize <- function(lst_data, path_ld_mat, path_ld_meta) {
  names_risk_factor <- names(lst_data$risk_factor)
  names_outcome <- names(lst_data$outcome)

  lst_data_combined <- c(lst_data$risk_factor, lst_data$outcome)
  lst_data_combined <- fill_rsID(lst_data_combined)

  count_snps(lst_data_combined, dir_output)
  if (any(sapply(lst_data_combined, nrow) < 2)) {
    return(NULL)
  }

  ld_meta <- vroom::vroom(path_ld_meta, delim = "\t", comment = "#")
  ld_meta <- with(ld_meta, data.frame(snp = rsid, effect = allele1, other = allele2, idx = 1:nrow(ld_meta)))

  ids_to_keep <- Reduce("intersect", lapply(lst_data_combined, function(x) x$snp))

  n_samples <- sapply(lst_data_combined, function(x) mean(x$nsample))
  max_sample <- max(n_samples)
  selected_factor <- sort(names(which(n_samples == max_sample)))[1]

  res <- list()
  ids_do_not_match <- NULL
  for (this_factor in setdiff(names(lst_data_combined), selected_factor)) {
    merged <- lst_data_combined[[selected_factor]] %>%
      inner_join(lst_data_combined[[this_factor]], by = "snp", suffix = c("_1", "_2")) %>%
      dplyr::filter(snp %in% ids_to_keep) %>%
      dplyr::mutate(
        effect_2_ori = effect_2,
        beta_2 = ifelse(other_2 == effect_1, -beta_2, beta_2),
        MAF_2 = ifelse(other_2 == effect_1, 1 - MAF_2, MAF_2),
        effect_2 = ifelse(other_2 == effect_1, other_2, effect_2),
        other_2 = ifelse(other_2 == effect_1, effect_2_ori, other_2)
      )

    ids_do_not_match <- unique(c(ids_do_not_match, merged %>% dplyr::filter(effect_1 != effect_2) %>% pull(snp)))

    if (length(res) == 0) {
      res[[1]] <- merged %>% dplyr::select(snp, grep("_1", colnames(merged), value = T))
      colnames(res[[1]]) <- sub("_1$", "", colnames(res[[1]]))
    }
    dat <- merged %>% dplyr::select(snp, grep("_2", colnames(merged), value = T))
    colnames(dat) <- sub("_2$", "", colnames(dat))
    res[[length(res) + 1]] <- dat
  }
  names(res) <- c(selected_factor, setdiff(names(lst_data_combined), selected_factor))

  # filter LD SNPs
  ld_merged <- as_tibble(res[[selected_factor]]) %>% left_join(ld_meta, by = "snp", suffix = c("_1", "_2"))
  ld_merged <- ld_merged %>% dplyr::filter(!is.na(effect_2) & (effect_1 == effect_2 | effect_1 == other_2))
  ids_do_not_match_with_ld <- setdiff(ids_to_keep, ld_merged$snp)

  ids_to_remove <- unique(c(ids_do_not_match_with_ld, ids_do_not_match))

  ld_merged <- ld_merged %>% dplyr::filter(!snp %in% ids_to_remove)
  ld_merged <- ld_merged[!duplicated(ld_merged$snp), ]
  for (id in names(res)) {
    res[[id]] <- res[[id]] %>%
      dplyr::filter(!duplicated(snp)) %>%
      as_tibble() %>%
      dplyr::filter(!snp %in% ids_to_remove) %>%
      as.list()
  }

  for (id in c(names_risk_factor, names_outcome)) {
    res[[id]]$type <- unique(res[[id]]$type)
    res[[id]]$N <- mean(res[[id]]$nsample)
    if (is.na(res[[id]]$s[1])) res[[id]]$s <- NULL
  }

  message(paste0("SNPs removed: ", paste(ids_to_remove, collapse = ",")))
  write.table(paste("Final_nSNPs", nrow(ld_merged)),
    file.path(dir_output, "00_n_snps.txt"),
    quote = F, row.names = F, col.names = F, append = T
  )

  # read and correct LD
  ld_mat <- load_ld_mat_UKBB(path_ld_mat, ld_merged)

  res$ld <- harmonize_ld(ld_mat, ld_merged)
  res$names$risk_factors <- names_risk_factor
  res$names$outcome <- names_outcome
  return(res)
}


run_propcoloc <- function(res, dir_output, RF1, RF2) {
  this_prune <- 0.4
  this_J <- 10
  gene_of_interest <- basename(dirname(dirname(dir_output)))
  # if (gene_of_interest == 'HMGCR' & RF1 == "LDL-C" & RF2 == "Aldosterone") {
  #   this_prune = 0.2;
  #   this_J = 5;
  # }
  lst_eQTL <- c("Brain caudate", "Small intestine", "Liver", "Pancreas", "Muscle", "Blood", "Fat", "Adrenal gland", "Testis", "Ovary", "Skin sun exposed", "Skin not sun exposed", "Adipose subcutaneous", "Adipose visceral", "Brain cerebellum", "Brain cortex", "Brain hippocampus", "Brain putamen", "Brain spinal cord", "Brain substantia nigra")
  lst_lipids <- c("LDL-C", "HDL-C", "Triglyceride")

  prop.coloc.res <- NULL
  n_try <- 0
  while (is.null(prop.coloc.res)) {
    if ((RF1 %in% lst_eQTL & RF2 %in% lst_eQTL) || (RF1 %in% lst_lipids & RF2 %in% lst_lipids)) {
      prop.coloc.res <- tryCatch(
        {
          prop.coloc::prop.coloc(
            b1 = res[[RF1]]$beta, se1 = res[[RF1]]$se, b2 = res[[RF2]]$beta, se2 = res[[RF2]]$se,
            n = res[[RF1]]$N,
            ld = res$ld, figs = TRUE, traits = c(RF1, RF2),
            prune = this_prune, J = this_J
          )
        },
        error = function(e) {
          return(NULL)
        }
      )
    } else {
      prop.coloc.res <- tryCatch(
        {
          prop.coloc::prop.coloc(
            b1 = res[[RF1]]$beta, se1 = res[[RF1]]$se, b2 = res[[RF2]]$beta, se2 = res[[RF2]]$se,
            n = c(res[[RF1]]$N, res[[RF2]]$N),
            ld = res$ld, figs = TRUE, traits = c(RF1, RF2),
            prune = this_prune, J = this_J
          )
        },
        error = function(e) {
          return(NULL)
        }
      )
    }
    if (is.null(prop.coloc.res) & n_try == 0) {
      this_prune <- 0.2
      this_J <- 5
      n_try <- n_try + 1
    } else if (is.null(prop.coloc.res) & n_try == 1) {
      prop.coloc.res <- list(p_full = NA, p_cond = NA, LM_full = NA, LM_cond = NA)
    }
  }
  prop.coloc.res$J <- this_J
  prop.coloc.res$prune <- this_prune
  prop.coloc.res$rsIDs <- res[[1]]$snp
  saveRDS(prop.coloc.res, file = paste0(dir_output, "/03_prop.coloc_", RF1, "_", RF2, ".RDS"))
  write.table(as.data.frame(prop.coloc.res[sapply(prop.coloc.res, length) == 1]),
    paste0(dir_output, "/03_prop.coloc_", RF1, "_", RF2, ".txt"),
    quote = F, row.names = F, col.names = T, sep = "\t"
  )
}

select_candidate <- function(gene_of_interest, this_window) {
  lst_files <- list.files(paste0("results/window_", this_window, "/", gene_of_interest, "/propcoloc"), full.names = T)
  lst_files <- lst_files[base::sapply(lst_files, function(x) file.info(x)$isdir)]

  df_res <- NULL
  for (this_dir in lst_files) {
    this_pair <- unlist(strsplit(basename(this_dir), "_"))
    this_group <- paste(this_pair, collapse = "_")
    fname1 <- paste0(this_dir, "/03_prop.coloc_", this_pair[1], "_", this_pair[2], ".RDS")
    fname2 <- paste0(this_dir, "/03_prop.coloc_", this_pair[2], "_", this_pair[1], ".RDS")
    if (file.exists(fname1)) {
      df1 <- with(readRDS(fname1), data.frame(group = this_group, fac1 = this_pair[1], fac2 = this_pair[2], p_full, p_cond, LM_full, LM_cond))
      df2 <- with(readRDS(fname2), data.frame(p_full_rev = p_full, p_cond_rev = p_cond, LM_full_rev = LM_full, LM_cond_rev = LM_cond))
      df_res <- rbind(df_res, cbind(df1, df2))
    }
  }
  df_res_subset <- df_res %>%
    filter((p_cond < .05 & LM_cond < .05) | (p_cond_rev < .05 & LM_cond_rev < .05))

  write.table(df_res_subset, paste0("results/window_", this_window, "/", gene_of_interest, "/selected_pairs_based_on_propcoloc_union.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
  return(df_res_subset)
}

get_table_from_pca_res <- function(pca_res, n_risk_factors, non_robust = F) {
  if (non_robust) {
    ci <- with(pca_res, data.frame(
      estimate = liml.nr,
      lower = liml.nr - qnorm(1 - 0.05 / 2) * se.liml.nr,
      upper = liml.nr + qnorm(1 - 0.05 / 2) * se.liml.nr,
      `p-value` = 2 * (1 - pnorm(abs(liml.nr / se.liml.nr))),
      PCs = rep(factors, n_risk_factors),
      `OID test p-value` = rep(1 - pchisq(Q.nr, factors - (n_risk_factors + 1)), n_risk_factors)
    ))
  } else {
    ci <- with(pca_res, data.frame(
      estimate = liml,
      lower = liml - qnorm(1 - 0.05 / 2) * se.liml,
      upper = liml + qnorm(1 - 0.05 / 2) * se.liml,
      `p-value` = 2 * (1 - pnorm(abs(liml / se.liml))),
      PCs = rep(factors, n_risk_factors),
      `OID test p-value` = rep(1 - pchisq(Q, factors - (n_risk_factors + 1)), n_risk_factors)
    ))
  }
  return(round(ci, 3))
}

compute_Fstat <- function(MRInputObj, nx, ny, thres, n_PCs, multivariate) {
  use_n_PCs <- F
  mr_pcgmm_Fstat <- NULL
  n_try <- 0
  while (is.null(mr_pcgmm_Fstat) & n_try <= 3) {
    if (multivariate) {
      if (use_n_PCs) {
        mr_pcgmm_Fstat <- tryCatch(
          {
            Fstat1 <- mr_mvpcgmm(MRInputObj, nx = nx, ny = ny, r = n_PCs, robust = T)@CondFstat
            Fstat2 <- mr_mvpcgmm(MRInputObj, nx = nx, ny = ny, r = n_PCs, robust = F)@CondFstat
            df_Fstat <- data.frame(mvpcgmm = Fstat1, mvpcgmm_nr = Fstat2)
            df_Fstat
          },
          error = function(e) {
            return(NULL)
          }
        )
      } else {
        mr_pcgmm_Fstat <- tryCatch(
          {
            Fstat1 <- mr_mvpcgmm(MRInputObj, nx = nx, ny = ny, thres = thres, robust = T)@CondFstat
            Fstat2 <- mr_mvpcgmm(MRInputObj, nx = nx, ny = ny, thres = thres, robust = F)@CondFstat
            df_Fstat <- data.frame(mvpcgmm = Fstat1, mvpcgmm_nr = Fstat2)
            df_Fstat
          },
          error = function(e) {
            return(NULL)
          }
        )
      }
    } else {
      if (use_n_PCs) {
        mr_pcgmm_Fstat <- tryCatch(
          {
            Fstat1 <- mr_pcgmm(MRInputObj, nx = nx, ny = ny, r = n_PCs, robust = T)@Fstat
            Fstat2 <- mr_pcgmm(MRInputObj, nx = nx, ny = ny, r = n_PCs, robust = F)@Fstat
            c(Fstat1, Fstat2)
          },
          error = function(e) {
            return(NULL)
          }
        )
      } else {
        mr_pcgmm_Fstat <- tryCatch(
          {
            Fstat1 <- mr_pcgmm(MRInputObj, nx = nx, ny = ny, thres = thres, robust = T)@Fstat
            Fstat2 <- mr_pcgmm(MRInputObj, nx = nx, ny = ny, thres = thres, robust = F)@Fstat
            c(Fstat1, Fstat2)
          },
          error = function(e) {
            return(NULL)
          }
        )
      }
    }
    if (!is.null(mr_pcgmm_Fstat)) {
      return(mr_pcgmm_Fstat)
    } else {
      use_n_PCs <- T
      if (n_try != 0) {
        n_PCs <- n_PCs + 1
      }
      n_try <- n_try + 1
    }
  }
  return(c(NA, NA))
}


run_PCA_liml <- function(res, dir_output, cor.x = NULL) {
  res_PCA_liml <- NULL
  ld <- res$ld
  names_risk_factors <- res$names$risk_factors
  names_outcome <- res$names$outcome

  bx <- sapply(res[names_risk_factors], function(x) x$beta)
  sx <- sapply(res[names_risk_factors], function(x) x$se)
  nx <- sapply(res[names_risk_factors], function(x) mean(x$nsample))

  by <- sapply(res[names_outcome], function(x) x$beta)
  sy <- sapply(res[names_outcome], function(x) x$se)
  ny <- sapply(res[names_outcome], function(x) mean(x$nsample))

  pca.no <- function(thres) {
    a <- 1 / ((nx[1] * sx[, 1]^2) + bx[, 1]^2)
    A <- (sqrt(a) %*% t(sqrt(a))) * ld
    return(which(cumsum(prcomp(A, scale = FALSE)$sdev^2 / sum((prcomp(A, scale = FALSE)$sdev^2))) > thres)[1])
  }

  pca_thres <- NULL
  for (thres in c(0.9999, 0.999, 0.99)) {
    n_PCs <- pca.no(thres)
    pca_thres <- rbind(pca_thres, data.frame(pca_thres = thres, n_PCs = n_PCs))
  }
  res_PCA_liml$pca_thres <- pca_thres

  # unconditional correlation between exposures set to 0
  if (is.null(cor.x)) cor.x <- diag(length(names_risk_factors))
  get_pca_res <- function(n_PCs) {
    PCA_liml(
      bx = bx, sx = sx,
      by = by, sy = sy,
      rho0 = ld,
      cor.x = cor.x,
      nx = nx, ny = ny,
      r = n_PCs
    )
  }

  if (ncol(bx) == 1) {
    multivariate <- F
    MRInputObj <- MendelianRandomization::mr_input(bx = as.vector(bx), bxse = as.vector(sx), by = as.vector(by), byse = as.vector(sy), correlation = ld)
  } else {
    multivariate <- T
    MRInputObj <- MendelianRandomization::mr_mvinput(bx = bx, bxse = sx, by = as.vector(by), byse = as.vector(sy), correlation = ld)
  }

  res_PCA_liml$pca_res <- res_PCA_liml$ci <- res_PCA_liml$ci.nr <- list()
  for (idx_pca in 1:nrow(pca_thres)) {
    thres <- pca_thres[idx_pca, "pca_thres"]
    n_PCs <- pca_thres[idx_pca, "n_PCs"]

    pca_res <- NULL
    n_try <- 0
    while (is.null(pca_res) & n_try <= 3) {
      ID <- paste0("PCA_thres_", thres * 100, "%_", n_PCs, "PCs")
      print(ID)
      pca_res <- tryCatch(
        {
          get_pca_res(n_PCs)
        },
        error = function(e) {
          return(NULL)
        }
      )

      if (!is.null(pca_res)) {
        res_PCA_liml$pca_res[[ID]] <- pca_res
        # robust PCA-GMM results using principal components explaining 99.99% of genetic variation
        ci <- get_table_from_pca_res(pca_res, n_risk_factors = length(names_risk_factors))
        # non-robust PCA-GMM results using principal components explaining 99.99% of genetic variation
        ci.nr <- get_table_from_pca_res(pca_res,
          n_risk_factors = length(names_risk_factors),
          non_robust = T
        )

        mr_pcgmm_Fstat <- compute_Fstat(MRInputObj, nx, ny, thres, n_PCs, multivariate)

        ci$Fstat <- unlist(mr_pcgmm_Fstat[1])
        ci.nr$Fstat <- unlist(mr_pcgmm_Fstat[2])
        rownames(ci) <- rownames(ci.nr) <- names_risk_factors
        res_PCA_liml$ci[[ID]] <- ci
        res_PCA_liml$ci.nr[[ID]] <- ci.nr
      } else {
        n_PCs <- n_PCs + 1
        n_try <- n_try + 1
      }
    }
    if (is.null(pca_res)) stop("Error in PCA_liml.")
  }

  write.table(res_PCA_liml$pca_thres, file.path(dir_output, "02_pca_thres.txt"), row.names = F, quote = F, col.names = T, sep = "\t")
  do_append <- F
  incl_col.names <- T
  for (this_type in c("ci", "ci.nr")) {
    for (group in names(res_PCA_liml[[this_type]])) {
      df <- data.frame(Type = this_type, group, res_PCA_liml[[this_type]][[group]])
      df$factor <- rownames(df)
      write.table(df,
        file.path(dir_output, "02_MVMR_PCA_liml.txt"),
        row.names = F, quote = F, col.names = incl_col.names, sep = "\t", append = do_append
      )
      do_append <- T
      incl_col.names <- F
    }
  }
  return(res_PCA_liml)
}

forestplot_main_combined <- function(name_outcome, this_window, eQTLonly = F) {
  p_df <- NULL
  lst_eQTL <- c("Brain caudate", "Small intestine", "Liver", "Pancreas", "Muscle", "Blood", "Fat", "Adrenal gland", "Testis", "Ovary", "Skin sun exposed", "Skin not sun exposed", "Adipose subcutaneous", "Adipose visceral", "Brain cerebellum", "Brain cortex", "Brain hippocampus", "Brain putamen", "Brain spinal cord", "Brain substantia nigra")
  lst_of_genes <- c("HMGCR", "PCSK9", "LPL", "ANGPTL4")
  for (gene_of_interest in lst_of_genes) {
    DIR <- paste0("results/window_", this_window, "/", gene_of_interest, "/MVMR/", name_outcome)
    lst_files <- list.files(DIR, full.names = T)
    lst_files <- lst_files[base::sapply(lst_files, function(x) file.info(x)$isdir)]

    for (this_file in lst_files) {
      this_pair <- basename(this_file)
      names_risk_factors <- unlist(strsplit(this_pair, "_"))

      if (eQTLonly) {
        if (!all(names_risk_factors %in% lst_eQTL)) next
      } else {
        if (sum(names_risk_factors %in% lst_eQTL) != 0) next
      }
      if (name_outcome %in% names_risk_factors) next
      fname <- file.path(this_file, "02_MVMR_PCA_liml.txt")
      if (!file.exists(fname)) next

      df <- read.delim(fname) %>%
        dplyr::mutate(
          gene = gene_of_interest,
          estimate = as.numeric(estimate),
          lower = as.numeric(lower),
          upper = as.numeric(upper),
          p.value = as.numeric(p.value),
          nlogP = -log10(p.value),
          Thres = factor(sapply(strsplit(group, "_"), function(x) x[3]), levels = c("99%", "99.9%", "99.99%")),
          Fstat = ifelse(Fstat > 10, 10, Fstat),
          Direction = factor(ifelse(estimate > 0, "Pos", "Neg"), levels = c("Pos", "Neg")),
          Pairs = this_pair,
          factor = factor
        ) %>%
        filter(Type == "ci", Thres == "99%")
      p_df <- rbind(p_df, df)
    }
  }
  if (is.null(p_df)) {
    return(NULL)
  }
  p_df_filt <- p_df %>%
    group_by(gene) %>%
    mutate(
      adjP = p.adjust(p.value, method = "fdr"),
      sig = p.value < 0.05,
      col = ifelse(adjP < .05, "sig1", ifelse(p.value < .05, "sig2", "notsig"))
    ) %>%
    ungroup() %>%
    group_by(gene, Pairs) %>%
    filter(sum(sig) > 0)

  if (nrow(p_df_filt) == 0) {
    return(NULL)
  }
  p1 <- ggplot(p_df_filt, aes(x = factor, y = estimate, fill = Thres, colour = col)) +
    geom_hline(yintercept = 0, lty = 2, col = "#bdb9b9") +
    geom_point(shape = 15, size = 1.25, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.5), width = 0.2, size = 0.5) +
    theme_classic() +
    ggh4x::facet_grid2(Pairs ~ gene, scales = "free") +
    ggtitle(paste0(name_outcome, "")) +
    scale_color_manual(values = c(notsig = "black", sig1 = "#aa0202", sig2 = "#0f7eba")) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none", strip.text.y = element_blank()) +
    xlab("") +
    ylab("Estimate") +
    coord_flip()
  n_pairs <- length(unique(p_df_filt$Pairs))
  this_width <- max(4.5, length(unique(p_df$p_df_filt)) * 2)
  this_height <- max(3, n_pairs / 2)
  pdf(paste0("results/window_", this_window, "/06_MVMR_PCA_liml_LM_", name_outcome, "_eQTLonly_", eQTLonly, ".pdf"), width = this_width, height = this_height)
  plot(p1)
  dev.off()
  return(p_df_filt)
}

forestplot_main_combined2 <- function(gene_of_interest, this_window, eQTLonly = F) {
  p_df <- NULL
  lst_eQTL <- c("Brain caudate", "Small intestine", "Liver", "Pancreas", "Muscle", "Blood", "Fat", "Adrenal gland", "Testis", "Ovary", "Skin sun exposed", "Skin not sun exposed", "Adipose subcutaneous", "Adipose visceral", "Brain cerebellum", "Brain cortex", "Brain hippocampus", "Brain putamen", "Brain spinal cord", "Brain substantia nigra")
  lst_of_outcomes <- c("CAD", "T2D", "Alzheimer", "Pancreatitis", "Hepatitis", "Myositis", "Myalgia", "Parkinson", "Atherosclerosis", "Amyotrophic lateral sclerosis", "Nonalcoholic fatty liver disease")
  for (name_outcome in lst_of_outcomes) {
    DIR <- paste0("results/window_", this_window, "/", gene_of_interest, "/MVMR/", name_outcome)
    lst_files <- list.files(DIR, full.names = T)
    lst_files <- lst_files[base::sapply(lst_files, function(x) file.info(x)$isdir)]

    for (this_file in lst_files) {
      this_pair <- basename(this_file)
      names_risk_factors <- unlist(strsplit(this_pair, "_"))

      if (eQTLonly) {
        if (!all(names_risk_factors %in% lst_eQTL)) next
      } else {
        if (sum(names_risk_factors %in% lst_eQTL) != 0) next
      }
      if (name_outcome %in% names_risk_factors) next
      fname <- file.path(this_file, "02_MVMR_PCA_liml.txt")
      if (!file.exists(fname)) next

      df <- read.delim(fname) %>%
        dplyr::mutate(
          gene = gene_of_interest,
          outcome = name_outcome,
          estimate = as.numeric(estimate),
          lower = as.numeric(lower),
          upper = as.numeric(upper),
          p.value = as.numeric(p.value),
          nlogP = -log10(p.value),
          Thres = factor(sapply(strsplit(group, "_"), function(x) x[3]), levels = c("99%", "99.9%", "99.99%")),
          Fstat = ifelse(Fstat > 10, 10, Fstat),
          Direction = factor(ifelse(estimate > 0, "Pos", "Neg"), levels = c("Pos", "Neg")),
          Pairs = this_pair,
          factor = factor
        ) %>%
        filter(Type == "ci", Thres == "99%")
      p_df <- rbind(p_df, df)
    }
  }
  if (is.null(p_df)) {
    return(NULL)
  }
  p_df_filt <- p_df %>%
    mutate(
      sig = p.value < 0.05,
      col = ifelse(p.value < .05, "sig", "notsig"),
      outcome = sapply(strsplit(p_df$outcome, " "), function(x) paste(x, collapse = "\n"))
    ) %>%
    group_by(outcome, Pairs) %>%
    filter(sum(sig) > 0) %>%
    group_by(Pairs) %>%
    filter(length(Pairs) > 2)

  if (nrow(p_df_filt) == 0) {
    return(NULL)
  }
  p1 <- ggplot(p_df_filt, aes(x = factor, y = estimate, colour = col)) +
    geom_hline(yintercept = 0, lty = 2, col = "#bdb9b9") +
    geom_point(shape = 15, size = 1.25, position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.5), width = 0.2, size = 0.5) +
    theme_classic() +
    ggh4x::facet_grid2(Pairs ~ outcome, scales = "free") +
    ggtitle(paste0(gene_of_interest, "")) +
    scale_color_manual(values = c(notsig = "black", sig = "#aa0202")) +
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none", strip.text.y = element_blank()) +
    xlab("") +
    ylab("Estimate") +
    coord_flip()
  n_pairs <- length(unique(p_df_filt$Pairs))
  this_width <- max(4.5, length(unique(p_df_filt$outcome)) * 2)
  this_height <- max(3, n_pairs / 2)
  pdf(paste0("results/window_", this_window, "/08_MVMR_PCA_liml_LM_", gene_of_interest, "_eQTLonly_", eQTLonly, ".pdf"), width = this_width, height = this_height)
  plot(p1)
  dev.off()
  return(p_df_filt)
}

plot_propcoloc_barplots_pairwise <- function(this_gene, lst_pairs, this_window, eQTLonly) {
  DIR <- paste0("results/window_", this_window, "/", this_gene, "/propcoloc")
  lst_factors <- unique(unlist(strsplit(lst_pairs, "_")))

  df_res <- NULL
  for (fac1 in lst_factors) {
    for (fac2 in lst_factors) {
      this_group <- paste0(fac1, "_", fac2)
      this_dir <- paste0(DIR, "/", this_group)
      if (!dir.exists(this_dir)) {
        this_group <- paste0(fac2, "_", fac1)
        this_dir <- file.path(DIR, this_group)
        if (!dir.exists(this_dir)) next
      }
      fname1 <- paste0(this_dir, "/03_prop.coloc_", fac1, "_", fac2, ".RDS")
      fname2 <- paste0(this_dir, "/03_prop.coloc_", fac2, "_", fac1, ".RDS")
      if (file.exists(fname1)) {
        fac1_ <- paste(unlist(strsplit(fac1, " ")), collapse = "\n")
        fac2_ <- paste(unlist(strsplit(fac2, " ")), collapse = "\n")
        df1 <- with(readRDS(fname1), data.frame(group = paste0(fac1, "_", fac2), fac1 = fac1_, fac2 = fac2_, p_full, p_cond, LM_full, LM_cond))
        df2 <- with(readRDS(fname2), data.frame(group = paste0(fac2, "_", fac1), fac1 = fac2_, fac2 = fac1_, p_full, p_cond, LM_full, LM_cond))
        df_res <- rbind(df_res, rbind(df1, df2))
      }
    }
  }
  df_res <- df_res %>%
    mutate(
      `Prop-coloc-cond P` = ifelse(-log10(p_cond) > 10, 10, -log10(p_cond)),
      `LM-cond P` = ifelse(-log10(LM_cond) > 10, 10, -log10(LM_cond)),
      sig = p_cond < .05 & LM_cond < .05
    )
  group_selected <- df_res %>%
    group_by(group) %>%
    summarise(selected = all(sig)) %>%
    filter(selected) %>%
    pull(group)

  p_df <- df_res %>%
    select(group, fac1, fac2, "Prop-coloc-cond P", "LM-cond P") %>%
    reshape2::melt() %>%
    mutate(Color = factor(ifelse(value < -log10(0.05), "P > 0.05", as.character(variable)), levels = c("Prop-coloc-cond P", "LM-cond P", "P > 0.05")))

  p <- easyGgplot2::ggplot2.barplot(
    data = p_df,
    xName = "variable", yName = "value",
    facetingVarNames = c("fac1", "fac2"),
    groupName = "Color",
    faceting = TRUE
  ) +
    scale_fill_manual(values = c("#00c5c5", "#19bb04", "#888888")) +
    theme_bw() +
    geom_hline(yintercept = -log10(0.05), color = "red", lty = 2) +
    ggtitle(this_gene) + ylab("-log10(P)") + xlab("") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    geom_rect(
      data = (p_df %>% filter(group %in% group_selected)),
      fill = NA, colour = "red", xmin = -Inf, xmax = Inf,
      lwd = 2, ymin = -Inf, ymax = Inf
    )
  this_width <- max(5, length(lst_factors) * 1.6)
  this_height <- max(5, length(lst_factors) * 1.6) - 1.5
  pdf(paste0("results/window_", this_window, "/07_plot_pairs_", this_gene, "_eQTLonly_", eQTLonly, ".pdf"), width = this_width, height = this_height)
  plot(p)
  dev.off()
}



plot_associations <- function(this_gene, lst_pairs, this_window, eQTLonly) {
  DIR <- paste0("results/window_", this_window, "/", this_gene, "/propcoloc")
  if (length(lst_pairs) == 0) next
  this_width <- 4.5
  this_height <- 3.5
  pdf(paste0("results/window_", this_window, "/09_plot_associations_", this_gene, "_eQTLonly_", eQTLonly, ".pdf"), width = this_width, height = this_height)
  for (this_pair in lst_pairs) {
    this_dir <- paste0(DIR, "/", this_pair)
    fac1 <- strsplit(this_pair, "_")[[1]][1]
    fac2 <- strsplit(this_pair, "_")[[1]][2]

    fname1 <- paste0(this_dir, "/03_prop.coloc_", fac1, "_", fac2, ".RDS")
    fname2 <- paste0(this_dir, "/03_prop.coloc_", fac2, "_", fac1, ".RDS")
    res1 <- readRDS(fname1)
    res2 <- readRDS(fname2)
    plot(res1$fig_uni)
    plot(res1$fig_multi)
    plot(res2$fig_uni)
    plot(res2$fig_multi)
  }
  dev.off()
}
