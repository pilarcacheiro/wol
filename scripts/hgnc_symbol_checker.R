###############################################################################################################################
###############################################################################################################################
### Script: hgnc_symbol_checker.R 
### Purpose: check gene symbols and retrieve hgnc id 
### Description: this function returns hgnc ids for protein coding genes symbols / (only protein coding genes /input file could 
### be modified)
### Input: "gene_with_protein_product.txt" (see README) file and a vector of gene symbols to check 
### Output: dataframe with 3 columns: "HGNC.ID":corresponding hgnc id ("-" , if no hgnc id was found): "Gene.Symbol": gene 
### symbol provided; "Type": mapping type (Approved.Symbol,Synonym.Symbol,Notfound.ProteinCoding.Symbol,...) 
###############################################################################################################################
###############################################################################################################################


hgnc.checker <- function(gene.symbols,gene.file){
  
  library(tidyr); library(dplyr); library(data.table)
  
  check.approved <- function(input.genes,database) {
    
    return(database %>% dplyr::select(hgnc_id,symbol) %>% mutate_if(is.factor,as.character) %>%
             filter(symbol!="") %>%
             filter(symbol %in% input.genes) %>% dplyr::rename (Gene.Symbol = symbol,HGNC.ID = hgnc_id) %>% 
             mutate_if(is.factor,as.character) %>% mutate(Type = "Approved.Symbol"))
    
  }
  
  
  check.synonyms <- function(input.genes,database){
    
    
    return(database %>% dplyr::select(hgnc_id,alias_symbol) %>% mutate_if(is.factor,as.character) %>%
             filter(alias_symbol!="") %>%
             tidyr::separate_rows(alias_symbol,sep="\\|") %>% dplyr::rename(HGNC.ID = hgnc_id) %>%
             mutate(Gene.Symbol = trimws(alias_symbol)) %>%  
             dplyr::select(HGNC.ID,Gene.Symbol)  %>% filter(Gene.Symbol %in% input.genes) %>%
             mutate_if(is.factor,as.character) %>% mutate(Type = "Synonym.Symbol"))
    
  }
  
  check.previous <- function(input.genes,database) {
    
    return(database %>% dplyr::select(hgnc_id,prev_symbol) %>% mutate_if(is.factor,as.character) %>%
             filter(prev_symbol!="") %>%
             tidyr::separate_rows(prev_symbol,sep="\\|") %>% dplyr::rename(HGNC.ID = hgnc_id) %>%
             mutate(Gene.Symbol = trimws(prev_symbol)) %>%  
             dplyr::select(HGNC.ID,Gene.Symbol)  %>% filter(Gene.Symbol %in% input.genes) %>%
             mutate_if(is.factor,as.character) %>% mutate(Type = "Previous.Symbol"))
    
  }
  
  
  
  check.duplicates.symbol <- function(file.to.check.symbols,duplicates.symbol){
    
    if(!length(duplicates.symbol)) { return(file.to.check.symbols)}
    
    else{ 
      
      final.nodup.symbol <- file.to.check.symbols %>% filter(!Gene.Symbol %in% duplicates.symbol)
      duplicate.symbols.df <- data.frame(HGNC.ID = rep("-",length(duplicates.symbol)), Gene.Symbol = duplicates.symbol,
                                         Type = "Ambiguous.Symbol")
      final.dups.symbol <- rbind(final.nodup.symbol,duplicate.symbols.df)
      
      return(final.dups.symbol)
      
    }
  }
  
  
  check.duplicates.id <- function(file.to.check.ids,duplicates.id){
    
    if(!length(duplicates.id)) { return(file.to.check.ids)}
    
    else{ 
      
      final.nodup.id <- file.to.check.ids %>% filter(!HGNC.ID %in% duplicates.id)
      duplicate.ids.df <- file.to.check.ids %>% filter(HGNC.ID %in% duplicates.id) %>% mutate(HGNC.ID = "-",Type="Ambiguous.Symbol")
      final.dups.id <- rbind(final.nodup.id,duplicate.ids.df)
      
      return(final.dups.id)
      
    }
  }
  
  
  genes <- trimws(gene.symbols)
  
  hgnc <- gene.file
  
  hgnc.approved <- check.approved(genes,hgnc)
  
  hgnc.synonyms <- check.synonyms(genes[! genes %in% hgnc.approved$Gene.Symbol],hgnc)
  
  hgnc.previous <- check.previous(genes[! genes %in% c(hgnc.approved$Gene.Symbol,hgnc.synonyms$Gene.Symbol)],hgnc)
  
  genes.not.found <- genes[! genes %in% c(hgnc.approved$Gene.Symbol,hgnc.synonyms$Gene.Symbol,hgnc.previous$Gene.Symbol)]
  
  hgnc.notfound <- data.frame(HGNC.ID = rep("-",length(genes.not.found)),Gene.Symbol = genes.not.found) %>%  mutate_if(is.factor,as.character) %>%
    mutate(Type = "Notfound.ProteinCoding.Symbol")
  
  hgnc.all <- rbind(hgnc.approved,hgnc.synonyms,hgnc.previous,hgnc.notfound)
  
  duplicates.symbol <- hgnc.all$Gene.Symbol[duplicated(hgnc.all$Gene.Symbol)]
  
  results.noduplicated.symbol <- check.duplicates.symbol(hgnc.all,duplicates.symbol)
  
  duplicates.id <- results.noduplicated.symbol$HGNC.ID[duplicated(results.noduplicated.symbol$HGNC.ID)& results.noduplicated.symbol$HGNC.ID!="-"]
  
  results.noduplicated.id <- check.duplicates.id(results.noduplicated.symbol,duplicates.id)
  
  results.final <- results.noduplicated.id
  
  return(results.final)
  
} 