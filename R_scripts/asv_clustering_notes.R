load("output/physeqlists.Rdata")
BiocManager::install("DECIPHER")

devtools::install_github("https://github.com/alexpiper/seqateurs/")

devtools::install_github("https://github.com/RIVM-IIV-Microbiome/biomeUtils/")
library(biomeUtils)
biomeUtils::clus(physeqlists$MAPP_16s)


library(seqateurs)


test <- physeqlists$MAPP_23s@otu_table
test2 <- physeqlists$MAPP_23s@tax_table
test2 <- as.data.frame(test2)
test2$sequence

test <- as.matrix(test)
colnames(test) <- test2$sequence
test <- as.matrix(test)
what <- cluster_otus(test,similarity=0.97, cores=1)
refseq(physeqlists$MAPP_23s)



test <- physeqlists$MAPP_16s@otu_table
test2 <- physeqlists$MAPP_16s@tax_table
test2 <- as.data.frame(test2)

test <- as.matrix(test)
colnames(test) <- test2$sequence
test <- as.matrix(test)
what <- cluster_otus(test,similarity=0.97, cores=1)


nrow(meco_16s$tax_table)
nrow(what)
length(unique(what$cluster))

cluster_otus <- function(x, similarity=0.97, cores=1) {
    
    if(is(x, "matrix")| is(x, "data.frame")){
        asv_sequences <- colnames(x)
    } else if(is(x, "phyloseq") & !is.null(phyloseq::refseq(x, errorIfNULL = FALSE))){
        asv_sequences <- as.vector(phyloseq::refseq(x))
    } else if(is(x, "phyloseq") & is.null(phyloseq::refseq(x, errorIfNULL = FALSE))){
        message("refseq() not found. Using tax_table rownames for sequences.")
        if (sum(grepl("[^ACTG]", rownames(phyloseq::tax_table(x)))) > 0) {
            stop("Error: Taxa do not appear to be DNA sequences.")
        }
        asv_sequences <- colnames(phyloseq::get_taxa(x))
    } else if(is(x, "DNAStringSet") ){
        asv_sequences <- as.character(x)
    } else if(is(x, "DNAbin")){
        asv_sequences <- taxreturn::DNAbin2char(x)
    } else{
        stop("Error: Taxa do not appear to be DNA sequences.")
    }
    seqs <- Biostrings::DNAStringSet(asv_sequences)
    
    # define cutoffs for clustering
    if(!dplyr::between(similarity, 0, 1)){
        stop("similarity must be a number between 0 and 1")
    }
    cutoff <- 1 - similarity
    
    ## Find clusters of ASVs to form the new OTUs
    otus <- DECIPHER::Clusterize(
        seqs,
        cutoff = cutoff,  # use `cutoff = 0.03` for a 97% OTU
        method = "overlap",
        includeTerminalGaps = FALSE,
        penalizeGapLetterMatches = NA,
        minCoverage = 0.5,
        maxAlignments = 100,
        invertCenters = FALSE,
        processors = cores) %>%
        dplyr::mutate(sequence = asv_sequences)  %>%
        dplyr::group_by(cluster) %>%
        dplyr::mutate(cluster_size = n_distinct(sequence)) %>%
        dplyr::ungroup()
    
    # Return cluster memberships
    return(otus)
}
