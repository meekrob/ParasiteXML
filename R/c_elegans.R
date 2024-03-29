#' Return C. elegans protein coding genes as a Genomic Ranges object.
#'
#' Run a preformatted query to get C. elegans protein coding gene locations, with options to subset the
#' query by gene name or wormbase gene IDs.
#' @param mart A 'mart' object as returned by 'getParamart'.
#' @param gene_list Retrieve only genomic locations for gene names specified here, rather than the whole genome.
#' @param WBID Retrieve only genomic locations for Wormbase IDs specified here, rather than the whole genome.
#' @return A Genomic Ranges object.
#' @seealso \code{\link{getCElegansPromoters}}, \code{\link{getParamart}}, \code{\link[GenomicRanges]{GenomicRanges}}
#' @export
getCElegansGeneLocs = function(mart, gene_list=NULL, WBID=NULL) {
  QUERY='<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "biotype" value = "protein_coding"/>
		<Filter name = "species_id_1010" value = "caelegprjna13758"/>
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "external_gene_id" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "strand" />
	</Dataset>
</Query>'
  library(dplyr)
  R_query = format_BM_from_XML(QUERY)

  if (! is.null(WBID)) {
    R_query = BM_addFilters(R_query, wbps_gene_id = WBID)
  }
  else if (! is.null(gene_list)) {
    R_query = BM_addFilters(R_query, gene_name = gene_list)
  }

  genes = runWithMart(R_query, mart)
  genes$chromosome_name = paste("chr",genes$chromosome_name, sep='')
  genes %>% 
    mutate(chromosome_name = replace(chromosome_name, chromosome_name == "chrMtDNA", "chrM")) -> genes
  genes$strand = ifelse(genes$strand==1, '+', '-')
  gr = makeGRangesFromMartDataFrame(genes)
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::Seqinfo(genome='ce11')
  gr
}
#' Return promoter ranges based on arguments, applied to the result of getCElegansGeneLocs.
#' 
#' 
#' @export
getCElegansPromoters = function(mart, upstream=2000,downstream=500,gene_list=NULL, WBID=NULL) {
  genes_df = getCElegansGeneLocs(mart, gene_list=NULL, WBID=NULL)
  GenomicRanges::promoters(genes_df, upstream, downstream)
}
