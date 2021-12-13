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

  R_query = format_BM_from_XML(QUERY)

  if (! is.null(wbps_gene_id)) {
    BM_addFilters(R_query, wbps_gene_id = wbps_gene_id)
  }
  else if (! is.null(gene_list)) {
    BM_addFilters(R_query, gene_list = gene_list)
  }

  genes = runWithMart(R_query, mart)
  genes$chromosome_name = paste("chr",genes$chromosome_name, sep='')
  genes$strand = ifelse(genes$strand==1, '+', '-')
  GenomicRanges::makeGRangesFromMartDataFrame(genes)
}

getCElegansPromoters = function(mart, upstream=2000,downstream=500,gene_list=NULL, WBID=NULL) {
  genes_df = getCElegansGeneLocs(mart, gene_list=NULL, WBID=NULL)
  GenomicRanges::promoters(genes_df, upstream, downstream)
}
