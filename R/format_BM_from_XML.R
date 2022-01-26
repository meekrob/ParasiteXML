# 1. Go to: https://parasite.wormbase.org/biomart
# 2. Configure a query and
# <?xml version="1.0" encoding="UTF-8"?>
# <!DOCTYPE Query>
# <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
#   <Dataset name = "wbps_gene" interface = "default" >
#     <Filter name = "only_hsapiens_homologue" excluded = "0"/>
#     <Filter name = "gene_name" value = "ELT-2,ELT-7"/>
#     <Filter name = "species_id_1010" value = "caelegprjna13758"/>
#     	<Attribute name = "production_name_1010" />
#     	<Attribute name = "wbps_gene_id" />
#     	<Attribute name = "hsapiens_gene" />
#     	<Attribute name = "hsapiens_gene_name" />
#   </Dataset>
# </Query>'

#' @export
AllCElegansProteinCoding = '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "biotype" value = "protein_coding"/>
		<Filter name = "species_id_1010" value = "caelegprjna13758"/>
		<Attribute name = "wbps_gene_id" />
	</Dataset>
</Query>'

#' @export
AllCElegansLINC = '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "biotype" value = "lincRNA"/>
		<Filter name = "species_id_1010" value = "caelegprjna13758"/>
		<Attribute name = "wbps_gene_id" />
	</Dataset>
</Query>'

#' @export
AllCElegansPseudogene = '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "biotype" value = "pseudogene"/>
		<Filter name = "species_id_1010" value = "caelegprjna13758"/>
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "external_gene_id" />
	</Dataset>
</Query>'

#' @export
IDConversionQuery = '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "gene_name" value = "ELT-2,ELT-7"/>
		<Filter name = "species_id_1010" value = "caelegprjna13758"/>
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "external_gene_id" />
	</Dataset>
</Query>'

#' Translate XML Query to R Code
#'
#' Take a text object returned from the XML formatted output at
#' \url{https://parasite.wormbase.org/biomart} and return a function call
#' to biomaRt::getBM with the filter names, values, and requested output
#' fields.
#' @param xmltext XML text copied from \url{https://parasite.wormbase.org/biomart}
#' @param printCode TRUE/FALSE Print the R code to the console in addition to returning the call.
#' @return A 'call' object that can be run with runWithMart().
#' @seealso \code{\link[biomaRt]{getBM}}, \code{\link{runWithMart}}
#' @examples
#' # A query that extracts all C. elegans long noncoding RNAs.
#' query = format_BM_from_XML('<?xml version="1.0" encoding="UTF-8"?>
#' <!DOCTYPE Query>
#' <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
#'   <Dataset name = "wbps_gene" interface = "default" >
#'      <Filter name = "biotype" value = "lincRNA"/>
#'      <Filter name = "species_id_1010" value = "caelegprjna13758"/>
#'      <Attribute name = "wbps_gene_id" />
#'   </Dataset>
#' </Query>')
#'
#' mart = getParamart()
#' lincs = runWithMart(query, mart)
#' head(lincs)
#' @export
format_BM_from_XML = function(xmltext, printCode=TRUE) {

  query=xml2::xml_find_all(xml2::read_xml(xmltext), "/Query")
  query_dataset = xml2::xml_find_first(query, "Dataset")
  filter_nodes = xml2::xml_find_all(query_dataset, "Filter")
  attribute_nodes = xml2::xml_find_all(query_dataset, "Attribute")

  filter_specs = split(xml2::xml_attr(filter_nodes, "value"),
                       xml2::xml_attr(filter_nodes, "name"))

  # handle booleans
  filter_specs[is.na(filter_specs)]<- T

  output_specs = xml2::xml_attr(attribute_nodes, "name")

  cl = call("biomaRt::getBM", filter=names(filter_specs), value=filter_specs, attributes=output_specs)
  if (printCode) {print(cl)}
  invisible(cl)
}

#' Add Query Attributes To biomaRt Call
#'
#' Add a set of attributes (output fields) requested by the query.
#' @param cll biomaRt call as returned by format_BM_from_XML()
#' @param atts String vector of names of attributes to add.
#' @return The updated biomaRt call object
#' @export
BM_addAttributes = function(cll, atts)
{
  cll$attributes = c(cll$attributes,atts)
  cll
}

#' Remove Query Attributes From biomaRt Call
#'
#' Remove a set of attributes (output fields) requested by the query.
#' @param cll biomaRt call as returned by format_BM_from_XML()
#' @param atts String vector of names of attributes to remove.
#' @return The updated biomaRt call object
#' @export
#'
BM_removeAttributes = function(cll, atts)
{
  matches = which(cll$attributes %in% atts)
  cll$attributes = cll$attributes[-matches]
  cll
}

#' Remove Filters From a biomaRt Call.
#'
#' Remove filters and their values from the biomaRt call.
#' @return The updated biomaRt call object
#' @seealso \code{\link{BM_addFilters}}, \code{\link{BM_substituteFilterValues}}
#' @export
BM_removeFilters = function(cll, removefilter)
{
  matches = which(cll$filter %in% removefilter)
  cll$filter = cll$filter[-matches]
  cll$value[removefilter] <- NULL
  cll
}


#' Add filters to a biomaRt call.
#'
#' Add a set of filters to a biomaRt call.
#' @param cll biomaRt call as returned by format_BM_from_XML()
#' @param ... key=value arguments to add to the biomaRt call.
#' @return The updated biomaRt call object
#' @seealso \code{\link{BM_removeFilters}}, \code{\link{BM_substituteFilterValues}}
#' @export
BM_addFilters = function(cll, ...)
{
  filterList = list(...)
  # The syntax of filter, value is different for the call than for the input
  # to 
  cll$filter = c(cll$filter, names(filterList))
  cll$value = c(cll$value, filterList)
  cll
}

#' Override filters in the biomaRt call.
#'
#' Replace filter values in a biomaRt call.
#' @param cll biomaRt call as returned by format_BM_from_XML()
#' @param ... key=value arguments to replace in the biomaRt call.
#' @return The updated biomaRt call object
#' @seealso \code{\link{BM_removeFilters}}, \code{\link{BM_addFilters}}
#' @export
BM_substituteFilterValues = function(cll, ...)
{
  overrides = list(...)
  for (name in names(overrides)) {
    if (name %in% cll$filter)
    {
      cll$value[name] = overrides[name]
    }
  }
  cll
}

#' Run a pre-formatted query with the provided "Mart" object.
#'
#' Run a pre-formatted query, as returned by format_BM_from_XML() to biomaRt
#' with a "Mart" object as returned by getParamart() or useMart().
#' @param cll biomaRt call as returned by format_BM_from_XML()
#' @param mart a "Mart" object as returned by getParamart() or useMart().
#' @return The query result as a dataframe.
#'
#' @seealso \code{\link[biomaRt]{getBM}}, \code{\link{getParamart}}
#' @examples
#' query = format_BM_from_XML('<?xml version="1.0" encoding="UTF-8"?>
#' <!DOCTYPE Query>
#' <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
#'   <Dataset name = "wbps_gene" interface = "default" >
#'      <Filter name = "biotype" value = "lincRNA"/>
#'      <Filter name = "species_id_1010" value = "caelegprjna13758"/>
#'      <Attribute name = "wbps_gene_id" />
#'   </Dataset>
#' </Query>')
#'
#' library(biomaRt)
#' mart = getParamart()
#' lincs = runWithMart(query, mart)
#' @export
runWithMart = function(cll, mart)
{
  cll$mart = mart
  output=eval(cll)
  cll$mart = NULL
  output
}

#' Return a database object initialized with paraSite (Wormbase)
#'
#' This function returns a "mart" object with paraSite options as default.
#' See biomaRt::useMart for alternative marts and the required arguments for them.
#' @return A Mart object from package biomaRt.
#' @examples
#' query = format_BM_from_XML('<?xml version="1.0" encoding="UTF-8"?>
#' <!DOCTYPE Query>
#' <Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
#'   <Dataset name = "wbps_gene" interface = "default" >
#'      <Filter name = "biotype" value = "lincRNA"/>
#'      <Filter name = "species_id_1010" value = "caelegprjna13758"/>
#'      <Attribute name = "wbps_gene_id" />
#'   </Dataset>
#' </Query>')
#'
#' library(biomaRt)
#' mart = getParamart()
#' lincs = runWithMart(query, mart)
#' @export
getParamart = function(biomart="parasite_mart", dataset="wbps_gene", host="https://parasite.wormbase.org", port = 443,...) {
  paramart <-
    biomaRt::useMart(biomart,
            dataset = dataset,
            host = host,
            port = port,
            ...)

  if (! is.null(paramart))
  {
    library(crayon)
    cat(green("Database connected ✓"),"\n")
    cat(silver("biomart      ...      " %+% (paramart@biomart)),"\n")
    cat(silver("host         ...      " %+% (paramart@host)),"\n")
    cat(silver("dataset      ...      " %+% (paramart@dataset)),"\n")
  }
  
  paramart
}

makeGRangesFromMartDataFrame = function(df,...) {
  GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = T,
                           seqnames.field = "chromosome_name",
                           start.field = "start_position",
                           end.field = "end_position",
                           ...)
}
