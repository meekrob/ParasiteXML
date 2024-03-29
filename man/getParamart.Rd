% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_BM_from_XML.R
\name{getParamart}
\alias{getParamart}
\title{Return a database object initialized with paraSite (Wormbase)}
\usage{
getParamart(
  biomart = "parasite_mart",
  dataset = "wbps_gene",
  host = "https://parasite.wormbase.org",
  port = 443,
  ...
)
}
\value{
A Mart object from package biomaRt.
}
\description{
This function returns a "mart" object with paraSite options as default.
See biomaRt::useMart for alternative marts and the required arguments for them.
}
\examples{
query = format_BM_from_XML('<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
  <Dataset name = "wbps_gene" interface = "default" >
     <Filter name = "biotype" value = "lincRNA"/>
     <Filter name = "species_id_1010" value = "caelegprjna13758"/>
     <Attribute name = "wbps_gene_id" />
  </Dataset>
</Query>')

library(biomaRt)
mart = getParamart()
lincs = runWithMart(query, mart)
}
