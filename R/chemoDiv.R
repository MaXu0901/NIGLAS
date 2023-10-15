#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param chem_dataset PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' #EXAMPLE1
#' # Example data of a DOM compositional table (50 samples by 200 DOM molecules)
#' data(data)
#' # Calculate selected types of α-diversity and evenness measures
#' chemoDiv(data)
#' }
#' @rdname chemoDiv
#' @export 
chemoDiv <- function(chem_dataset) {
  library(vegan)
  data.frame(Sample = rownames(chem_dataset),
             shannon = diversity(chem_dataset,index = "shannon"),
             simpson = diversity(chem_dataset,index = "simpson"),
             invsimpson = diversity(chem_dataset,index = "invsimpson"))
}
