#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param chem_dataset PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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
