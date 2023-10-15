#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dom_comm PARAM_DESCRIPTION
#' @param MARGIN PARAM_DESCRIPTION, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[vegan]{decostand}}
#' @rdname molRA
#' @export 
#' @importFrom vegan decostand
molRA <- function(dom_comm,MARGIN = 1) {
  x <- vegan::decostand(dom_comm,MARGIN = MARGIN,method = "total")
  return(x)
}



