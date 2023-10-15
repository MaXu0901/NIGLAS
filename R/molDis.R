#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dom_comm PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION, Default: 'euclidean'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[vegan]{vegdist}}
#' @rdname molDis
#' @export 
#' @importFrom vegan vegdist
molDis <- function(dom_comm,method = "euclidean") {
  x <- vegan::vegdist(dom_comm,method = "bray")
  return(x)
}


