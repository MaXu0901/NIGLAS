#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dom_comm PARAM_DESCRIPTION
#' @param dom_trait PARAM_DESCRIPTION
#' @param metric PARAM_DESCRIPTION, Default: c("MolForm")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{case_when}}
#' @rdname molDK
#' @export 
#' @importFrom dplyr select mutate case_when
molDK <- function(dom_comm,dom_trait,metric = c("MolForm")){
  
  if(identical(colnames(dom_comm),rownames(dom_trait)) == F){
    stop("Species labels in 'mol' and 'trait' need to be identical and ordered alphabetically (or simply in the same order).", "\n")
  }else{
    mol.Form = dom_trait %>% 
      dplyr::select(all_of(metric)) %>% 
      dplyr::mutate(Dark_or_Known = dplyr::case_when(is.na(.)~"Dark",
                                                             !is.na(.)~"Known"))
  }
  return(mol.Form)
}

