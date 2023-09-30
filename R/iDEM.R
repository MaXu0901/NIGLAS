#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param comm_dom PARAM_DESCRIPTION
#' @param know_seq PARAM_DESCRIPTION
#' @param n PARAM_DESCRIPTION, Default: length(know_seq)
#' @param bootstrap PARAM_DESCRIPTION, Default: 100
#' @param threshold PARAM_DESCRIPTION, Default: 0.3
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname iDEM
#' @export 
iDEM <- function(comm_dom,know_seq,n=length(know_seq),bootstrap = 100,threshold=0.3) {
  dark_seq <- colnames(comm_dom)[!colnames(comm_dom) %in% know_seq]
  know_sel <- sample(know_seq,n)
  comm_know <- comm_dom[,know_sel]
  comm_dk <- comm_dom[,c(sample(know_sel,0.5*n),sample(dark_seq,0.5*n))]
  list_net_know <- list()
  list_net_dk <- list()
  for (i in 1:bootstrap) {
    net_know <- SpiecEasi::sparcc(comm_know)
    net_dk <- SpiecEasi::sparcc(comm_dk)
    assign(paste("net_know",i),net_know)
    assign(paste("net_dk",i),net_dk)
    list_net_know[[i]] <- get(paste("net_know",i))
    list_net_dk[[i]] <- get(paste("net_dk",i))
  }
  NET_node.metric <- function(sparcc.r,threshold){
    sparcc.r.all <- ifelse(abs(sparcc.r) > threshold, sparcc.r, 0)
    g.all <- igraph::graph_from_adjacency_matrix(as.matrix(sparcc.r.all), mode = 'undirected', weighted = TRUE, diag = FALSE)
    node.all <- data.frame(
      nodes_id = V(g.all)$name,Degree = degree(g.all)
      )
    return(mean(node.all$Degree))
  }
  degree_net_know <- lapply(list_net_know, function(net) NET_node.metric(net, threshold))
  degree_net_dk <- lapply(list_net_dk, function(net) NET_node.metric(net, threshold))
  mean_degree_net_know <- mean(unlist(degree_net_know))
  mean_degree_net_dk <- mean(unlist(degree_net_dk))
  return((mean_degree_net_dk/mean_degree_net_know-1)*100)
}
