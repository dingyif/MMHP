clusterAnalysis <- function(input_matrix){
  if(sum(input_matrix)==0){
    return(0)
  }
  
  g_community <- graph_from_adjacency_matrix(input_matrix,mode="undirected")
  cluster_temp <- handelExceptionCluster(g_community)
  membership <- cluster_temp$membership #cluster_leading_eigen
  number_cluster <- length(unique(membership))
  cluster_member_count <- numeric(number_cluster)
  
  within_count <- 0
  between_count <- 0
  for(cluster in c(1:number_cluster)){
    within_count <- within_count + sum(input_matrix[membership==cluster,membership==cluster])
    cluster_member_count[cluster] <- sum(membership==cluster)
    for(another_cluster in c(1:number_cluster)[-cluster]){
      between_count <- between_count + sum(input_matrix[membership==cluster,membership==another_cluster])
    }
  }
  
  between_count_normalizor <- 0
  if(number_cluster == 1){
    between_count_normalizor <- 0
  }else if(number_cluster == 2){
    between_count_normalizor <- cluster_member_count[1]*cluster_member_count[2]
  }else{
    for(cluster in c(1:(number_cluster-1))){
      for(cluster2 in c((cluster+1):number_cluster)){
        between_count_normalizor <- between_count_normalizor + cluster_member_count[cluster]*cluster_member_count[cluster2]
      }
    }
  }
  
  normalized_between_count <- between_count/between_count_normalizor
  normalized_within_count <- within_count/ifelse(sum(sapply(cluster_member_count,choose,k=2))==0,1,sum(sapply(cluster_member_count,choose,k=2)))
  normalized_all_count <- (within_count+between_count)/(between_count_normalizor+sum(sapply(cluster_member_count,choose,k=2)))
  return(normalized_between_count/normalized_within_count)
}

handelExceptionCluster <- function(g_community, max_iter = 10){
  return_result <- try(log("error"), silent = TRUE)
  i <- 1
  while(class(return_result) == "try-error" & i<max_iter){
    return_result <- try(cluster_leading_eigen(g_community))
    i <- i+1
  }
  return(return_result)
}
