laplacianEigen <- function(start, end, state="social"){
  if(state == "social"){
    if(start == 0){
      input_matrix <- adjm[round(end),,]-utility_state_cumcount[round(end),,]
    }else{
      input_matrix <- adjm[round(end),,]-utility_state_cumcount[round(end),,]-adjm[round(start),,]+utility_state_cumcount[round(start),,]
    }
  }else if(state == "utility"){
    if(start == 0){
      input_matrix <- utility_state_cumcount[round(end),,]
    }else{
      input_matrix <- utility_state_cumcount[round(end),,]-utility_state_cumcount[round(start),,]
    }
  }else if(state == "total"){
    if(start == 0){
      input_matrix <- adjm[round(end),,]
    }else{
      input_matrix <- adjm[round(end),,]-adjm[round(start),,]
    }
  }else{
    return("ERROR")
  }
  
  if(sum(input_matrix)==0){
    return(0)
  }
  
  g_community <- graph_from_adjacency_matrix(input_matrix,mode="undirected",weighted=TRUE)
  eigen_values <- eigen(laplacian_matrix(g_community))$values
  return(list(max_eigen = eigen_values[1],
              second_smallest_eigen = eigen_values[nrow(input_matrix)-1],
              zero_eigen = sum(eigen_values==0),
              smallest_nonzero = eigen_values[nrow(input_matrix)-sum(eigen_values==0)]))
}
