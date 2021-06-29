parseSlurmNodeList<- function(arg1) {
  #parses SLURM_NODELIST environment variable into a list of nodes (i.e. to be passed to makePSOCKcluster() or h2o.init())
  #only works for a list from a single partition
  #assumes node indices are 4 digits long (i.e. cn-0004, not cn-004 or cn-00004)
  
  #if only 1 node, no parsing needed
  if (!grepl("(\\[)", arg1)) {
    return(arg1)
  }
  
  ##get node prefix
  node.prefex <- sub("\\[.*", "", arg1)
  
  #remove the opening square bracket and everything before it
  nodes <- sub(".*\\[", "", arg1)
  
  #remove trailing square bracket
  nodes <- sub("]", "", nodes)
  
  #split list
  nodes <- unlist(strsplit(nodes, ","))
  
  ##expand ranges denoted by "-"
  which.collapsed <- grep("-", nodes)
  
  for (i in which.collapsed) {
    collapsed <- nodes[i]
    range.bounds <- as.numeric(unlist(strsplit(collapsed, "-")))
    expanded <- paste0(range.bounds[1]:range.bounds[2])
    for (j in 1:length(expanded)){
      while (nchar(expanded[j]) < 4){
        expanded[j] <- paste0("0", expanded[j])
      }
    }
    expanded <- paste(expanded, collapse = " ")
    nodes[i] <- expanded
  }
  
  #expanded ranges were entered as space-delimited single elements in the list of nodes
  #fix nodes so that each element contains a single node address
  nodes <- unlist(strsplit(nodes, "\\s"))
  
  #concatenate to required form: xx-####
  nodes <- paste0(node.prefex, nodes)
  return(nodes)
}
