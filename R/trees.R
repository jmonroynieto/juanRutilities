
#' Get Parent Node 
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' returns the parent of the node.
#'
#' @param phylotable Phylo object or tibble from phylo  that is being queried
#' @param node The node you wish to find the parent for
#' @return An integer corresponding to the parent node
#' @export
get_parent <- function(phylotable, node) {
  if (any(class(phylotable) == "phylo")) {
    phylotable <- tibble::as_tibble(phylotable)
  }
  return(phylotable[phylotable$node==node,'parent'][[1]])
}

#' Get Children Nodes
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' returns an integer vector with the direct children of the node.
#'
#' @param phylotable Phylo object or tibble from phylo  that is being queried
#' @param parent The node you wish to find the children for
#' @return An integer vector corresponding to the children nodes
#' @export
get_children <- function(phylotable, parent) {
  if (any(class(phylotable) %in% c("phylo"))) {
    phylotable <- tibble::as_tibble(phylotable)
  }
  children_vector <- dplyr::pull(phylotable[phylotable$parent==parent,],node,name = label)
  if (length(children_vector)==0) {
    return(parent)
  }
  return(children_vector)
}

#' Get sibling Nodes
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' returns an integer vector with the direct siblings of the node.
#'
#' @param phylotable Phylo object or tibble from phylo  that is being queried
#' @param node The node you wish to find the siblings for
#' @return An integer vector corresponding to the sibling nodes; empty if the clade has no siblings
#' @export
get_siblings <- function(phylotable, node) {
  if (any(class(phylotable) %in% c("phylo"))) {
    phylotable <- tibble::as_tibble(phylotable)
  }
  sibling_vector <- unlist(sapply(get_children(phylotable, get_parent(phylotable, node)), function(x) {if (x!=node) {return(x)}}))
  if (length(sibling_vector)==0) {
    return(c())
  }
  return(sibling_vector)
}
#' Checks If Node Is Leaf 
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' a node to query. It returns a boolean indicating if the node is a leaf.
#' 
#'
#' @param phylotable Phylo object or tibble from phylo that is being queried
#' @param query The node you want to investigate
#' @return A boolean set to TRUE if the node has no further descendents
#' @export
is_leaf <- function(phylotable, query){
  if (any(class(phylotable) %in% c("phylo"))) {
    phylotable <- tibble::as_tibble(phylotable)
  }
  if (all(query == get_children(phylotable = phylotable, parent = query))) {
    return(TRUE)
  }
  return(FALSE)
}

#' Get the root of a phylogeny 
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' it returns the integer of the node with no further parents
#' 
#'
#' @param phylotable Phylo object or tibble from phylo that is being queried
#' @return An integer representing the root node of a phylogeny
#' @export
find_root <- function(phylotable) {
  rootX <- NULL
  node <- 1
  while (is.null(rootX)) {
    node2 <- get_parent(phylotable, node)
    if (node == node2) {
      rootX <- node
    }
    node <- node2
  }
  return(rootX)
}

#' Finds All DesendentsOf a Node
#'
#' This function takes in a phylo object or a tibble from a phylo object, an environment
#'  to populate, the expectation being that it is empty, and a node to query. It returns a
#'  NULL and modifies the environment that was passed. If no environment is passed, it assumes
#'  that there is an environment called answer_env. Because this function could not initialize
#'  and return this environment, there is a wrapper funtion that is exposed under the same name
#'  without the starting dot. The environment is set to contain a range of variables with the
#'  name of the node and a boolean assignment indicating if the node is a leaf. When testing,
#'  uncomment the showme() calls in the code bellow and use capture.output(type="messages").
#' 
#'
#' @param phylotable Phylo object or tibble from phylo that is being queried
#' @param ancestor The node you want to investigate
#' @param env The environment to which to write the answers
#' @return A boolean set to TRUE if the node has no further descendents
.get_offspring <-  function(phylotable, ancestor, env) {
  if (any(class(phylotable) %in% c("phylo"))) {
    phylotable <- tibble::as_tibble(phylotable)
  }
  #[tag-vEargp] This assignment happens for leaf nodes only and overwrites the previous FALSE value in the environment
  if (is_leaf(phylotable = phylotable, ancestor)) {
    child <- ancestor
    assign(envir = env, x= as.character(child), value = TRUE)
    #showme(sprintf("%s\t%s", child, TRUE))
    return()
  }
  for (child in get_children(phylotable = phylotable, ancestor)) {
    #[tag-vEargp] This assignment happens for all nodes
    if (ancestor!=child) {
      #showme(sprintf("%s\t%s", child, FALSE))
      assign(envir = env, x= as.character(child), value = FALSE)
      .get_offspring(phylotable = phylotable, ancestor = child, env = env)
      next
    }
  }
  return()
}
#' Get All Descendants of a Node
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' it returns the integer of the node with no further parents
#' 
#'
#' @param phylotable Phylo object or tibble from phylo that is being queried
#' @param ancestor The node crowning the clade you seek
#' @return An integer vector with all the nodes in the clade. This vector has the attribute "leaf" which describes with a boolean the status of each vector as a tip.
#' @export
get_offspring <- function(phylotable, ancestor) {
  answer_env <- new.env()
  .get_offspring(phylotable = phylotable, ancestor = ancestor, answer_env)
  answer_vector <- as.numeric(names(unlist(as.list(answer_env))))
  attr(answer_vector, which = "leaf") <- unname(unlist(as.list(answer_env)))
  rm(answer_env)
  return(answer_vector)
}
#' Prints nwk format recursively
#'
#'  depth is the key! -1 to print everything, 0 prints nothing
#' 
#'
#' @param phylotable Phylo object or tibble from phylo that is being queried
#' @param parent The node you want to investigate
#' @param depth OPTIONAL normally -1. This will return all. 0 returns nothing, one returns parent.
#' @return A character string in nwk format
.printClade_wLength <- function(phylotable,parent, depth=-1) {
  if (depth == 0) {return(c(""))}
  if (is_leaf(phylotable, parent)) { depth = 1 }
  children <- get_children(phylotable = phylotable, parent)
  names(children)[names(children)==""] <- as.character(children[names(children)==""])
  children_result <- unlist(lapply(children, function(x) {.printClade_wLength(phylotable,x, depth=ifelse(depth==-1,-1,depth-1))}))
  parent_label <- phylotable$label[phylotable$node == parent]
  parent_label <- ifelse(parent_label=="", parent, parent_label)
  lengthVal <- phylotable$branch.length[phylotable$node==parent]
  formatString <- ifelse(depth == 1, "%.0s%s:%f", "(%s)%s:%f")
  return(sprintf(formatString,paste0(children_result, collapse = ","), parent_label, lengthVal))
}
#' Prints nwk format recursively
#'
#'  depth is the key! -1 to print everything, 0 prints nothing
#' 
#'
#' @param phylotable Phylo object or tibble from phylo that is being queried
#' @param parent The node you want to investigate
#' @param depth OPTIONAL normally -1. This will return all. 0 returns nothing, one returns parent.
#' @return A string in newick format
.printClade <- function(phylotable,parent, depth=-1) {
  if (depth == 0) {return(c(""))}
  if (is_leaf(phylotable, parent)) { depth = 1 }
  if ("branch.length" %in% colnames(phylotable)) {
    return(.printClade_wLength(phylotable, parent, depth))
  }
  children <- get_children(phylotable = phylotable, parent)
  names(children)[names(children)==""] <- as.character(children[names(children)==""])
  children_result <- unlist(lapply(children, function(x) {.printClade(phylotable,x, depth=ifelse(depth==-1,-1,depth-1))}))
  parent_label <- phylotable$label[phylotable$node == parent]
  parent_label <- ifelse(parent_label=="", parent, parent_label)
  formatString <- ifelse(depth == 1, "%.0s%s", "(%s)%s")
  return(sprintf(formatString,paste0(children_result, collapse = ","), parent_label))
}



#' Transforms a table into a nwk string
#'
#' This function takes in a phylo object or a tibble from a phylo object and 
#' it returns its string representation in newick format.
#' The nodenumber is overwriten by the labels, if present.
#' 
#'
#' @param phylotable Phylo object or tibble from phylo. fmt as: parent,node,\<label\>
#' @return A string with a rooted nwk string
#' @export
transform_phylotable <- function(phylotable) {
  rootNode <- find_root(phylotable)
  clades <- get_children(phylotable = phylotable, parent = rootNode)
  results <- unlist(lapply(clades, function(x) {
    if (x==rootNode) {
      return()
    }
    return(.printClade(phylotable, x, depth = -1))
  }))
  sprintf("(%s)%s;",paste0(results, collapse = ","),rootNode)
}
#' Makes the necessary changes to a phylotable to customize nodes
#'
#' @param phylotable Phylo object or tibble from phylo. fmt as: parent,node,\<labe\>
#' @param old Int representing old node number
#' @param replacement Int value to place in 
#' @param checkForConflicts Bool value to limit the recursive calls. Do not use. this function operates under the expectation that once the nodes are replaced with random numbers, there should not be more conflicts.
#' @return a  phylo tibble with the requested changes
#' @export
replaceNode <- function(phylotable, old, replacement,.checkForConflicts = TRUE) {
  if (any(class(phylotable) %in% c("phylo"))) {
    phylotable <- tibble::as_tibble(phylotable)
  }
  #warning("You hardcoded values [tag-jXEPmG]")
  obstacles <- phylotable[phylotable$node==replacement,]
  present_already <- c(as.integer(phylotable$label[grepl(phylotable$label, perl = T, pattern ="\\d+")]),phylotable$node)
  unused_nodeIndexes <- Filter(x = sample(x = 1000:9000, replace = F, size = nrow(phylotable)*3),
                               f = function(x){!(x %in%  present_already)})
  if (.checkForConflicts == TRUE && nrow(obstacles) > 0) {
    for (i in 1:nrow(obstacles)) {
      phylotable <- replaceNode(phylotable,
                                old = obstacles[[i,'node']],
                                replacement = pop(unused_nodeIndexes),
                                .checkForConflicts = FALSE)
    }
  }
  #message(sprintf("changed %s with %s", phylotable[phylotable$node == old, 'node'], replacement))
  phylotable[phylotable$node == old, 'node'] <- replacement
  if (nrow(phylotable[phylotable$parent == old, 'parent'])>0) {
    phylotable[phylotable$parent == old, 'parent'] <- replacement
  }
  return(phylotable)
}

#' Returns The Last Common Ancestor of a Set of Tips
#' 
#' At least one of the 
#'
#' @param phylotable Phylo object or tibble from phylo. fmt as: parent,node,\<labe\>
#' @param tiplist Char vectorwith node lables to include.
#' @param nodelist Int vector with nodes to include.
#' @return Int representing last common ancestor
#' @export
#' 
#' 

find_LCA <-  function(phylotable, tiplist=NULL, nodelist=NULL) {
  #warning("you have hardcoded values.tag-LsHEpS"); phylotable <- set1; nodelist=NULL
  if (any(class(phylotable) %in% c("phylo"))) {
    phylo <- phylotable
    phylotable <- tibble::as_tibble(phylotable)
  } else {
    phylo <- ape::as.phylo(phylotable)
  }
  #this relies on NULL objects and integer(0) being ignored by append
  #warning("You have handcoded values.tag-SZwPUf"); tiplist <- eta
  tiplist <- phylotable$node[phylotable$label %in% tiplist]
  requestlist <- append(nodelist, tiplist)
  
  #Find root (perhaps change to do ntip+1)
  rootNode <- find_root(phylotable)
  #select n tips and find their moderates
  sampleSize <- ifelse(length(requestlist)<5,
                          length(requestlist),
                          ifelse(length(requestlist)>5 & length(requestlist)<21,
                                 length(requestlist) %/% 3,
                                 length(requestlist) %/% 10
                                 )
                      )
  sampleSize <- ifelse( 25 >= sampleSize, sampleSize, 20)
  nodesample  <- sample(x = requestlist,size = sampleSize, replace = F)
  #use the first tip to find path to root `rootWay`
  rootWay <- ape::nodepath(phy = as.phylo(phylotable), from = pop(nodesample), to = rootNode)
  #Find paths between the rest
  topStep <- 2
  while(length(nodesample) > 0){
    if (rootWay[topStep] == rootNode) {return(rootWay[topStep])}
    stepptingStones <- ape::nodepath(phy = phylo, from = pop(nodesample), to = rootWay[topStep])
    challenge <- utils::tail(which(rootWay %in% stepptingStones), n=1)
    if (challenge > topStep) { topStep <- challenge}
  }
  #find topmost node in 'rootWay' that is crossed
  #test if descendants grabs all the requested nodes If not find a path to the missing nodes and the root, return the topmost intersection.
  cadidate_offspring <- get_offspring(phylotable, ancestor = rootWay[topStep])
  missingNodes <- requestlist[!(requestlist %in% cadidate_offspring)]
  if (length(missingNodes) == 0) {
    return(rootWay[topStep])
  } else {
    return(find_LCA(phylotable,nodelist = c(rootWay[topStep], missingNodes)))
  }
}

#' Default coloring scale for all things CRYPTOCOCCUS
#'
#' @export
CryptoScale <- ggplot2::scale_color_manual(values = c(
  "VGI" = "#a87be0" ,
  "VGIIa" = "#654EF6" ,
  "VGIIb" = "#85A9FF" ,
  "VGIIc" = "#00CBF9" ,
  "VGII" = "#2C32A3" ,
  "VGIII" = "#167820" ,
  "VGIV" = "#e773ab" ,
  "VGV" = "#2c8c82" ,
  "VGVI" = "#F7703B" ,
  "VNI" = "#ec6f72" ,
  "VNII" = "#fded5e" ,
  "VNBI" = "#f6a979" ,
  "VNBII" = "#cb8d67" ,
  "VNIII" = "#bfb25f" ,
  "VNIV" = "#f57c4c"
), na.value = "black", name="Mol_Type", guide = ggplot2::guide_legend(override.aes = list(shape = "c")))
