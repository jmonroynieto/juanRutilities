
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
  children_vector <- pull(phylotable[phylotable$parent==parent,'node'],node)
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
