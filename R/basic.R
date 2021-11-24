#' @export
push <- function(x, values) (assign(as.character(substitute(x)), c(x, values), parent.frame()))
#' @export
pop <- function(x) {last  = x[-1:-length(x)+1]; assign(as.character(substitute(x)), x[-length(x)], parent.frame()); return(last)}
#' @export
showme <- function(x) { message(paste0(x, collapse = '\n'))}
#' @export
random_char <- function(length=6) { result <- rawToChar(as.raw(sample(c(65:90,97:122), length, replace=T))); return(result)}
#' @export
add_path<- function(){showme('.libPaths( c("/Users/jmonroy-nieto/test/R/4/trees", .libPaths()))')}
#' @export
df.fromNAMES <- function(rows, columns) {
  x <- data.frame(matrix(NA, nrow = length(rows), ncol=length(columns)))
  row.names(x) <- rows
  colnames(x) <- columns
  return(x)
}
