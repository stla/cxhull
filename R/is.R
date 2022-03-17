isBoolean <- function(x){
  is.atomic(x) && is.logical(x) && length(x) == 1L && !is.na(x)
}