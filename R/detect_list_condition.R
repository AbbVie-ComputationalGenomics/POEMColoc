#' Title
#'
#' @param dataset 
#' POEMColoc input dataset
#' @return
#' TRUE if a list of datasets, FALSE if a single dataset
#'
#' @examples
detect_list_condition <- function(dataset) {
  if (!is.list(dataset)) {
    return(FALSE)
  } else {
    if (length(dataset)==0) {
      return(TRUE)
    } else {
      if (is.list(dataset[[1]])) {
        if (!all(c("pos", "chr") %in% names(dataset[[1]]))) {
          stop("Input is ill-formated")
        } else {
          return(TRUE)
        }
      } else {
        return(FALSE)
      }
    }
  }
}
