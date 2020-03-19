#' see if input is from function or dataset
#'
#' @param dataset_input 
#' the corresponding item from the dataset
#' @param function_input 
#' the corresponding item from the function
#' @return
#' 'both' if in both dataset and function 'function' if function only, 'dataset' if dataset only, 'none' if neither

check_input_condition <- function(dataset_input, function_input) {
  in_dataset <- length(dataset_input > 0)
  in_function <- !is.null(function_input)
  if (in_dataset & in_function) {
    return("both")
  }
  if (in_dataset & !in_function) {
    return("dataset")
  }
  if (in_function & !in_dataset) {
    return("function")
  }
  if (!in_function & !in_dataset) {
    return("none")
  }
}
