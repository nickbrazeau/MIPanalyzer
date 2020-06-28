
#' @title Identify Inbreeding Coefficients of Locations (Demes) in Continuous Space
#' @param K_gendist_geodist dataframe; The genetic-geographic data by deme (K)
#' @param start_params named numeric vector; vector of start parameters. 
#' @param learningrate numeric; alpha parameter for how much each "step" is weighted in the gradient descent
#' @param m_lowerbound numeric; a lower bound for the m parameter
#' @param m_upperbound numeric; an upper bound for the m parameter
#' @param steps numeric; the number of "steps" as we move down the gradient
#' @param report_progress boolean; whether or not a progress bar should be shown as you iterate through steps
#' @details The gen.geo.dist dataframe must be named with the following columns: 
#'          "smpl1"; "smpl2"; "locat1"; "locat2"; "gendist"; "geodist"; which corresponds to: 
#'          Sample 1 Name; Sample 2 Name; Sample 1 Location; Sample 2 Location; 
#'          Pairwise Genetic Distance; Pairwise Geographpic Distance. Note, the order of the 
#'          columns do not matter but the names of the columns must match.
#' @details The start_params vector names must match the cluster names (i.e. clusters must be 
#'          have a name that we can match on for the starting relatedness paramerts). In addition,
#'          you must provide a start parameter for "m". 
#' @description The purpose of this statistic is to identify an inbreeding coefficient, or degree of 
#'              relatedness, for a given location in space. We assume that locations in spaces can be 
#'              represented as "clusters," such that multiple individuals live in the same cluster 
#'              (i.e. samples are sourced from the same location). The expected pairwise relationship
#'              between two individuals, or samples, is dependent on the each sample's cluster's inbreeding
#'              coefficient and the geographic distance between the clusters.  
#' @export
#'          

deme_inbreeding_coef <- function(K_gendist_geodist, 
                                 start_params = c(), 
                                 m_lowerbound = 0,
                                 m_upperbound = 1,
                                 f_learningrate = 1e-4,
                                 m_learningrate = 1e-10,
                                 steps = 1e3,
                                 report_progress = TRUE){
  
  #..............................................................
  # Assertions & Catches
  #..............................................................
  if (!all(colnames(K_gendist_geodist) %in% c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))) {
    stop("The K_gendist_geodist must contain columns with names: smpl1, smpl2, locat1, locat2, gendist, geodist")
  }
  locats <- names(start_params)[!grepl("m", names(start_params))]
  if (!all(unique(c(K_gendist_geodist$locat1, K_gendist_geodist$locat2)) %in% locats)) {
    stop("You have cluster names in your K_gendist_geodist dataframe that are not included in your start parameters")
  }
  if (!any(grepl("m", names(start_params)))) {
    stop("A m start parameter must be provided (i.e. there must be a vector element name m in your start parameter vector)")
  }
  assert_dataframe(K_gendist_geodist)
  assert_vector(start_params)
  assert_single_numeric(f_learningrate)
  assert_single_numeric(m_learningrate)
  assert_single_int(steps)
  assert_single_logical(report_progress)
  
  #..............................................................
  # setup and create progress bars
  #..............................................................
  # note, samples names there just for sanity
  pb <- txtProgressBar(min = 0, max = steps, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  
  #............................................................
  # R manipulations before C++
  #...........................................................
  expand_pairwise <- function(y){
    yexpand <- y
    colnames(yexpand) <- c("smpl2", "smpl1", "locat2", "locat1", "gendist", "geodist")
    yexpand <- rbind.data.frame(y, yexpand) # now have all pairwise possibilities
    yexpand <- yexpand[!duplicated(yexpand), ] # remove duplicate selfs
    return(yexpand)
  }
  # use efficient R functions to group pairs and wrangle data for faster C++ manipulation
  # get deme names and lift over sorted names for i and j
  demes <- sort(unique(c(K_gendist_geodist$locat1, K_gendist_geodist$locat2)))
  keyi <- data.frame(locat1 = demes, i = 1:length(demes))
  keyj <- data.frame(locat2 = demes, j = 1:length(demes))
  
  # get genetic data by pairs through efficient nest
  gendist <- K_gendist_geodist %>% 
    expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("locat1", "locat2", "gendist")) %>% 
    dplyr::group_by_at(c("locat1", "locat2")) %>% 
    tidyr::nest(.) %>% 
    dplyr::arrange_at(c("locat1", "locat2")) %>% 
    dplyr::left_join(., keyi, by = "locat1") %>% 
    dplyr::left_join(., keyj, by = "locat2")
  
  
  # put gendist into an array
  n_Kpairmax <- max(purrr::map_dbl(gendist$data, nrow))
  gendist_arr <- array(data = -1, dim = c(length(locats), length(locats), n_Kpairmax))
  for (i in 1:nrow(gendist)) {
    gendist_arr[gendist$i[i], gendist$j[i], 1:nrow(gendist$data[[i]])] <- unname(unlist(gendist$data[[i]]))
  }
  
  # put geo information into distance matrix
  geodist <- K_gendist_geodist %>% 
    expand_pairwise(.) %>% # get all pairwise for full matrix
    dplyr::select(c("locat1", "locat2", "geodist")) %>% 
    dplyr::group_by_at(c("locat1", "locat2")) %>% 
    tidyr::nest(.) %>% 
    dplyr::arrange_at(c("locat1", "locat2")) 
  geodist$data <- purrr::map_dbl(geodist$data, function(x){unique(x[[1]])})
  geodist <- geodist %>% 
    dplyr::left_join(., keyi, by = "locat1") %>% 
    dplyr::left_join(., keyj, by = "locat2")
  
  # upper tri
  geodist_mat <- matrix(data = -1, nrow = length(locats), ncol = length(locats))
  for (i in 1:nrow(geodist)) {
    geodist_mat[geodist$i[i], geodist$j[i]] <- geodist$data[i]
  }
  diag(geodist_mat) <- 0
  
  #..............................................................
  # run efficient C++ function
  #..............................................................
  
  args <- list(gendist = as.vector(gendist_arr),
               geodist = as.vector(geodist_mat),
               fvec = unname( start_params[!grepl("m", names(start_params))] ),
               n_Kpairmax = n_Kpairmax,
               m = unname(start_params["m"]),
               m_lowerbound = m_lowerbound,
               m_upperbound = m_upperbound,
               f_learningrate = f_learningrate,
               m_learningrate = m_learningrate,
               steps = steps,
               report_progress = report_progress
  )
  
  args_functions <- list(update_progress = update_progress)
  output_raw <- cluster_inbreeding_coef_cpp(args, args_functions, args_progress)
  
  # process output
  colnames(keyi) <- c("Deme", "key")
  if (Sys.getenv("DEME_INBREED_DEBUG") == "TRUE") {
    output <- list(
      deme_key = keyi,
      m_run = output_raw$m_run,
      fi_run = output_raw$fi_run,
      cost = output_raw$cost,
      Final_Fis = output_raw$Final_Fis,
      Final_m = output_raw$Final_m)
  } else {
    output <- list(
      deme_key = keyi,
      cost = output_raw$cost,
      Final_Fis = output_raw$Final_Fis,
      Final_m = output_raw$Final_m)
  }
  
  # return list
  return(output)
  
}
