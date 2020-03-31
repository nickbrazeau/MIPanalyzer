
#' @title Identify Inbreeding Coefficients of Locations in Continuous Space
#' @param gen.geo.dist dataframe; The genetic-geographic ata
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

cluster_inbreeding_coef <- function(clst_gendist_geodist, 
                                    start_params = c(), 
                                    m_lowerbound = 0,
                                    m_upperbound = 1,
                                    learningrate = 0.05, steps = 100,
                                    report_progress = TRUE){
  
  #..............................................................
  # Assertions & Catches
  #..............................................................
  if (!all(colnames(clst_gendist_geodist) %in% c("smpl1", "smpl2", "locat1", "locat2", "gendist", "geodist"))) {
    stop("The clst_gendist_geodist must contain columns with names: smpl1, smpl2, locat1, locat2, gendist, geodist")
  }
  locats <- names(start_params)[!grepl("m", names(start_params))]
  if (!all(unique(c(clst_gendist_geodist$locat1, clst_gendist_geodist$locat2)) %in% locats)) {
    stop("You have cluster names in your clst_gendist_geodist dataframe that are not included in your start parameters")
  }
  if (!any(grepl("m", names(start_params)))) {
    stop("A m start parameter must be provided (i.e. there must be a vector element name m in your start parameter vector)")
  }
  assert_dataframe(clst_gendist_geodist)
  assert_vector(start_params)
  assert_single_numeric(learningrate)
  assert_single_int(steps)
  assert_single_logical(report_progress)
  
  #..............................................................
  # setup and create progress bars
  #..............................................................
  # note, samples names there just for sanity
  pb <- txtProgressBar(min = 0, max = steps, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  #..............................................................
  # run efficient C++ function
  #..............................................................
  args <- list(clst1 = as.character(clst_gendist_geodist[, "locat1"]),
               clst2 = as.character(clst_gendist_geodist[, "locat2"]),
               gendist = clst_gendist_geodist[, "gendist"],
               geodist = clst_gendist_geodist[, "geodist"],
               locat_names = names(start_params)[!grepl("m", names(start_params))],
               locat_fs = unname( start_params[!grepl("m", names(start_params))] ),
               m = unname(start_params[length(start_params)]),
               m_lowerbound = m_lowerbound,
               m_upperbound = m_upperbound,
               learningrate = learningrate,
               steps = steps,
               report_progress = report_progress
  )
  
  args_functions <- list(update_progress = update_progress)
  output_raw <- cluster_inbreeding_coef_cpp(args, args_functions, args_progress)
  
  # process output
  output <- output_raw
  
  # return list
  return(output)
  
}
