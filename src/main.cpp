
#include "main.h"
#include "misc_v2.h"

using namespace std;

//------------------------------------------------
// estimate f by maximum likelihood
// [[Rcpp::export]]
Rcpp::List inbreeding_mle_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // extract inputs
  vector<vector<int>> x = rcpp_to_matrix_int(args["x"]);
  vector<double> f = rcpp_to_vector_double(args["f"]);
  vector<double> p = rcpp_to_vector_double(args["p"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // get dimensions
  int n = x.size();
  int L = x[0].size();
  
  // create vector q = 1-p
  vector<double> q(L);
  for (int j=0; j<L; ++j) {
    q[j] = 1.0 - p[j];
  }
  
  // create lookup tables
  int nf = int(f.size());
  vector<vector<double>> lookup_homo1(nf, vector<double>(L));
  vector<vector<double>> lookup_homo2(nf, vector<double>(L));
  vector<vector<double>> lookup_het(nf, vector<double>(L));
  for (int k=0; k<nf; ++k) {
    for (int j=0; j<L; ++j) {
      lookup_homo1[k][j] = log((1-f[k])*p[j]*p[j] + f[k]*p[j]);
      lookup_homo2[k][j] = log((1-f[k])*q[j]*q[j] + f[k]*q[j]);
      lookup_het[k][j] = log((1-f[k])*2*p[j]*q[j]);
    }
  }
  
  // create objects for storing results
  vector<double> loglike_vec(nf);
  vector<vector<vector<double>>> ret_all(n, vector<vector<double>>(n, vector<double>(nf)));
  vector<vector<double>> ret_ml(n, vector<double>(n));
  
  // loop through all pairwise samples
  for (int i1=0; i1<(n-1); ++i1) {
    
    // report progress
    if (report_progress) {
      update_progress(args_progress, "pb", i1, n-1);
    }
    
    for (int i2=(i1+1); i2<n; ++i2) {
      
      // calculate loglike for every value of f
      for (int k=0; k<nf; ++k) {
        double loglike = 0;
        for (int j=0; j<L; ++j) {
          if (x[i1][j] == -1 || x[i2][j] == -1) {
            continue;
          }
          if (x[i1][j] == 1 && x[i2][j] == 1) {
            loglike += lookup_homo1[k][j];
          } else if (x[i1][j] == 0 && x[i2][j] == 0) {
            loglike += lookup_homo2[k][j];
          } else {
            loglike += lookup_het[k][j];
          }
        }
        loglike_vec[k] = loglike;
      }
      
      // store all values
      ret_all[i1][i2] = loglike_vec;
      
      // store maximum likelihood f
      double best_loglike = loglike_vec[0];
      for (int k=1; k<nf; ++k) {
        if (loglike_vec[k] > best_loglike) {
          best_loglike = loglike_vec[k];
          ret_ml[i1][i2] = f[k];
        }
      }
      
    }  // end i2 loop
  }  // end i1 loop
  
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("ret_all") = ret_all, Rcpp::Named("ret_ml") = ret_ml);
  
}


//------------------------------------------------
// Perform gradient descient to calculate cluster Fi's
// [[Rcpp::export]]
Rcpp::List cluster_inbreeding_coef_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // extract inputs
  // to get Fis
  vector<string> indclst1 = rcpp_to_vector_string(args["clst1"]);
  vector<string> indclst2 = rcpp_to_vector_string(args["clst2"]);
  vector<string> locat_names = rcpp_to_vector_string(args["locat_names"]);
  vector<double> locat_fs = rcpp_to_vector_double(args["locat_fs"]);
  
  // data
  vector<double> gendist = rcpp_to_vector_double(args["gendist"]);
  vector<double> geodist = rcpp_to_vector_double(args["geodist"]);
  
  // int params
  double m = rcpp_to_double(args["m"]);
  double m_lowerbound = rcpp_to_double(args["m_lowerbound"]);
  double m_upperbound = rcpp_to_double(args["m_upperbound"]);
  int steps = rcpp_to_int(args["steps"]);
  double learningrate = rcpp_to_double(args["learningrate"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // dimensions
  int n = indclst1.size(); // n-1 samples
  int c = locat_names.size(); // n-1 clusters
  
  // items of interest to keep track of
  vector<double> cost(steps);
  fill(cost.begin(), cost.end(), 0);
  vector<double> m_run(steps);
  vector<vector<double>> fi_run(steps, vector<double>(c));
  
  // loop through steps
  for (int step = 0; step < steps; step++) {
    
    // report progress
    if (report_progress) {
      update_progress(args_progress, "pb", step, steps);
    }
    //-------------------------------
    // setup
    //-------------------------------
    // need to store sample level Fi's for indclst1
    vector<double> indclst1_fis(n);
    for (int ind = 0; ind < n; ind++) { // for each ind
      for (int clst = 0; clst < c; clst++){ // for each clust
        if (indclst1[ind] == locat_names[clst]) { // match ind  to clst param new fi
          indclst1_fis[ind] = locat_fs[clst];
        }
      }
    }
    // need to store sample level Fi's for indclst2
    vector<double> indclst2_fis(n);
    for (int ind = 0; ind < n; ind++) {
      for (int clst = 0; clst < c; clst++){
        if (indclst2[ind] == locat_names[clst]) { // match ind  to clst param new fi
          indclst2_fis[ind] = locat_fs[clst];
        }
      }
    }
    
    //-------------------------------
    // get current cost for given F and M
    //  NB this is being pulled from locat_fs and m
    //-------------------------------
    for (int ind = 0; ind < n; ind++) {
      cost[step] += pow( (gendist[ind] - ((indclst1_fis[ind] + indclst2_fis[ind])/2) *
        exp(-m*geodist[ind])), 2);
    }
    
    //-------------------------------
    // Run calculations to update F and M
    //-------------------------------
    //-------------------------------
    // F gradient
    //-------------------------------
    //  # N.B. expansion have to expand matrix in order for all sample i's to be accounted for in the gradient (fi + fj where j can be i)
    // clear results from previous step
    vector<double> fgrad(c);
    fill(fgrad.begin(), fgrad.end(), 0);
    double fpartial = 0;
    double floss = 0;
    
    // step through partial deriv and loss for gradient
    for (int clst = 0; clst < c; clst++){
      for (int ind = 0; ind < n; ind++) {
        if (indclst1[ind] == locat_names[clst] | indclst2[ind] == locat_names[clst]) {
          // gradient
          fpartial = -gendist[ind] * exp(-geodist[ind] * m) +
            ((indclst1_fis[ind] + indclst2_fis[ind])/2) *
            exp(-2*geodist[ind] * m);
          //loss
          floss = pow( (gendist[ind] - ((indclst1_fis[ind] + indclst2_fis[ind])/2) *
            exp(-m*geodist[ind])), 2);
          // gradient
          fgrad[clst] = fgrad[clst] + (fpartial * floss);
        }
      }
    }
    //-------------------------------
    // M gradient
    //-------------------------------
    //  # N.B. all terms included here, easier sum -- longer partial derivative
    // clear previous step results
    double mgrad = 0;
    
    // step through partial deriv and loss for gradient
    for (int ind = 0; ind < n; ind++) {
      // gradient
      double mpartial = 2 * gendist[ind] * geodist[ind] * ((indclst1_fis[ind] + indclst2_fis[ind])/2) *
        exp(-geodist[ind] * m) -
        2 * geodist[ind] *
        ((pow(indclst1_fis[ind], 2) + 2 * indclst1_fis[ind] * indclst2_fis[ind] + pow(indclst2_fis[ind], 2))/4) *
        exp(-2 * geodist[ind] * m);
      // loss
      double mloss = pow( (gendist[ind] - ((indclst1_fis[ind] + indclst2_fis[ind])/2) *
        exp(-m*geodist[ind])), 2);
      // gradient
      mgrad += mpartial * mloss;
    }
    
    //-------------------------------
    // Update F and M
    //-------------------------------
    // update F
    for (int clst = 0; clst < c; clst++){
      // update fs
      locat_fs[clst] = locat_fs[clst] - learningrate * fgrad[clst];
      // hard bounds on f
      if (locat_fs[clst] < 0) {
        locat_fs[clst] = 0;
      }
      if (locat_fs[clst] > 1) {
        locat_fs[clst] = 1;
      }
      // store for out
      fi_run[step][clst] = locat_fs[clst];
    }
    // update M
    m = m - learningrate * mgrad;
    // hard bounds for M
    if (m < m_lowerbound) {
      m = m_lowerbound;
    }
    if (m > m_upperbound) {
      m = m_upperbound;
    }
    // store for out
    m_run[step] = m;
    
  } // end steps
  
  //-------------------------------
  // Out
  //-------------------------------
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("m_run") = m_run,
                            Rcpp::Named("fi_run") = fi_run,
                            Rcpp::Named("cost") = cost,
                            Rcpp::Named("Final_Fis") = locat_fs,
                            Rcpp::Named("Final_m") = m
  );
  
}
























// sanity
