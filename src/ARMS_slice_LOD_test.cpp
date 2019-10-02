// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "arms/arms.c"
#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;
using namespace arma;

/* ************************ Create structures for normal log density ********************** */

struct norm_parm {
  double y, sigma_2_y;
  arma::vec mean_x_preds, beta, x_preds;
  arma::mat sigma_x_preds;
  int missing_pred;
};

double norm(double x, norm_parm d)
{
  arma::uvec missing_preds = find_nonfinite(d.x_preds);
  arma::vec obs_preds = d.x_preds;
  obs_preds.shed_row(as_scalar(missing_preds));
  arma::vec beta_obs = d.beta;
  beta_obs.shed_row(as_scalar(missing_preds));
  double mean_y = dot(obs_preds,beta_obs)+as_scalar((d.beta)(as_scalar(missing_preds)))*x;
  arma::vec x_preds = d.x_preds;
  x_preds(as_scalar(missing_preds)) = x;
  x_preds.shed_row(1); // remove intercept
  arma::vec x_minus_mu = x_preds - (d.mean_x_preds);
  return -0.5*(pow((d.y) - mean_y, 2.0)/(d.sigma_2_y))-0.5*as_scalar(x_minus_mu.t()*inv(d.sigma_x_preds)*x_minus_mu);
}

// Likelihood for 1+ LOD Biomarkers; will replace norm() for all numbers of LOD once function is finished and tested
double norm_multi(double x, norm_parm d)
{
  arma::vec obs_preds = d.x_preds;
  obs_preds.shed_row(d.missing_pred);
  arma::vec beta_obs = d.beta;
  beta_obs.shed_row(d.missing_pred);
  double mean_y = dot(obs_preds,beta_obs)+as_scalar((d.beta)(d.missing_pred))*x;
  arma::vec x_preds = d.x_preds;
  x_preds(d.missing_pred) = x;
  x_preds.shed_row(1); // remove intercept
  arma::vec x_minus_mu = x_preds - (d.mean_x_preds);
  return -0.5*(pow((d.y) - mean_y, 2.0)/(d.sigma_2_y))-0.5*as_scalar(x_minus_mu.t()*inv(d.sigma_x_preds)*x_minus_mu);
}

// Do same for ARMS sampling (need to refer to d as pointer)
double norm_arms(double x, void *norm_data)
{
  norm_parm *d;
  d = (norm_parm *)norm_data;
  arma::uvec missing_preds = find_nonfinite(d->x_preds);
  arma::vec obs_preds = d->x_preds;
  obs_preds.shed_row(as_scalar(missing_preds));
  arma::vec beta_obs = d->beta;
  beta_obs.shed_row(as_scalar(missing_preds));
  double mean_y = dot(obs_preds,beta_obs)+as_scalar((d->beta)(as_scalar(missing_preds)))*x;
  arma::vec x_preds = d->x_preds;
  x_preds(as_scalar(missing_preds)) = x;
  x_preds.shed_row(1); // remove intercept
  arma::vec x_minus_mu = x_preds - (d->mean_x_preds);
  return -0.5*(pow((d->y) - mean_y, 2.0)/(d->sigma_2_y))-0.5*as_scalar(x_minus_mu.t()*inv(d->sigma_x_preds)*x_minus_mu);
}

double norm_multi_arms(double x, void *norm_data)
{
  norm_parm *d;
  d = (norm_parm *)norm_data;
  arma::vec obs_preds = d->x_preds;
  obs_preds.shed_row(d->missing_pred);
  arma::vec beta_obs = d->beta;
  beta_obs.shed_row(d->missing_pred);
  double mean_y = dot(obs_preds,beta_obs)+as_scalar((d->beta)(d->missing_pred))*x;
  arma::vec x_preds = d->x_preds;
  x_preds(d->missing_pred) = x;
  x_preds.shed_row(1); // remove intercept
  arma::vec x_minus_mu = x_preds - (d->mean_x_preds);
  return -0.5*(pow((d->y) - mean_y, 2.0)/(d->sigma_2_y))-0.5*as_scalar(x_minus_mu.t()*inv(d->sigma_x_preds)*x_minus_mu);
}

/* ************************** Slice Sampling Code ***************/
double slice_fcs(	double (*myfunc)(double x, norm_parm d),
                  norm_parm d,
                  double init_x0,
                  double lower = -1000,
                  double upper = 1000,
                  int burn_in = 50) 
{
  
  double f0, x_1, f_x1, h0, L, R, U;
  double x_0 = init_x0;
  int m = 0;
  while(m<burn_in){
    m=m+1;
    // calculate current full conditional value;
    f0 = myfunc(x_0,d);
    
    // calculate height of the horizontal slice;
    h0 = f0 - ::Rf_rexp(1.0);		
    
    // Calculate initial horizontal interval;
    L = lower;
    R = upper;  
    
    // Step in;
    U = ::Rf_runif(0.0,1.0);
    x_1 = L+U*(R-L);
    f_x1 = myfunc(x_1,d);
    
    int stop=0;
    
    while(stop==0)
    {
      /* Rcpp::Rcout << "L=" << L << " - R=" << R << std::endl; */
      U = ::Rf_runif(0.0,1.0);
      x_1 = L+U*(R-L);
      f_x1 = myfunc(x_1,d);
      
      if(f_x1 >= h0){
        stop=1;
      }
      else{
        if(x_1<x_0){
          stop=0;
          L=x_1;
        }
        else{
          stop=0;
          R=x_1;
        }
      }
    }
    x_0 = x_1;
    // cout << x_1 << endl;
  }
  return x_1;
}

/* ************************** Slice Sampling Examples ***************/
// [[Rcpp::export]]
arma::vec slice_sample_normal_ex(arma::vec y_data, arma::vec x_data,
                              arma::vec mean_x_preds, arma::vec beta, double sigma_2_y, arma::mat sigma_x_preds,
                              double LOD_lower, double LOD_upper,
                              int sample_size){
  arma::vec sampled_data = x_data;
  norm_parm d;
  d.y = as_scalar(y_data);
  d.x_preds = x_data;
  d.mean_x_preds = mean_x_preds;
  d.beta = beta;
  d.sigma_2_y = sigma_2_y;
  d.sigma_x_preds = sigma_x_preds;
  
  double xprev;
    if(LOD_lower<0){
      xprev=LOD_upper*sqrt(2);
    }else{
      xprev=LOD_upper/sqrt(2);
    }
  
  double burn_in_initial = slice_fcs(norm, d, xprev, LOD_lower, LOD_upper);
  arma::vec sim_sample = vec(sample_size);
    for(int i=0; i<sample_size; i++){
      sim_sample(i)=slice_fcs(norm, d, burn_in_initial, LOD_lower, LOD_upper, 1);
    }
   return sim_sample;
  }

/* ************************** LOD Regression Fit ****************/
// [[Rcpp::export]]
List LOD_fit_test(arma::vec y_data, arma::mat x_data, 
                  arma::vec mean_x_preds, arma::vec beta, double sigma_2_y, 
                  arma::mat sigma_x_preds,
                  int no_of_samples,
                  double threshold, int max_iterations,
                  arma::vec LOD_u_l,
                  double ystart)
{
  /* Set starting point for slice sampler parameters */ 
  double xprev;
  if(LOD_u_l(1)<0){ 
    xprev=LOD_u_l(1)*sqrt(2);
  }else{
    xprev=LOD_u_l(1)/sqrt(2);
  }
  
  /* Set objects to index iterations and stop when finished */ 
  int iteration_index = 0;
  int betas_coverged = 1;
  int iterations_finished = 0;
  
  arma::uvec ids_LOD = find_nonfinite(x_data);
  arma::vec y_expand=zeros(no_of_samples*ids_LOD.n_elem+x_data.n_rows-ids_LOD.n_elem);
  arma::vec diff = zeros(beta.n_elem)+threshold+1.0;
  arma::mat x_data_return=mat(no_of_samples*ids_LOD.n_elem+x_data.n_rows-ids_LOD.n_elem,x_data.n_cols+1,fill::zeros);
  
  arma::vec weighted_mean_x_preds = mean_x_preds;
  arma::mat weighted_cov_x_preds = sigma_x_preds;
  
  arma::mat beta_estimates=mat(max_iterations, x_data.n_cols, fill::zeros);
  beta_estimates.row(0) = trans(beta);
  
  while(iteration_index < (max_iterations-1) && iterations_finished == 0){
    iteration_index = iteration_index + 1;
    unsigned row_index_count = 0;
    
    for(int i=0; i <x_data.n_rows ;i=i+1){
      norm_parm d;
      d.y = as_scalar(y_data(i));
      d.x_preds = trans(x_data.row(i));
      d.beta = trans(beta_estimates.row(iteration_index-1));
      d.sigma_2_y = sigma_2_y;
      
      d.sigma_x_preds = weighted_cov_x_preds;
      d.mean_x_preds = weighted_mean_x_preds;
      
      if( (d.x_preds).has_nan() == FALSE){
        row_index_count = row_index_count + 1;
        y_expand(row_index_count-1) = d.y;
        x_data_return(arma::uvec {row_index_count-1}, find_finite(d.x_preds)) = trans(d.x_preds);
        x_data_return(row_index_count-1, x_data_return.n_cols-1) = 1;
        
      } else {
        
        for(int s=0; s<no_of_samples; s=s+1){
          row_index_count = row_index_count + 1;
          
          y_expand(row_index_count-1) = d.y;
          
          /* Sample using slice sampler */ 
          double xsamp=slice_fcs(norm, d, xprev, LOD_u_l(0), LOD_u_l(1));

          x_data_return(arma::uvec {row_index_count-1},find_finite(d.x_preds))=trans(d.x_preds(find_finite(d.x_preds)));
          x_data_return(row_index_count-1,as_scalar(find_nonfinite(d.x_preds)))=xsamp;
          x_data_return(row_index_count-1, x_data_return.n_cols-1) = (1.0)/no_of_samples;
          
        }
      }
    }
    
    arma::mat x_data_analysis = x_data_return;
    arma::vec weight_vec = x_data_return.col(x_data_return.n_cols-1);
    x_data_analysis.shed_col(x_data_return.n_cols-1);
    arma::vec beta_hat_wls = inv(trans(x_data_analysis)*diagmat(weight_vec)*x_data_analysis)*trans(x_data_analysis)*diagmat(weight_vec)*y_expand;
    beta_estimates.row(iteration_index) = trans(beta_hat_wls);
    
    // Update mean and covariance matrix
    arma::mat weighted_data = x_data_analysis;
    weighted_data.shed_col(0);
    for(int k=0;k<weighted_data.n_rows; k=k+1){
      weighted_data.row(k) = as_scalar(weight_vec(k))*weighted_data.row(k)*(1.0)/x_data.n_rows;
    }
    
    weighted_mean_x_preds = trans(sum(weighted_data, 0));
    
    // Now covariance
    x_data_analysis.shed_col(0);
    arma::mat weighted_mean_x_mat = mat(x_data_analysis.n_rows,x_data_analysis.n_cols,fill::zeros);
    
    for(int k=0;k<weighted_mean_x_mat.n_rows; k=k+1){
      weighted_mean_x_mat.row(k) = trans(weighted_mean_x_preds);
    }
    
    weighted_cov_x_preds = ((1.0)/x_data.n_rows)*trans(diagmat(pow(weight_vec, 0.5))*(x_data_analysis-weighted_mean_x_mat))*(diagmat(pow(weight_vec, 0.5))*(x_data_analysis-weighted_mean_x_mat));
    
    for(int k=0; k<beta.n_elem; k=k+1){
      diff(k) = abs(beta_estimates(iteration_index,k)-beta_estimates(iteration_index-1,k));
      if(as_scalar(diff(k))>=threshold){
        betas_coverged=0;
      } else {
        betas_coverged=1;
      }
    }
    if(betas_coverged==1){
      iterations_finished=1;
    } else{
      iterations_finished=0;
    }
    
  }
  
  return List::create(Named("beta_estimates") = beta_estimates,
                      Named("beta_estimate_last_iteration") = beta_estimates.row(iteration_index));
}

// [[Rcpp::export]]
List bootstrap_test(int num_of_boots,
                    arma::vec y_data, arma::mat x_data, 
                    int no_of_samples, 
                    double threshold, int max_iterations,
                    arma::vec LOD_u_l,
                    double ystart){
  
  List boot_results(num_of_boots);

  // Now implement bootstrap
  for(int j=0; j<num_of_boots; j=j+1){
    // Resample data with replacement
    arma::uvec ids=conv_to<uvec>::from(linspace(0, y_data.n_elem-1, 
                                          y_data.n_elem));
    arma::uvec boot_ids = sample(ids, ids.n_elem, true);
    arma::vec y_data_boot = y_data.elem(boot_ids);
    arma::mat x_data_boot = x_data.rows(boot_ids);
    // Run analysis on bootstrap sample
    // First, need new set of inital values of parameters based on bootstrap sample
    // Do this using complete cases from bootstrap sample
    
    // Cycle through columns of X to find subjects with LOD values
    
    arma::uvec no_LOD_by_column = {x_data_boot.n_rows+1};
    
    for(int k=0; k<x_data_boot.n_cols; k=k+1){
      no_LOD_by_column.insert_rows(no_LOD_by_column.n_rows,
                                   find_nonfinite(x_data_boot.col(k)));
    }
    
    no_LOD_by_column.shed_row(0);
    
    // Create complete data objects
    arma::mat x_data_boot_complete = x_data_boot;
    arma::vec y_data_boot_complete = y_data_boot;
    
    if(no_LOD_by_column.n_elem>0){
      
      for(int k=0; k<no_LOD_by_column.n_elem; k++){
        x_data_boot_complete.shed_row(no_LOD_by_column(no_LOD_by_column.n_elem-k-1));
        y_data_boot_complete.shed_row(no_LOD_by_column(no_LOD_by_column.n_elem-k-1));
      }
    }
    
    
    // Calculate beta and sigma_2 (residual variance)
    arma::vec beta = inv(trans(x_data_boot_complete)*x_data_boot_complete)*trans(x_data_boot_complete)*y_data_boot_complete;
    double sum_squared_residuals = sum(pow(y_data_boot_complete-x_data_boot_complete*beta, 2.0));
    
    double sigma_2_y = (1.0/(x_data_boot_complete.n_rows-x_data_boot_complete.n_cols))*sum_squared_residuals;
    
    x_data_boot_complete.shed_col(0); //removing intercept
    
    arma::vec mean_x_preds = conv_to<vec>::from(mean(x_data_boot_complete, 0));
    arma::mat sigma_x_preds = cov(x_data_boot_complete, 0);
    boot_results[j] =
      LOD_fit_test(y_data_boot, x_data_boot,
                   mean_x_preds, beta, sigma_2_y,
                   sigma_x_preds,
                   no_of_samples, threshold, max_iterations,
                   LOD_u_l,
                   ystart)("beta_estimate_last_iteration");
  }
  
  return boot_results;
}

// LOD Multiple Covariate Testing
// [[Rcpp::export]]
List LOD_fit_multiple(arma::vec y_data, arma::mat x_data, 
                     arma::vec mean_x_preds, arma::vec beta, double sigma_2_y, 
                     arma::mat sigma_x_preds,
                     int no_of_samples,
                     double threshold, int max_iterations,
                     arma::mat LOD_u_l,
                     int sampler){
  
  /* Set starting point for slice sampler parameters */ 
  arma::vec xprev=vec(LOD_u_l.n_rows);
  for(int i=0; i<LOD_u_l.n_rows; i=i+1){
    if(LOD_u_l(i,1)<0){ 
      xprev(i)=as_scalar(LOD_u_l(i,1))*sqrt(2);
    }else{
      xprev(i)=as_scalar(LOD_u_l(i,1))/sqrt(2);
    }
  }
    
  /* Set objects to index iterations and stop when finished */
  int iteration_index = 0;
  int betas_coverged = 1;
  int iterations_finished = 0;

  arma::vec ids_LOD=zeros(x_data.n_rows);
   for(int i=0; i<x_data.n_rows; i=i+1){
     if((x_data.row(i)).has_nan()){
       ids_LOD(i)=1;
     }
   }
  
  arma::vec y_expand=zeros(no_of_samples*accu(ids_LOD)+x_data.n_rows-accu(ids_LOD));
  arma::vec diff = zeros(beta.n_elem)+threshold+1.0;
  arma::mat x_data_return=mat(no_of_samples*accu(ids_LOD)+x_data.n_rows-accu(ids_LOD),x_data.n_cols+1,fill::zeros);
  
  arma::vec weighted_mean_x_preds = mean_x_preds;
  arma::mat weighted_cov_x_preds = sigma_x_preds;
   
  arma::mat beta_estimates=mat(max_iterations, x_data.n_cols, fill::zeros);
  beta_estimates.row(0) = trans(beta);
   
  while(iteration_index < (max_iterations-1) && iterations_finished == 0){
     iteration_index = iteration_index + 1;
     unsigned row_index_count = 0;
     
     for(int i=0; i <x_data.n_rows ;i=i+1){
       norm_parm d;
       d.y = as_scalar(y_data(i));
       d.x_preds = trans(x_data.row(i));
       d.beta = trans(beta_estimates.row(iteration_index-1));
       d.sigma_2_y = sigma_2_y;
       
       d.sigma_x_preds = weighted_cov_x_preds;
       d.mean_x_preds = weighted_mean_x_preds;
       
       if( (d.x_preds).has_nan() == FALSE){
         row_index_count = row_index_count + 1;
         y_expand(row_index_count-1) = d.y;
         x_data_return(arma::uvec {row_index_count-1}, find_finite(d.x_preds)) = trans(d.x_preds);
         x_data_return(row_index_count-1, x_data_return.n_cols-1) = 1;
         
       } 
        else {
          
          // create vector to hold indicies of missing biomarkers
         arma::uvec missing_indicies = find_nonfinite(d.x_preds);

          for(int s=0; s<no_of_samples; s=s+1){
            row_index_count = row_index_count + 1;
            
            y_expand(row_index_count-1) = d.y;
             (d.x_preds)(missing_indicies) = xprev(missing_indicies); 

            /* Sample using slice sampler */ 
              double xsamp;
              if(sampler==0){
                  if(s==0){
                    for(int l=0; l<missing_indicies.n_elem; l=l+1){
                      d.missing_pred = missing_indicies(l); 
                      double burn_in_inital = slice_fcs(norm_multi, d, xprev(d.missing_pred), 
                                                 as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {0})), as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {1})));
                      xprev(d.missing_pred) = burn_in_inital;
                      (d.x_preds)(d.missing_pred)= burn_in_inital;
                    }
                    
                    for(int l=0; l<missing_indicies.n_elem; l=l+1){
                      d.missing_pred = missing_indicies(l); 
                      xsamp=slice_fcs(norm_multi, d, xprev(d.missing_pred), 
                                      as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {0})), as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {1})), 1);
                      (d.x_preds)(d.missing_pred)=xsamp;
                    }
                  }else{
                    for(int l=0; l<missing_indicies.n_elem; l=l+1){
                      d.missing_pred = missing_indicies(l);
                      xsamp=slice_fcs(norm_multi, d, xprev(d.missing_pred), 
                                             as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {0})), as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {1})), 1);
                      (d.x_preds)(d.missing_pred)=xsamp;
                    }
                  }
              }else{
                  /* Set ARMS parameters */ 
                for(int l=0; l<missing_indicies.n_elem; l=l+1){
                    
                  /* Sample using slice sampler */ 
                  d.missing_pred = missing_indicies(l);  
                  int ninit = 4;
                  double xl = as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {0})), xr=as_scalar(LOD_u_l(arma::uvec {d.missing_pred}, arma::uvec {1}));
                  int dometrop = 1;
                  double xprev;
                  if(xr<0){ 
                    xprev=xr*sqrt(2);
                  }else{
                    xprev=xr/sqrt(2);
                  }
                  double convex = 1.0;
                  int npoint = 100;
                  
                  int ncent = 0;
                  double qcent[20] = {.25, .5, .75}; //Random, not used
                  double xcent[20]; //Not actually used, see ncent
                  
                  int neval;
                  int nsamples = 1;
                  
                  double xinit_vec[4] = {0.2*(xr-xl) + xl,
                                         0.4*(xr-xl) + xl,
                                         0.6*(xr-xl) + xl,
                                         0.8*(xr-xl) + xl};
                  double err;
                  
                  err = arms(xinit_vec, ninit, &xl, &xr, &norm_multi_arms, &d, &convex, npoint, dometrop, 
                             &xprev, &xsamp, nsamples, qcent, xcent, ncent, &neval);
                  (d.x_preds)(d.missing_pred)=xsamp;
                }
              }
            
            x_data_return(arma::uvec {row_index_count-1}, find_finite(d.x_preds))=trans(d.x_preds);
            x_data_return(row_index_count-1, x_data_return.n_cols-1) = (1.0)/no_of_samples;
            
          }
        }
     }
     
     arma::mat x_data_analysis = x_data_return;
     arma::vec weight_vec = x_data_return.col(x_data_return.n_cols-1);
     x_data_analysis.shed_col(x_data_return.n_cols-1);
     arma::vec beta_hat_wls = inv(trans(x_data_analysis)*diagmat(weight_vec)*x_data_analysis)*trans(x_data_analysis)*diagmat(weight_vec)*y_expand;
     beta_estimates.row(iteration_index) = trans(beta_hat_wls);

    // Update mean and covariance matrix
    arma::mat weighted_data = x_data_analysis;
    weighted_data.shed_col(0);
    for(int k=0;k<weighted_data.n_rows; k=k+1){
      weighted_data.row(k) = as_scalar(weight_vec(k))*weighted_data.row(k)*(1.0)/x_data.n_rows;
    }

    weighted_mean_x_preds = trans(sum(weighted_data, 0));

    // Now covariance
    x_data_analysis.shed_col(0);
    arma::mat weighted_mean_x_mat = mat(x_data_analysis.n_rows,x_data_analysis.n_cols,fill::zeros);

    for(int k=0;k<weighted_mean_x_mat.n_rows; k=k+1){
      weighted_mean_x_mat.row(k) = trans(weighted_mean_x_preds);
    }

    weighted_cov_x_preds = ((1.0)/x_data.n_rows)*trans(diagmat(pow(weight_vec, 0.5))*(x_data_analysis-weighted_mean_x_mat))*(diagmat(pow(weight_vec, 0.5))*(x_data_analysis-weighted_mean_x_mat));

    for(int k=0; k<beta.n_elem; k=k+1){
      diff(k) = abs(beta_estimates(iteration_index,k)-beta_estimates(iteration_index-1,k));
      if(as_scalar(diff(k))>=threshold){
        betas_coverged=0;
      } else {
        betas_coverged=1;
      }
    }
    if(betas_coverged==1){
      iterations_finished=1;
    } else{
      iterations_finished=0;
    }

   }
  // 
  // return List::create(Named("beta_estimates") = beta_estimates,
  //                     Named("beta_estimate_last_iteration") = beta_estimates.row(iteration_index));
  
  return List::create(Named("y_expand_last_int") = y_expand,
                      Named("x_data_return_last_int") = x_data_return,
                      Named("beta_estimates") = beta_estimates,
                      Named("beta_estimate_last_iteration") = beta_estimates.row(iteration_index));
}

// [[Rcpp::export]]
List bootstrap_multi_test(int num_of_boots,
                    arma::vec y_data, arma::mat x_data, 
                    int no_of_samples, 
                    double threshold, int max_iterations,
                    arma::mat LOD_u_l,
                    int sampler){
  List boot_results(num_of_boots);
  
  // Now implement bootstrap
  for(int j=0; j<num_of_boots; j=j+1){
    // Resample data with replacement
    arma::uvec ids=conv_to<uvec>::from(linspace(0, y_data.n_elem-1, 
                                          y_data.n_elem));
    arma::uvec boot_ids = sample(ids, ids.n_elem, true);
    arma::vec y_data_boot = y_data.elem(boot_ids);
    arma::mat x_data_boot = x_data.rows(boot_ids);
    // Run analysis on bootstrap sample
    // First, need new set of inital values of parameters based on bootstrap sample
    // Do this using complete cases from bootstrap sample
    
    // Cycle through columns of X to find subjects with LOD values
    
    arma::uvec no_LOD_by_column = {x_data_boot.n_rows+1};
    
    for(int k=0; k<x_data_boot.n_cols; k=k+1){
      no_LOD_by_column.insert_rows(no_LOD_by_column.n_rows,
                                   find_nonfinite(x_data_boot.col(k)));
    }
    
    no_LOD_by_column.shed_row(0);
    arma::uvec no_LOD_by_column_unique = unique(no_LOD_by_column);
    
    // Create complete data objects
    arma::mat x_data_boot_complete = x_data_boot;
    arma::vec y_data_boot_complete = y_data_boot;

    if(no_LOD_by_column_unique.n_elem>0){

      for(int k=0; k<no_LOD_by_column_unique.n_elem; k++){
        x_data_boot_complete.shed_row(no_LOD_by_column_unique(no_LOD_by_column_unique.n_elem-k-1));
        y_data_boot_complete.shed_row(no_LOD_by_column_unique(no_LOD_by_column_unique.n_elem-k-1));
      }
    }


    // Calculate beta and sigma_2 (residual variance)
    arma::vec beta = inv(trans(x_data_boot_complete)*x_data_boot_complete)*trans(x_data_boot_complete)*y_data_boot_complete;
    double sum_squared_residuals = sum(pow(y_data_boot_complete-x_data_boot_complete*beta, 2.0));

    double sigma_2_y = (1.0/(x_data_boot_complete.n_rows-x_data_boot_complete.n_cols))*sum_squared_residuals;

    x_data_boot_complete.shed_col(0); //removing intercept

    arma::vec mean_x_preds = conv_to<vec>::from(mean(x_data_boot_complete, 0));
    arma::mat sigma_x_preds = cov(x_data_boot_complete, 0);
    boot_results[j] =
      LOD_fit_multiple(y_data_boot, x_data_boot,
                   mean_x_preds, beta, sigma_2_y,
                   sigma_x_preds,
                   no_of_samples, threshold, max_iterations,
                   LOD_u_l,
                   sampler)("beta_estimate_last_iteration");
    // boot_results[j] = List::create(no_LOD_by_column,
    //                                no_LOD_by_column_unique);
  }
  
  return boot_results;
}