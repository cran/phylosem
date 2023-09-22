
#define TMB_LIB_INIT R_init_phylosem
#include <TMB.hpp>

// SparseMatrix for Ornstein-Uhlenbeck network correlations
//template<class Type>
//Eigen::SparseMatrix<Type> Q_network( Type log_alpha,
//                                     int n_s,
//                                     vector<int> parent_s,
//                                     vector<int> child_s,
//                                     vector<Type> dist_s ){
//
//  Eigen::SparseMatrix<Type> Q( n_s, n_s );
//  Type alpha = exp( log_alpha );
//  for(int s=0; s<n_s; s++){
//    Q.coeffRef( s, s ) = Type(1.0);
//  }
//  for(int s=1; s<parent_s.size(); s++){
//    if( exp(-dist_s(s))!=0 ){
//      Q.coeffRef( parent_s(s), child_s(s) ) = -exp(-alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//      Q.coeffRef( child_s(s), parent_s(s) ) = Q.coeffRef( parent_s(s), child_s(s) );
//      Q.coeffRef( parent_s(s), parent_s(s) ) += exp(-2*alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//      Q.coeffRef( child_s(s), child_s(s) ) += exp(-2*alpha*dist_s(s)) / (1-exp(-2*alpha*dist_s(s)));
//    }
//  }
//  return Q;
//}

// ICAR for BM with sum-to-zero constraint
//  See: C:\Users\James.Thorson\Desktop\Work files\AFSC\2022-09 -- ICAR specification

// Precision of evolutionary covariance
//template<class Type>
//Eigen::SparseMatrix<Type> Q_sem( vector<Type> beta_z,
//                                 matrix<int> RAM,
//                                 int n_vars ){
//
//  // Define temporary objects
//  Eigen::SparseMatrix<Type> Q_vv( n_vars, n_vars );
//  // SEM
//  Eigen::SparseMatrix<Type> Linv_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Rho_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Gamma_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> Gammainv_vv(n_vars, n_vars);
//  Eigen::SparseMatrix<Type> I_vv( n_vars, n_vars );
//  Rho_vv.setZero();
//  Gamma_vv.setZero();
//  I_vv.setIdentity();
//  for(int zI=0; zI<RAM.rows(); zI++){
//    if(RAM(zI,0)==1) Rho_vv.coeffRef( RAM(zI,1)-1, RAM(zI,2)-1 ) = beta_z(RAM(zI,3)-1);
//    if(RAM(zI,0)==2) Gamma_vv.coeffRef( RAM(zI,1)-1, RAM(zI,2)-1 ) = beta_z(RAM(zI,3)-1); // Cholesky of covariance, so -Inf to Inf;
//  }
//  Gammainv_vv = atomic::matinv( Gamma_vv );
//  Linv_vv = Gammainv_vv * ( I_vv - Rho_vv );
//  Q_vv = Linv_vv.transpose() * Linv_vv;
//  return Q_vv;
//}

// Evolutionary covariance
//template<class Type>
//matrix<Type> V_sem( vector<Type> beta_z,
//                                 matrix<int> RAM,
//                                 vector<Type> RAMstart,
//                                 int n_vars ){
//
//  // Define temporary objects
//  matrix<Type> V_vv(n_vars, n_vars);
//  // SEM
//  matrix<Type> L_vv(n_vars, n_vars);
//  matrix<Type> Rho_vv(n_vars, n_vars);
//  matrix<Type> Gamma_vv(n_vars, n_vars);
//  matrix<Type> I_vv( n_vars, n_vars );
//  Rho_vv.setZero();
//  Gamma_vv.setZero();
//  I_vv.setIdentity();
//  Type tmp;
//  for(int r=0; r<RAM.rows(); r++){
//    // Extract estimated or fixed value
//    if(RAM(r,3)>=1){
//      tmp = beta_z(RAM(r,3)-1);
//    }else{
//      tmp = RAMstart(r);
//    }
//    // Assign to proper matrix
//    if(RAM(r,0)==1){
//      Rho_vv( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }else{
//      Gamma_vv( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
//    }
//  }
//  L_vv = I_vv - Rho_vv;
//  L_vv = atomic::matinv( L_vv );
//  L_vv = L_vv * Gamma_vv;
//  V_vv = L_vv * L_vv.transpose();
//  return V_vv;
//}

template<class Type>
Type objective_function<Type>::operator() ()
{
  //using namespace Eigen;
  using namespace density;

  // Data
  DATA_INTEGER( n_tip );
  DATA_IMATRIX( edge_ez );
  DATA_VECTOR( length_e );
  DATA_IMATRIX( RAM );
  DATA_VECTOR( RAMstart );
  DATA_INTEGER( estimate_ou );
  DATA_INTEGER( estimate_lambda );
  //DATA_INTEGER( estimate_kappa );
  DATA_VECTOR( height_v );
  DATA_MATRIX( y_ij );
  DATA_UPDATE( y_ij ); // Experiment with cAIC
  DATA_IVECTOR( v_i );
  DATA_IVECTOR( familycode_j );

  // Parameters
  PARAMETER_VECTOR( beta_z );
  PARAMETER_VECTOR( lnsigma_j );
  PARAMETER( lnalpha );
  PARAMETER( logitlambda );
  PARAMETER( lnkappa );
  PARAMETER_MATRIX( x_vj );
  PARAMETER_VECTOR( xbar_j );

  // Indices
  int n_e = edge_ez.rows();   // edges
  int n_v = x_vj.rows();      // vertices = edges + 1
  int n_j = x_vj.cols();      // variables
  int n_i = y_ij.rows();      // data
  int vparent, vchild;
  int vroot = n_tip;          // vroot = n_tip+1 - 1, where latter -1 is converting from R to CPP indexing

  // JNLL
  Type jnll = 0;
  vector<Type> jnll_v( n_v );
  matrix<Type> jnll_ij( n_i, n_j );
  jnll_ij.setZero();
  jnll_v.setZero();

  // Global vars
  Type alpha = exp( lnalpha );
  Type lambda = invlogit( logitlambda );
  Type kappa = exp( lnkappa );
  matrix<Type> yhat_ij( n_i, n_j );
  vector<Type> rho_v( n_v );
  vector<Type> var_v( n_v );
  vector<Type> sigma_j( n_j );
  sigma_j = exp( lnsigma_j );
  matrix<Type> Vtmp( n_j, n_j );
  vector<Type> xtmp_j( n_j );

  // Evolutionary inverse-covariance
  matrix<Type> V_jj( n_j, n_j );
  // vector<Type> beta_z, matrix<int> RAM, vector<Type> RAMstart, int n_vars
  //V_jj = V_sem( beta_z, RAM, RAMstart, n_j );

  // Assemble covariance
  // SEM
  matrix<Type> L_jj(n_j, n_j);
  matrix<Type> Rho_jj(n_j, n_j);
  matrix<Type> Gamma_jj(n_j, n_j);
  matrix<Type> I_jj( n_j, n_j );
  Rho_jj.setZero();
  Gamma_jj.setZero();
  I_jj.setIdentity();
  Type tmp;
  for(int r=0; r<RAM.rows(); r++){
    // Extract estimated or fixed value
    if(RAM(r,3)>=1){
      tmp = beta_z(RAM(r,3)-1);
    }else{
      tmp = RAMstart(r);
    }
    // Assign to proper matrix
    if(RAM(r,0)==1){
      Rho_jj( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    }else{
      Gamma_jj( RAM(r,1)-1, RAM(r,2)-1 ) = tmp;
    }
  }
  L_jj = I_jj - Rho_jj;
  L_jj = atomic::matinv( L_jj );
  L_jj = L_jj * Gamma_jj;
  V_jj = L_jj * L_jj.transpose();

  // Distribution of OU evolution -- Root
  // Correlation between i and parent(i) as distance -> INF
  if( estimate_ou==1 ){
    rho_v(vroot) = 0;
    // SD of Ornstein-Uhlenbeck process as distance -> INF
    var_v(vroot) = Type(1.0) / (2*alpha);
    // conditional probability
    Vtmp = V_jj * var_v(vroot);
    //xtmp_j = x_vj.row(root).array() - xbar_j;
    for(int j=0; j<n_j; j++) xtmp_j(j) = x_vj(vroot,j) - xbar_j(j);
    jnll_v(vroot) = MVNORM(Vtmp)( xtmp_j );
    // Optionally fix the root at the mean
  }else{
    rho_v(vroot) = NAN;
    var_v(vroot) = NAN;
  }
  // Distribution of OU evolution -- Edges
  for(int e=0; e<n_e; e++){
    vchild = edge_ez(e,1);
    vparent = edge_ez(e,0);
    if( estimate_ou==1 ){
      // Correlation between i and parent(i)
      rho_v(vchild) = exp( -alpha * pow(length_e(e),kappa) );
      // SD of O-U process
      var_v(vchild) = Type(1.0)/(2*alpha) * (Type(1.0)-exp( -2 * alpha * pow(length_e(e),kappa) ));
    }else{
      rho_v(vchild) = Type(1.0);
      var_v(vchild) = pow(length_e(e),kappa);
    }
    // conditional probability
    if( estimate_lambda==1 ){
      if( vchild < n_tip ){
        Vtmp = V_jj * ( lambda*var_v(vchild) + (1-lambda)*height_v(vchild) );
      }else{
        Vtmp = V_jj * ( lambda*var_v(vchild) );
      }
    }else{
      Vtmp = V_jj * var_v(vchild);
    }
    //xtmp_j = (x_vj.row(vchild).array()-xbar_j) - rho_v(vchild)*(x_vj.row(vparent).array()-xbar_j);
    for(int j=0; j<n_j; j++) xtmp_j(j) = (x_vj(vchild,j)-xbar_j(j)) - rho_v(vchild)*(x_vj(vparent,j)-xbar_j(j));
    jnll_v(vchild) = MVNORM(Vtmp)( xtmp_j );
  }
  jnll += jnll_v.sum();

  // Distribution for data
  for(int i=0; i<n_i; i++){
  for(int j=0; j<n_j; j++){
    if( !R_IsNA(asDouble(y_ij(i,j))) ){
      // familycode = 0 :  don't include likelihood
      // familycode = 1 :  normal
      if( familycode_j(j)==1 ){
        yhat_ij(i,j) = x_vj(v_i(i),j);
        jnll_ij(i,j) -= dnorm( y_ij(i,j), yhat_ij(i,j), sigma_j(j), true );
      }
      // familycode = 2 :  binomial
      if( familycode_j(j)==2 ){
        yhat_ij(i,j) = invlogit(x_vj(v_i(i),j));
        jnll_ij(i,j) -= dbinom( y_ij(i,j), Type(1.0), yhat_ij(i,j), true );
      }
      // familycode = 3 :  Poisson
      if( familycode_j(j)==3 ){
        yhat_ij(i,j) = exp(x_vj(v_i(i),j));
        jnll_ij(i,j) -= dpois( y_ij(i,j), yhat_ij(i,j), true );
      }
      // familycode = 4 :  Gamma:   shape = 1/CV^2; scale = mean*CV^2
      if( familycode_j(j)==4 ){
        yhat_ij(i,j) = exp(x_vj(v_i(i),j));
        jnll_ij(i,j) -= dgamma( y_ij(i,j), pow(sigma_j(j),-2), yhat_ij(i,j)*pow(sigma_j(j),2), true );
      }
    }
  }}
  jnll += jnll_ij.sum();

  // Calculate intercept
  vector<Type> root_j = x_vj.row(vroot);
  vector<Type> intercept_j( n_j );
  intercept_j = (I_jj - Rho_jj) * root_j.matrix();
  REPORT( root_j );
  REPORT( intercept_j );
  ADREPORT( intercept_j );

  // Reporting
  REPORT( rho_v );
  REPORT( var_v );
  REPORT( V_jj );
  REPORT( Rho_jj );
  REPORT( Gamma_jj );
  REPORT( jnll );
  REPORT( jnll_v );
  REPORT( jnll_ij );
  REPORT( alpha );
  REPORT( x_vj );
  REPORT( yhat_ij );  // Testing for cAIC
  ADREPORT( Rho_jj );
  ADREPORT( alpha );
  ADREPORT( lambda );
  ADREPORT( kappa );
  return jnll;
}
