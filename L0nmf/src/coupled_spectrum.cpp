#include <Rcpp.h>
#include <L0glm/scalar_families.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
// Pooled shared-spectrum score and Fisher information on the original scale.
// [[Rcpp::export]]
List slcn_pooled_terms(NumericMatrix y,NumericMatrix mean,NumericVector design,
 int family,NumericVector phi,NumericVector theta,NumericMatrix weights,int nthreads=1,double censor=.5) {
 const int n=y.nrow(),m=y.ncol();
 if(mean.nrow()!=n||mean.ncol()!=m||design.size()!=n||weights.nrow()!=n||weights.ncol()!=m)
   stop("Incompatible pooled-spectrum dimensions.");
 if(!(phi.size()==1||phi.size()==m)||!(theta.size()==1||theta.size()==m)||
    nthreads<1||family<0||family>6)stop("Invalid pooled-spectrum parameters.");
 NumericVector gradient(m),information(m);
 #ifdef _OPENMP
 #pragma omp parallel for schedule(static) num_threads(nthreads)
 #endif
 for(int j=0;j<m;++j) {
   long double g=0,I=0;
   for(int i=0;i<n;++i)if(weights(i,j)>0&&design[i]>0) {
     auto z=L0glmScalar::derivatives(y(i,j),std::max(1e-15,mean(i,j)),family,phi[phi.size()==1?0:j],theta[theta.size()==1?0:j],censor);
     g+=weights(i,j)*design[i]*z.gradient;
     I+=weights(i,j)*design[i]*design[i]*z.information;
   }
   gradient[j]=(double)g;information[j]=(double)I;
 }
 return List::create(_["gradient"]=gradient,_["information"]=information);
}
