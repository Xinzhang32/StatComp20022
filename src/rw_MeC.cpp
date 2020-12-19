#include <Rcpp.h>
using namespace Rcpp;

//' @title A markov sampler using Rcpp
//' @description a random walk Metropolis sampler for generating the standard Laplace distribution
//' @param sigma the variance
//' @param x0 the Initial value
//' @param N Number of random numbers
//' @return Returns random sequence and rejection probability
//' @examples
//' \dontrun{
//' rw_MeC(5,7,500)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rw_MeC(double sigma, double x0, int N)
{
  //Metropolis Randomwalk using C
  NumericVector x(N);
  x[0] = x0;
  double u, y;
  int k = 0;
  for (int i = 1; i < N; i++)
  {
    y = rnorm(1, x[i-1], sigma)[0];
    u = runif(1)[0];
    if (u <= exp(-(abs(y) - abs(x[i-1]))))
    {
      x[i] = y;
    }
    else
    {
      x[i] = x[i-1];
      k++;
    }
  }
  return x;
}
