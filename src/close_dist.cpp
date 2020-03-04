#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector dist_to_parents(std::string comp, StringVector parents) {
  IntegerVector dist(parents.size());
  for (size_t i = 0; i < parents.size(); i++) {
    for (size_t j = 0; j < comp.size(); j++) {
      dist[i] += (comp[j] != parents[i][j]);
    }
  }
  return dist;
}

// [[Rcpp::export]]
List calcDist(CharacterVector dna, int i_start, int i_end, bool progress=false) {
  // pairwise iterate over x
  std::vector<int> st1, st2, dist;
  int64_t count = 0, total = (dna.length()-i_start) * (dna.length()-i_start+1) / 2 -
                             (dna.length()-i_end-1) * (dna.length()-i_end) / 2;
  for (size_t i = (i_start-1); i < i_end; i++) {
    const String::StringProxy &dna1 = dna[i];
    for (size_t j = (i+1); j < dna.length(); j++) {
      const String::StringProxy &dna2 = dna[j];
      int d = 0;
      for (String::StringProxy::iterator s = dna1.begin(), t = dna2.begin();
           s != dna1.end() && t != dna2.end(); ++s, ++t) {
        if (*s != *t) {
          if (++d > 2)
            break;
        }
      }
      if (d <= 2) { // found!
        st1.push_back(i);
        st2.push_back(j);
        dist.push_back(d);
      }
      if (++count % 1000000 == 0) {
        if (progress) {
          Rcpp::Rcout << "Up to " << count << " of " << total << std::endl;
        }
        Rcpp::checkUserInterrupt();
      }
    }
  }
  if (progress) {
    Rcpp::Rcout << "Total:" << count << std::endl;
  }
  return List::create(Named("gST") = st1,
                      Named("gST2") = st2, 
                      Named("Distance") = dist);
}
