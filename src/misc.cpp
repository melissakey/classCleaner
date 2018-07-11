#include <Rcpp.h>
using namespace Rcpp;

Rcpp::IntegerVector find_string_locs(Rcpp::CharacterVector input_vector, std::string id){
  Rcpp::IntegerVector id_locs; // Setup storage for found IDs

  for(int i =0; i < input_vector.size(); i++) // Loop through input
    if(input_vector[i] == id) // check if input matches target
      id_locs.push_back(i);

    return id_locs; // send locations to R (c++ index shift!)
}
