// LICENSE: GPL (>= 2)

#include <Rcpp.h>
#include <list>



// [[Rcpp::export]]
Rcpp::NumericVector
not_able_knots(Rcpp::NumericMatrix all, int delta) {
    std::list<int> notAble;
    int k = 0;	
    for(int i = 0; i < all.ncol(); i++) {
	for (int j = 0; j < all.nrow() - 1; j++) {
	    if (all(j + 1, i) - all(j, i) <= delta) {
		notAble.push_back(i+1);
		break;
	    }
	}
    }
    Rcpp::NumericVector RnotAble = Rcpp::wrap(notAble);
    return RnotAble;
}
