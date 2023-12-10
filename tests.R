library(Matrix);
library(GPvecchia);
source('utils.R')


test_back_and_forwardsolve = function(verbose = FALSE)
{
   dim = 100

   A = rsparsematrix(dim, dim, 0.1) 
   ATA = A%*%t(A) + sparseMatrix(1:dim, 1:dim, x=1.0);
   L_upper = as(ichol(ATA), "triangularMatrix");
   L_lower = t(L_upper);

   b = runif(n=dim, min=0.01, max=1.0)

   z_actual = back_and_forwardsolve(L_upper, L_lower, b)
   z_expected = solve(L_upper %*% L_lower, b)

   cat("back_and_forwardsolve == solve: ", all.equal(z_actual, z_expected, tolerance=0.001), "\n")
}


test_pcg_ichol = function(verbose = FALSE)
{
   dim = 1000

   A = rsparsematrix(dim, dim, nnz = 10*dim)

   diag_corrector = sparseMatrix(1:dim, 1:dim, x=1.0)

   ATA <- (A %*% t(A)) + diag_corrector

   b = runif(n=dim, min=0.01, max=1.0)

   y = pcg_ichol(ATA, b, verbose=verbose)
   y_exact = solve(ATA, b)

   cat("pcg_ichol == solve: ", all.equal(y_exact, y, tolerance=0.001), "\n")
}

test_pcg_ichol()
test_back_and_forwardsolve()
