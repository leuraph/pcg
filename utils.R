library(Matrix);


#' Solves the linear system ABx = b by solving Ay=b and BX=y.
solve_twice <- function(A, B, b)
{
   y = solve(A, b)
   x = solve(B, y)
   return(x)
}


pcg_ichol <- function(A, b, initial_guess = rep(0, length(b)), maxiter=10000, relative_threshold=1e-5, verbose=FALSE)
{
   L_upper = as(ichol(A), "triangularMatrix");
   L_lower = t(L_upper);

   # initialize
   r_zero = as.vector(b - A %*% initial_guess)
   r_old = r_zero
   x_old = initial_guess
   x_new = initial_guess

   for(i in 1:maxiter)
   {
        z = solve_twice(L_upper, L_lower, r_old)
        rho_old = r_old %*% z
        if (i==1)
        {
            p_new = z
        }
        else
        {
            beta = c(rho_old/rho_old_old)
            p_new = z + beta*p_old
        }
        q = as.vector(A %*% p_new)
        alpha = c(rho_old/(p_new %*% q))
        x_new = x_old + alpha * p_new
        r_new = r_old - alpha * q

        # assignment updates
        rho_old_old = rho_old
        r_old = r_new
        p_old = p_new
        x_old = x_new

        # Check convergence, continue if necessery
        relative_residual = norm(r_new, type="2") / norm(r_zero, type="2")
        if(verbose)
        {
        cat("iteration: ", i, ", relative residual: ",relative_residual)
        }
        if(relative_residual < relative_threshold)
        {
            break;
        }
   }
   return(x_new)
}
