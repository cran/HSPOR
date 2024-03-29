% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HKSPOR_DynProg.R
\name{HKSPOR_DynProg}
\alias{HKSPOR_DynProg}
\title{Inference method for any number K of regimes using dynamic programming}
\usage{
HKSPOR_DynProg(X, Y, deg, K, constraint = 1, smoothing = TRUE,
  verbose = FALSE, plotG = TRUE)
}
\arguments{
\item{X}{A numerical vector corresponding to the explanatory variable. X must be sorted in ascending order
if this is not the case, X will be sorted in the function and the corresponding permutation will be applied to Y. The
user will be notified by a warning message. In addition, if X contains NAs, they will be deleted from the data and the user will be notified by a warning message.
Finally, if X contains duplicate data, the excess data will be deleted and the user will be notified by a warning message.}

\item{Y}{A numerical vector corresponding to the variable to be explain. It should contain at least two regimes that
could be modelled by polynomials. In addition, if Y contains NAs they will be deleted from the data and the
user will be notified by a warning message. Finally, if X contains dupplicate data, the excess data will be deleted and
the value of the remaining Y will become the average of the Ys, calculated for this value of X.}

\item{deg}{The degree of the polynomials. The size of X and Y must be greater than K(deg+2) + K.}

\item{K}{The number of regimes. The size of X and Y must be greater than K(deg+2) + K.}

\item{constraint}{Number that determines the regularity assumption that is applied for the parameters estimation.
By default, the variable is set to 1, i. e. the parameters estimation is done under continuity constraint.
If the variable is 0 or 2, the estimation of the parameters will be done without assumption of regularity
(constraint = 0) or under assumption of differentiability (constraint = 2). Warning, if the differentiability
assumption is not verified by the model, it is preferable not to use it to estimate the model parameters.
In addition, in this dynamic programming method, to ensure that the number of constraints is not greater
that the number of parameters to be estimated, the degree of the polynomials must be at least equal to
3 to be able to use the differentiability assumption.}

\item{smoothing}{A Boolean. If TRUE (default), the method will estimate the parameters of a piecewise polynomial regression
model with latent variable by maximizing the log-likelihood weighted by the probability of being in the
latent variable regime. If FALSE, the method will estimate the parameters of the piecewise polynomial regression
model.}

\item{verbose}{A Boolean. If FALSE (default) the HKSPOR_Dynprog function will return only one dataframe
containing the parameter estimates obtained for a model at K regimes. If TRUE, the function will return
all the results obtained for a model with 1 regime up to K regimes.}

\item{plotG}{A Boolean. If TRUE (default) the estimation results obtained by the HKSPOR_DynProg function are plotted.}
}
\value{
One or more dataframes depend on the verbose value. If verbose = False, the output table will
contain the estimated parameters of the polynomial regression model at K regimes: jump times, polynomial
coefficients and variances of K regimes. If verbose = True then there will be K dataframes in output.
Each table will contain the results of the estimated parameters obtained for each value of k (k=1,...,k=K).
If plotG = TRUE, the data (X,Y) and the estimated model(s) will
be plotted.
}
\description{
HKSPOR_DynProg is an inference method implemented in the form of a Bellman algorithm that estimates, under the assumption of regularity,
the parameters of a polynomial regression model for a number K of regimes given by the user..
}
\examples{
#generated data with three regimes
set.seed(1)
xgrid1 = seq(0,10,length.out=6)
xgrid2 = seq(10.2,20,length.out=6)
ygrid1 = xgrid1^2-xgrid1+1+ rnorm(length(xgrid1),0,4)
ygrid2 = rep(91,length(xgrid2))+ rnorm(length(xgrid2),0,4)
datX = c(xgrid1,xgrid2)
datY = c(ygrid1,ygrid2)

#Inference of a polynomial regression model with two regimes (K=2) on these data.
#The degree of the polynomials is fixed to 2 and the parameters are estimated
#under continuity constraint.
HKSPOR_DynProg(datX,datY,2,2)

\donttest{
set.seed(2)
xgrid1 = seq(0,10,by=0.2)
xgrid2 = seq(10.2,20,by=0.2)
xgrid3 = seq(20.2,30,by=0.2)
ygrid1 = xgrid1^2-xgrid1+1+ rnorm(length(xgrid1),0,3)
ygrid2 = rep(91,length(xgrid2))+ rnorm(length(xgrid2),0,3)
ygrid3 = -10*xgrid3+300+rnorm(length(xgrid3),0,3)
datX = c(xgrid1,xgrid2,xgrid3)
datY = c(ygrid1,ygrid2,ygrid3)

#Inference of a polynomial regression model with three (K=3) regimes on these data.
#The degree of the polynomials is fixed to 2 and the parameters are estimated
#under continuity constraint.
HKSPOR_DynProg(datX,datY,2,3)
#Executed time : 3.658121 mins (intel core i7 processor)
}
}
