\name{compute.stream}
\alias{compute.stream}
\title{
Calculates point of degeneration j0 into noise of the Idata, applying moderate deviation-based inference
}
\description{
The estimation of \eqn{\hat{j}_0} is achieved via a moderate deviation-based approach. The probability that an estimator, computed from a pilot sample size \eqn{\nu}, exceeds a value \emph{z}, the deviation above \emph{z} is said to be a moderate deviation if its associated probability is polynomially small as a function of \eqn{\nu}, and to be a large deviation if the probability is exponentially small in \eqn{\nu}.  The values of \eqn{z=z_\nu} that are associated with moderate deviations are
\eqn{z_\nu\equiv\bigl(C\,\nu^{-1}\,\log\nu\bigr)^{1/2}}, where \eqn{C>\frac{1}{4}}. The null hypothesis that \eqn{p_k=\frac{1}{2}} for \eqn{\nu} consecutive values of \emph{k}, versus the alternative hypothesis that \eqn{p_k>\frac{1}{2}} for at least one of the values of \emph{k}, is rejected when \eqn{\hat{p}_j^\pm-\frac{1}{2}>z_\nu}.  The probabilities \eqn{\hat{p}_j^+} and \eqn{\hat{p}_j^-} are estimates of \eqn{p_j} computed from the \eqn{\nu} data pairs \eqn{I_\ell} for which \eqn{\ell} lies immediately to the right of \emph{j}, or immediately to the left of \emph{j}, respectively. 

The iterative algorithm consists of an ordered sequence of "test stages" \eqn{s_1, s_2,\ldots}  In stage \eqn{s_k} an integer \eqn{J_{s_k}} is estimated, which is a potential lower bound to \eqn{j_0} (when \eqn{k} is odd), or a potential upper bound to \eqn{j_0} (when \eqn{k} is even).
}
\usage{
compute.stream(Idata, const=0.251, v, r=1.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Idata}{
Input data is a vector of 0s and 1s (see \code{prepare.idata})
}
  \item{const}{
Denotes the constant C of the moderate deviation bound, needs to be larger than 0.25 (default is 0.251)
}
  \item{v}{
Denotes the pilot sample size \eqn{\nu} related to the degree of randomness in the assignments. In each step the noise is estimated from the Idata as probability of 1 within the interval of size \eqn{\nu}, moving from \eqn{J_{s_{k-1}} -r \nu} if \eqn{k} is odd or \eqn{J_{s_{k-1}} +r \nu} if \eqn{k} is even, until convergence or break (see \code{r})
}

  \item{r}{
Denotes a technical constant determining the starting point from which the probability for \eqn{I=1} is estimated in a window of size \code{v} (see \code{v}, default is 1.2)
}
}
\value{
A named list containing:
 \item{j0_est}{Is the estimated index for which the \code{Idata} degenerate into noise}
 \item{k}{\eqn{k=j0_est-1}}
 \item{reason.break}{The reason why the computation has ended - convergence or break condition}
 \item{js}{Is the sequence of estimated \eqn{j_0} in each iteration run, also showing the convergence behaviour}
 \item{v}{Is the preselected value of the parameter \eqn{\nu}}
}
\author{
 Eva Budinska <budinska@iba.muni.cz>, Michael G. Schimek <michael.schimek@medunigraz.at>
}
\seealso{
  \code{\link{prepare.idata}}
}
\examples{
set.seed(465)
myhead <- rbinom(20, 1, 0.8)
mytail <- rbinom(20, 1, 0.5)
mydata <- c(myhead, mytail)
compute.stream(mydata, v=10)	
}

