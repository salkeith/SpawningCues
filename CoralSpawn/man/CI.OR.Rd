\name{CI.OR}
\alias{CI.OR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CI.OR}
\description{
Calculate 95 percent confidence intervals and odds ratios for regression coefficients.
}
\usage{
CI.OR(fit.avg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.avg}{
Output from averaged model created using MuMIn package.
}
}
\details{
Odds ratios are the exponentiated regression coefficients.
Odds ratio: change in "odds" of a positive result with a 1 unit change in the predictor variable.
}
\value{
\item{Confidence intervals }{Estimated coefficient and 95 percent confidence intervals}
\item{Odds ratios }{Odds raitos and 95 percent confidence intervals}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
S. A. Keith
}
\note{
-
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fit.avg) 
{
    cimod <- confint(fit.avg)
    comod <- coef(fit.avg)
    CImod <- cbind(comod, cimod)
    CImod
    CImod[!(CImod[, 2] < 0 & CImod[, 3] > 0), ]
    CIORmod <- exp(CImod)
    print(CImod)
    return(list(CI95 = data.frame(CImod), odds.ratios = data.frame(CIORmod)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
