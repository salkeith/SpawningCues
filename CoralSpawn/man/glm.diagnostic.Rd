\name{glm.diagnostic}
\alias{glm.diagnostic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GLM diagnostics}
\description{
Outputs diagnostic information for generalised linear model.
}
\usage{
glm.diagnostic(formula, error.distribution, dataset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
General formula for regression models e.g., x ~ y
}
  \item{error.distribution}{
Error distribution e.g., binomial, poisson
}
  \item{dataset}{
Dataframe that includes response and predictor variables
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
S. A. Keith
}
\note{
%%  ~~further notes~~
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
function (formula, error.distribution, dataset) 
{
    require(cars)
    fit <- glm(formula, family = error.distribution, data = dataset, 
        na.action = na.fail)
    summary(fit)
    par(mfcol = c(2, 2))
    plot(fit)
    hist(resid(fit))
    print(paste = "VIF =", vif(fit)[, 3])
    X2 <- sum(residuals(fit, type = "pearson")^2)
    phi <- X2/fit$df.residual
    if (phi > 1) 
        print(paste("Model is overdispersed. Phi =", round(phi, 
            3)))
    else print(paste("Model is NOT overdispersed. Phi =", round(phi, 
        3)))
    p.val <- pchisq(phi * fit$df.residual, fit$df.residual, lower = F)
    if (p.val < 0.05) 
        print(paste("Model is overdispersed. p.val =", round(p.val, 
            3)))
    else print(paste("Model is NOT overdispersed. p.val =", round(p.val, 
        3)))
    pp <- sum(resid(fit, type = "pearson")^2)
    gof <- 1 - pchisq(pp, fit$df.resid)
    if (gof < 0.05) 
        print(paste("Model fit is poor. pp =", round(gof, 3)))
    else print(paste("Model fit is OK. pp =", round(gof, 3)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
