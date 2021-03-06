\name{select.average}
\alias{select.average}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model selection and averaging
}
\description{
Models within 3 delta AICc of the best model are averaged. Full coefficients are best.
}
\usage{
select.average(mod.type = NULL, formula, error.distribution, dataset, msubset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod.type}{
Whether it is a "mixed" model.
}
  \item{formula}{
General formula for regression models e.g., x ~ y. Add random effects if it is a mixed model.
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
\item{dredge.res }{Result of the dredge function - list of all possible models with AICc values}
\item{averaged.model }{Averaged model generated by MuMIn model.average function}
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
function (mod.type = NULL, formula, error.distribution, dataset, 
    msubset) 
{
    require(lme4)
    if (mod.type == "mixed") 
        fit <- glmer(formula, family = error.distribution, data = dataset, 
            na.action = na.fail)
    else fit <- glm(formula, family = error.distribution, data = dataset, 
        na.action = na.fail)
    print("Patience... the next step can be slow...")
    dmod <- dredge(fit, subset = msubset)
    fit.avg <- summary(model.avg(dmod, subset = delta < 3))
    print(fit.avg)
    return(list(dredge.res = dmod, averaged.model = fit.avg))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
