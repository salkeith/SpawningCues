\name{mixed.mod.output}
\alias{mixed.mod.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixed effects model output from glmer}
\description{
Run a generalised linear mixed effects model using lme4 and view the summary. Diagnostic information also supplied.
}
\usage{
mixed.mod.output(formula, error.distribution, dataset, random)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
General formula for mixed effect regression models e.g., x ~ y + (1 | z).
}
  \item{error.distribution}{
Error distribution e.g., binomial, poisson
}
  \item{dataset}{
Dataframe that includes response, predictor variables and random effects.
}
  \item{random}{
String vector of the names of your random effects as written in the dataframe.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Printed outputs for summary, dotplots and variance partition coefficients.
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
function (formula, error.distribution, dataset, random) 
{
    require(lme4)
    fit <- glmer(formula, family = error.distribution, data = dataset, 
        na.action = na.fail)
    print(summary(fit))
    dp <- dotplot(ranef(fit, condVar = TRUE, whichel = random[1]))
    print(dp)
    if (length(random) > 1) 
        dp2 <- dotplot(ranef(fit, condVar = TRUE, whichel = random[2]))
    try(print(dp2))
    random.vars <- as.data.frame(VarCorr(fit))
    month.vpc <- random.vars[which(random.vars[, 1] == random[1]), 
        4]
    month.vpc/(month.vpc + 3.29)
    site.vpc <- random.vars[which(random.vars[, 1] == random[2]), 
        4]
    site.vpc/(site.vpc + 3.29)
    print(paste(random[1], " VPC = ", round(month.vpc, 3), ". ", 
        random[1], " accounts for ", round(month.vpc, 5) * 100, 
        " of the variance in the data", sep = ""))
    print(paste(random[2], " VPC = ", round(site.vpc, 3), ". ", 
        random[2], " accounts for ", round(site.vpc, 5) * 100, 
        " of the variance in the data", sep = ""))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
