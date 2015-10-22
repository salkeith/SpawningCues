\name{lin.quad}
\alias{lin.quad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Linear or quadratic terms?}
\description{
Tests linear and quadratic function regressed against response variable for each predictor variable.
}
\usage{
lin.quad(col.numbers, dataset, error.distribution)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{col.numbers}{
The column numbers in your dataframe that correspond to the predictor variables.
}
  \item{dataset}{
Dataframe that includes response and predictor variables
}
  \item{error.distribution}{
Error distribution e.g., binomial, poisson
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
\author{S.A. Keith}
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
function (col.numbers, dataset, error.distribution) 
{
    for (i in min(col.numbers):max(col.numbers)) {
        diff.AICc <- AICc(glm(Spawning ~ dataset[, i], data = dataset, 
            family = error.distribution)) - AICc(glm(Spawning ~ 
            dataset[, i] + I(dataset[, i]^2), data = dataset, 
            family = error.distribution))
        print(paste(colnames(dataset)[i], " difference in AICc =", 
            diff.AICc))
        if (diff.AICc > 3) 
            print("Use quadratic")
        else print("Use linear")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
