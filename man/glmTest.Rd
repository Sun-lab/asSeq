\name{glmTest}
\alias{glmTest}
\title{
  score test for generalized linear model
}
\description{
 Given the NULL model fitted by function glmFit, test for additional covariate. 
}
\usage{
glmTest(gm, Z)
}
\arguments{
  \item{gm}{
  	A model fitted by function glmFit, of class "glmFit"
	}
  \item{Z}{
	a set of additional covariate. Z must be a matrix. Each column corresponds 
	to a covariate and each row corresponds to an observation. 
}
}
\value{
 A list of chi-square test statistic, degree of freedom and p-value. 
}
\author{
Wei Sun weisun@email.unc.edu
}
\note{
Part of codes is adopted from R/CNVtools
}

\keyword{ methods }

