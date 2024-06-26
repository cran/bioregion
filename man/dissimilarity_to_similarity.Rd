% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity_dissimilarity_conversion.R
\name{dissimilarity_to_similarity}
\alias{dissimilarity_to_similarity}
\title{Convert dissimilarity metrics to similarity metrics}
\usage{
dissimilarity_to_similarity(dissimilarity, include_formula = TRUE)
}
\arguments{
\item{dissimilarity}{the output object from \code{\link[=dissimilarity]{dissimilarity()}} or
\code{\link[=similarity_to_dissimilarity]{similarity_to_dissimilarity()}}.}

\item{include_formula}{a \code{boolean} indicating if the metrics based on your
own formula(s) should be converted (see Details). This argument is set to
\code{TRUE} by default.}
}
\value{
A \code{data.frame} with additional class
\code{bioregion.pairwise.metric}, providing similarity
metric(s) between each pair of sites based on a dissimilarity object.
}
\description{
This function converts a \code{data.frame} of dissimilarity metrics (beta diversity)
between sites to similarity metrics.
}
\note{
\loadmathjax
The behavior of this function changes depending on column names. Columns
\code{Site1} and \code{Site2} are copied identically. If there are columns called
\code{a}, \code{b}, \code{c}, \code{A}, \code{B}, \code{C} they will also be copied identically. If there
are columns based on your own formula (argument \code{formula} in
\code{\link[=dissimilarity]{dissimilarity()}}) or not in the original list of dissimilarity metrics
(argument \code{metrics} in \code{\link[=dissimilarity]{dissimilarity()}}) and if the argument
\code{include_formula} is set to \code{FALSE}, they will also be copied identically.
Otherwise there are going to be converted like they other columns (default
behavior).

If a column is called \code{Euclidean}, the similarity will be calculated based
on the following formula:

\mjeqn{Euclidean similarity = 1 / (1 - Euclidean distance)}{Euclidean similarity = 1 / (1 - Euclidean distance)}

Otherwise, all other columns will be transformed into dissimilarity with the
following formula:

\mjeqn{similarity = 1 - dissimilarity}{similarity = 1 - dissimilarity}
}
\examples{
comat <- matrix(sample(0:1000, size = 50, replace = TRUE,
prob = 1 / 1:1001), 5, 10)
rownames(comat) <- paste0("Site", 1:5)
colnames(comat) <- paste0("Species", 1:10)

dissimil <- dissimilarity(comat, metric = "all")
dissimil

similarity <- dissimilarity_to_similarity(dissimil)
similarity

}
\seealso{
\code{\link[=similarity_to_dissimilarity]{similarity_to_dissimilarity()}} \code{\link[=similarity]{similarity()}} \code{\link[=dissimilarity]{dissimilarity()}}
}
\author{
Maxime Lenormand (\email{maxime.lenormand@inrae.fr}),
Boris Leroy (\email{leroy.boris@gmail.com}) and
Pierre Denelle (\email{pierre.denelle@gmail.com})
}
