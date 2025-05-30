% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeTree.R
\name{makeTree}
\alias{makeTree}
\title{makeTree}
\usage{
makeTree(data)
}
\arguments{
\item{data}{A dataframe object}
}
\value{
Returns a makeTree object
}
\description{
Assuming a specific structure, create a tree with the following columns:
 from (node label), to (node label), Estimate (+ integer), Total (+ integer),
 and Count (for terminal nodes with marginal counts).
 'from' and 'to' describe the edge for that row of data, where 'Estimate'
 and 'Total' are assumed to come from surveys of size 'Total' (a sample of
 the population at node 'from'), and observe 'Estimate' number of those
 individuals at 'Total' which move to the node described by 'to'.
 'Estimate' and 'Total' columns are used for branching probabilities only.
 'Count' column is NA for rows where 'to' nodes are not leaves; and also
 for all leaves without a marginal count.
 A Population (logical) column is not needed, but can be added if 'Estimate'
 and 'Total' come from population numbers, rather than samples.
 A 'Description' column (string) is also possible to include if particulars
 are desired on the tree diagram.
 'TerminalCount' (binary) will be created for functional purposes, where
 marginal counts are included on leaves.
}
\examples{
data(treeData1)
tree <- makeTree(treeData1)
}
