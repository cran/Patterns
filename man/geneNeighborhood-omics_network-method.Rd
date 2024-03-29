% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omics_network.R
\name{geneNeighborhood,omics_network-method}
\alias{geneNeighborhood,omics_network-method}
\alias{geneNeighborhood}
\alias{geneNeighborhood-methods}
\title{Find the neighborhood of a set of nodes.}
\usage{
\S4method{geneNeighborhood}{omics_network}(
  net,
  targets,
  nv = 0,
  order = length(net@time_pt) - 1,
  label_v = NULL,
  ini = NULL,
  frame.color = "white",
  label.hub = FALSE,
  graph = TRUE,
  names = F
)
}
\arguments{
\item{net}{a omics_network object}

\item{targets}{a vector containing the set of nodes}

\item{nv}{the level of cutoff. Defaut to 0.}

\item{order}{of the neighborhood. Defaut to `length(net@time_pt)-1`.}

\item{label_v}{vector defining the vertex labels.}

\item{ini}{using the ``position'' function, you can
fix the position of the nodes.}

\item{frame.color}{color of the frames.}

\item{label.hub}{logical ; if TRUE only the hubs are labeled.}

\item{graph}{plot graph of the network. Defaults to `TRUE`.}

\item{names}{return names of the neighbors. Defaults to `FALSE`.}
}
\value{
The neighborhood of the targeted genes.
}
\description{
Find the neighborhood of a set of nodes.
}
\examples{

data(Selection)
data(infos)
#Find probesets for EGR1
pbst_EGR1 = infos[infos$hgnc_symbol=="EGR1", "affy_hg_u133_plus_2"]

gene_IDs = infos[match(Selection@name, infos$affy_hg_u133_plus_), "hgnc_symbol"]

data(network)
#A nv value can chosen using the cutoff function
nv=.11 
EGR1<-which(is.element(Selection@name,pbst_EGR1))
P<-position(network,nv=nv)

geneNeighborhood(network,targets=EGR1,nv=nv,ini=P,
label_v=gene_IDs)
          
}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}
\keyword{methods}
