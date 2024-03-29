% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omics_array.R
\docType{methods}
\name{unionOmics-methods}
\alias{unionOmics-methods}
\alias{unionOmics}
\alias{unionOmics,omics_array,omics_array-method}
\alias{unionOmics,list,ANY-method}
\title{Makes the union between two omics_array objects.}
\usage{
\S4method{unionOmics}{omics_array,omics_array}(M1, M2)
}
\arguments{
\item{M1}{a omics-array or a list of omics-arrays}

\item{M2}{a omics-array or nothing if M1 is a list of omics-arrays}
}
\description{
Makes the union between two omics_array objects.
}
\section{Methods}{
 \describe{

\item{list("signature(M1 = \"omics_array\", M2 = \"omics_array\")")}{
Returns a omics_array object which is the union of M1 and M2.  }

\item{list("signature(M1 = \"list\", M2 = \"ANY\")")}{ Returns a omics_array
object which is the union of the elements of M1.  } }
}

\examples{

if(require(CascadeData)){
data(micro_S, package="CascadeData")
#Create another omicsarray object with 100 genes
Mbis<-M<-as.omics_array(micro_S[1:100,],1:4,6)
#Rename the 100 genes
Mbis@name<-paste(M@name,"bis")
rownames(Mbis@omicsarray) <- Mbis@name
#Union (merge without duplicated names) of the two omicsarrays. 
str(unionOmics(M,Mbis))
}

}
\author{
Bertrand Frederic, Myriam Maumy-Bertrand.
}
\keyword{methods}
