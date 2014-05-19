\name{MAD-package}
\alias{MAD-package}
\alias{MAD}
\docType{package}
\title{MAD}
\description{
%More about what it does (maybe more than one line)
%~~ DA concise (1-5 lines) description of the package ~~
This package is used to generate a moving average experiemental design for plant breeding. 
}
\details{
\tabular{ll}{
Package: \tab MAD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-17\cr
LiMcense: 
%What license is it under?
\tab GPL \cr
}
%~~ An overview of how to use the package, including the most important functions ~~
Curently, this package includes a function design.dma. To use this function, you need to give a accession list and a check list in your experiment. You also need to sed up the number or rows and number of blocks in your design. 
}
\author{
%Who wrote it#
Tyler Tiede in University of Minnesota wrote design.dma fucntion \cr
Aimin Yan in Cornell University reframed Tyler's reframed code to make it into a R package \cr
Maintainer: Aimin Yan <ay247@cornell.edu>
%~~ The author and/or maintainer of the package ~~
}
\references{
%~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ 
%package 
Moving Average Design
}

\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%\code{\link{agricolae}}
}
\examples{
%~~ simple examples of the most important functions ~~
#simple example to use design.dma
acc=c(paste("Acc",seq(1,330,1),sep="_")) 
chk=c(paste("chk",seq(1,4,1),sep="_"))
test.ma<-design.dma(entries=acc,chk.names= chk, num.rows=9, num.cols=NULL, num.sec.chk=3, designID="tester1", annoy=T)
#test.ma
}



