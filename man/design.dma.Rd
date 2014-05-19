\name{design.dma}
\alias{design.dma}
%- Also NEED an '\alis' for EACH other topic documented here.
\title{
Generate a moving average experiment design
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to generate a moving average experiment design given an accession list and a check list. 
}
\usage{
design.dma(enviro = "Eretz", entries = NULL, num.entries = NULL, chk.names = NULL, num.sec.chk = NULL, num.rows = NULL, num.cols = NULL, plot.start = 1001, designID = NULL, annoy = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enviro}{
%%     ~~Describe \code{enviro} here~~
Label for experiemnt
}
  \item{entries}{
%%     ~~Describe \code{entries} here~~
Accession list 
}
  \item{num.entries}{
%%     ~~Describe \code{num.entries} here~~
The number of accession
}
  \item{chk.names}{
%%     ~~Describe \code{chk.names} here~~
Check list
}
  \item{num.sec.chk}{
%%     ~~Describe \code{num.sec.chk} here~~

}
  \item{num.rows}{
%%     ~~Describe \code{num.rows} here~~
Number of rows
}
  \item{num.cols}{
%%     ~~Describe \code{num.cols} here~~
Number of columns
}
  \item{plot.start}{
%%     ~~Describe \code{plot.start} here~~
The start index for plot
}
  \item{designID}{
%%     ~~Describe \code{designID} here~~
Label for accession
}
  \item{annoy}{
%%     ~~Describe \code{annoy} here~~
Check if there exists a output directory in the current working directory
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
%%  ~~who you are~~
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

acc=c(paste("Acc",seq(1,330,1),sep="_")) 
chk=c(paste("chk",seq(1,4,1),sep="_"))

test.ma<-design.dma(entries=acc,chk.names= chk, num.rows=9, num.cols=NULL, num.sec.chk=3, designID="tester1", annoy=T)
#test.ma[[2]]  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Moving Average }
\keyword{Experiment design }% __ONLY ONE__ keyword per line
