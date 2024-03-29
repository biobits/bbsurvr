% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival.R
\name{bbggkmplot}
\alias{bbggkmplot}
\title{bbggkmplot - streamlined KM survival curve}
\usage{
bbggkmplot(
  daten,
  gruppe = NULL,
  time = time,
  status = status,
  xlab = "Time in months",
  cex.lab = 1,
  cex.axis = 1,
  watermark = TRUE,
  ylab = "",
  title = "",
  survtime = NULL,
  survtimetext = NULL,
  risk.table = TRUE,
  logrank = FALSE,
  xmax = 100,
  showmedian = T,
  median.dig = 2,
  ggtheme = theme_light(),
  conf.int = TRUE,
  legend = "right"
)
}
\arguments{
\item{daten}{the data.frame with survival data}

\item{gruppe}{optional: the factor the plot has to be grouped by}

\item{time}{optional: data holding the the time}

\item{status}{optional: data holding survival status (0/1)}

\item{xlab}{Label for x-axis (default="Monate (Anzeige bis max. 5 Jahre)")}

\item{cex.lab}{Fontsize for label (default=1)}

\item{cex.axis}{Fontsize for axis (default=1)}

\item{watermark}{if TRUE the biobits watermark will be printet on plot}

\item{ylab}{the label for y-axis}

\item{title}{the title of the plot}

\item{survtime}{the time intervall in month to show the survival rate for (e.g. 24 for 2-year survival-rate)}

\item{survtimetext}{the text to label the survival rate}

\item{risk.table}{if TRUE the the risk table is ploted beneath the graph. Defalut is true.}

\item{logrank}{if true a logrank test is performed  and the p-value will be printet on the plot (default=FALSE)}

\item{xmax}{MAx Value for X-axis}

\item{showmedian}{if true the median value is shown in plot /legend (for more than one group). Default is true.}

\item{median.dig}{to how many digits the median should be rounded. default is 2}

\item{ggtheme}{function, ggplot2 theme name. Default value is theme_light.
Allowed values include ggplot2 official themes: see theme}

\item{conf.int}{logical value. If TRUE, plots confidence interval}

\item{legend}{character specifying legend position. Allowed values are one of c("top", "bottom", "left", "right", "none").}
}
\value{
a survival plot
}
\description{
R Function to streamline the generation of survival plots
}
\details{
accepts a data frame containing survival data an
}
\examples{
\dontrun{
data("myeloma")
bbggkmplot(daten = myeloma ,time=time
          ,gruppe = molecular_group
          ,status=event,logrank=T
          ,watermark = T,risk.table = T
          ,showmedian = T
          ,survtime=60
          ,survtimetext="5-Y SR")
          }

}
\author{
Stefan Bartels, \email{email@biobits.eu}
}
