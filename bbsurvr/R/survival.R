
#########################################################################################################################################################################################################
##Funktion f?r die Berechnung der ?berlebenswahrscheinlichkeit f?r einen bestimmten Zeitpunkt
##########################################################################################################################################################################################################
#' Function to calculate the probability of survival for a point in time
#'
#' @param S xx
#' @param totaltimes xx
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' \dontrun{
#' x<-kmsurv(S,totaltimes)
#' }
#'
#'@export
kmsurv <- function(S, totaltimes) {
  f <- survfit.km(factor(rep(1,nrow(S))), S)
  tt <- c(0, f$totaltimes)
  ss <- c(1, f$surv) # add first point to survival curve
  approx(tt, ss, xout=totaltimes, method='constant', f=0)$y
}
##########################################################################################################################################################################################################
##Funktion fuer das extrahieren des p-Wertes aus der surfdiff funktion
##########################################################################################################################################################################################################
#' Function to exctract p-Value from surfdiff object
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param daten Data.Frame with survival datacontaining min. two columns with a) survival time b) survival status
#' @param group columnsname to be groubed by
#' @param time name of column containing the intervall data. default="time"
#' @param status name of colmn containig status data. default="status"
#'
#' @return a numeric vector representing the p-value
#'
#' @import survival
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' \dontrun{
#' pvalue<-p.value.survdiff(daten=df,group="Treatment",time="duration",status="survivalstatus")
#' }
#'
#'@export
p.value.survdiff<-function(daten,group=NULL, time="time",status="status")
{
  sdf<-survdiff(Surv(daten[[time]],daten[[status]])~daten[,group], data =daten)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  return(p.val)

}

##########################################################################################################################################################################################################
##Funktion f?r eine zum erstellen eines einheitlichen BoxPlots
##########################################################################################################################################################################################################
#' R Function to streamline the generation of (grouped) boxplots
#'
#' @title bbBoxPlot - streamlined BoxPlot
#'
#' @param daten the vector of numeric data
#' @param gruppe if given the factors to be grouped by
#' @param filename if given output will be delivered to file
#' @param ylab label of y axis
#' @param titel title of plot
#'
#' @return a base boxplot
#'
#' @import bbhelper
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' bp<-bbBoxPlot(iris$Petal.Width,gruppe=iris$Species,ylab="Petal Width")
#'
#'@export
bbBoxPlot<- function(daten,gruppe=NULL,filename=NULL,ylab,titel=NULL)
{

  mypalette<-bbhelper::getBBColors(length(levels(gruppe)))
  if (is.null(filename)==FALSE){png(file=paste(c(filename,".png"),collapse = ""),bg="white",res=300,width=1600,height=1600) }
  boxplot(daten ~ gruppe,ylab=ylab,col=mypalette,  horizontal=FALSE,cex.axis=0.7, las=3,show.names=FALSE)
  legend("topright",levels(gruppe),fill=mypalette,cex=0.8)
  if (is.null(titel)==FALSE)
  {
    title(main = list(titel,cex= 0.8, col="#60759B"))
  }
  text(0.2,0,paste(" created by biobits\n",Sys.Date()),col="grey",pos=4,cex=0.6)
  if (is.null(filename)==FALSE){dev.off()}
}

##########################################################################################################################################################################################################
##Funktion f?r das erzeugen eines angepassten Farbthemas
##########################################################################################################################################################################################################
#' R fuction to set a unified layout for base plots
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' \dontrun{
#' bbTheme()
#' }
#'
#'@export
bbTheme <- function() {
  par <- col.whitebg()
  par$strip.background$col <- rep("#60759B", 7)
  par$add.text$col <- "#eeeeaa"
  par$add.text$font <- 2
  par$background$col <- "#ffffff"
  par$superpose.line$lty <- rep(1, 7)
  par$superpose.line$col[1:2] <- c("#880000", "#008800")
  par$superpose.symbol$col[1:2] <- c("#880000", "#008800")
  par
}

##########################################################################################################################################################################################################
## Function to streamline survival plots on base of survminer
##########################################################################################################################################################################################################
#' R Function to streamline the generation of survival plots
#'
#' accepts a data frame containing survival data an
#'
#' @title bbggkmplot - streamlined KM survival curve
#'
#' @param daten the data.frame with survival data
#' @param gruppe optional: the factor the plot has to be grouped by
#' @param time optional: data holding the the time
#' @param status optional: data holding survival status (0/1)
#' @param title the title of the plot
#' @param survtime the time intervall in month to show the survival rate for (e.g. 24 for 2-year survival-rate)
#' @param survtimetext the text to label the survival rate
#' @param risk.table if TRUE the the risk table is ploted beneath the graph. Defalut is true.
#' @param showmedian if true the median value is shown in plot /legend (for more than one group). Default is true.
#' @param median.dig to how many digits the median should be rounded. default is 2
#' @param xmax MAx Value for X-axis
#' @param xlab Label for x-axis (default="Monate (Anzeige bis max. 5 Jahre)")
#' @param cex.lab Fontsize for label (default=1)
#' @param cex.axis Fontsize for axis (default=1)
#' @param watermark if TRUE the biobits watermark will be printet on plot
#' @param ylab the label for y-axis
#' @param logrank if true a logrank test is performed  and the p-value will be printet on the plot (default=FALSE)
#'
#'
#' @return a survival plot
#'
#' @import survminer survival bbhelper dplyr ggplot2 ggpubr
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' \dontrun{
#' data("myeloma")
#' bbggkmplot(daten = myeloma ,time=time
#'           ,gruppe = molecular_group
#'           ,status=event,logrank=T
#'           ,watermark = T,risk.table = T
#'           ,showmedian = T
#'           ,survtime=60
#'           ,survtimetext="5-Y SR")
#'           }
#'
#'@export
bbggkmplot<-function(daten,gruppe=NULL,time=time,status=status,xlab="Time in months"
                     ,cex.lab=1,cex.axis=1,watermark=TRUE,ylab="",title="",survtime=NULL,survtimetext=NULL
                     ,risk.table  = TRUE,logrank=FALSE,xmax=100,showmedian=T,median.dig=2){
  qTIME <- dplyr::enquo(time)                    # Create quosure
  qSTATUS  <- dplyr::enquo(status)               # Create quosure
  qLevels<-1
  qMedian<-"none"
  qLegend<-"right"
  groupnames<-NULL
  subtext<-""

  if(!missing(gruppe)){
    qGRUPPE <- dplyr::enquo(gruppe)
    Sdata<-daten%>%dplyr::select(time= !! dplyr::as_label(qTIME),status= !! dplyr::as_label(qSTATUS),group= !! dplyr::as_label(qGRUPPE))
    Sdata$group<-as.factor(Sdata$group)
    fit <- survival::survfit(survival::Surv(time=time, event=status) ~ group, data=Sdata)
    groupnames<-levels(droplevels(Sdata$group))
    qLevels<-length(groupnames)

    if (showmedian){
      gMedian<-as.vector(survminer::surv_median(fit)["median"])
      medianText<-paste0("\nMedian = ", as.vector(unlist(round(gMedian,median.dig))))
      groupnames<-paste0(groupnames,medianText)
    }
  }else{
    Sdata<-daten%>%dplyr::select(time= !! dplyr::as_label(qTIME),status= !! dplyr::as_label(qSTATUS))
    fit <- survival::survfit(survival::Surv(time=time, event=status)~1 , data=Sdata)
    groupnames<-""
    if(showmedian){
      qMedian<-"hv"
      gMedian<-as.vector(survminer::surv_median(fit)["median"])
      subtext<-paste("Median=",round(gMedian,median.dig))
    }
    qLegend<-"none"
  }

  mypalette<-bbhelper::getBBColors(qLevels)

  #survivalprobability at timepoint

  if(!is.null(survtime)){
    survtimepoint<-formatC((summary(fit, times=survtime,extend=TRUE)$surv),digits=2,format="f")
    # Wenn wir nur eine Gruppe haben Angabe als "caption"
    if(missing(gruppe)){
      subtext<-paste0(subtext,paste("\n",survtimetext," = ",survtimepoint))
    }else{
      groupnames<-paste0(groupnames,paste("\n",survtimetext," = ",survtimepoint,"\n"))
    }
  }
  survplot<-survminer::ggsurvplot(
    fit, # fitted survfit object
    data=Sdata,
    risk.table  = risk.table, # include risk table?
    conf.int    = TRUE, # add confidence intervals?
    pval        = logrank, # add p-value to the plot?
    pval.size=3,
    title=title,
    break.time.by = 12, #default because we mostly use months
    legend=qLegend,
    xlim = c(0, xmax),
    xlab = xlab,
    caption  =subtext,
    show.legend = FALSE,
    ylab=ylab,
    risk.table.col = "strata",
    risk.table.fontsize=3,
    palette=mypalette,
    ggtheme = theme_light(),
    legend.title = "",
    surv.median.line = qMedian, # median survival in plot
    conf.int.style = "step",  # customize style of confidence intervals
    risk.table.y.text = FALSE,# show bars instead of names in text annotations
    # in legend of risk table.
    legend.labs = groupnames,
    font.x       = c(8), #  font for the x axis
    font.legend = c(8),
    font.caption =c(7)
  )
  if (watermark){
    #
    survplot$plot<- survplot$plot + ggplot2::annotate("text", x = Inf, y = -Inf,
                                                      label = paste(" created by biobits\n",Sys.Date()),
                                                      hjust=1.1, vjust=-0.2, col="grey", cex=2.5, alpha = 0.8)
  }
  if(risk.table){
    survplot$table <- ggpubr::ggpar(
      survplot$table,
      font.title    = c(8),
      font.x        = c(8),
      font.xtickslab = c(8)
    )
  }
  survplot
}


##########################################################################################################################################################################################################
##Funktion zum erstellen eines einheitlichen KM-Plots
##########################################################################################################################################################################################################
#' R deprecated ! Function to streamline the generation of survival plots
#'
#' accepts a data frame containing survival data an
#' @title bbKMplot - streamlined KM survival curve
#'
#' @param daten the data.frame with survival data
#' @param gruppe optional: the factor the plot has to be grouped by
#' @param titel the title of the plot
#' @param filename the filename if the plot should be saved as file
#' @param survtime the time intervall in month to show the survival rate for (e.g. 24 for 2-year survival-rate)
#' @param survtimetext the text to label the survival rate
#' @param legendeout if TRUE the legend will be printet outside the plot margins
#' @param file.out if TRUE the plot will be saved to file
#' @param legfontsize fontsize of legend (default=1)
#' @param subtext if needed a subtext to be placed beneath the plot
#' @param xmax MAx Value for X-axis
#' @param xlab Label for x-axis (default="Monate (Anzeige bis max. 5 Jahre)")
#' @param cex.lab Fontsize for label (default=1)
#' @param cex.axis Fontsize for axis (default=1)
#' @param watermark if TRUE the biobits watermark will be printet on plot
#' @param ylab the label for y-axis
#' @param logrank if true a logrank test is performed for two(!) survival curves and the result will be printet on the plot (default=FALSE)
#'
#'
#' @return a survival plot
#'
#' @import survival bbhelper
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' \dontrun{
#' p<-bbKMPlot(daten,gruppe=NULL,titel,filename,survtime=NULL,survtimetext=NULL,legendeout=NULL,file.out=TRUE,legfontsize=NULL,subtext="",
#' xmax=61.25,xlab="Beobachtungszeit in Monaten",cex.lab=1,cex.axis=1,watermark=TRUE,ylab="",logrank=FALSE)
#' }
#'
#'@export
bbKMPlot <- function (daten,gruppe=NULL,titel,filename,survtime=NULL,survtimetext=NULL,legendeout=NULL,file.out=TRUE,legfontsize=NULL,subtext="",
                      xmax=61.25,xlab="Beobachtungszeit in Monaten",cex.lab=1,cex.axis=1,watermark=TRUE,ylab="",logrank=FALSE)
{
  #status als numeric definieren
  daten$status<-as.numeric(daten$status)
  blwd<-2
  #Schriftgroesse der Legende
  if (is.null(legfontsize)==TRUE){leg.cex<-0.6} else {leg.cex<-legfontsize}

  if ((is.null(gruppe)==TRUE) || ((length(levels(as.factor(daten[[gruppe]])))==1)))
  {
    #aucount <- length(daten$time)
    if(file.out==TRUE)
    {png(file=paste(c(filename,".png"),collapse = ""),bg="white",res=300,width=1600,height=1600)}
    fitg<-survival::survfit(survival::Surv(time, status)~1 ,data = daten,type="kaplan-meier")
    survmedian<-(summary(fitg)$table["median"][[1]])
    aucount<-formatC((summary(fitg)$table["records"][[1]]),digits=0,format="f")
    survmedian_text<-formatC(survmedian,digits = 1,format = "f")

    plot(fitg,xlab=xlab,ylab=ylab,xmax=xmax,cex.lab=cex.lab,cex.axis=cex.axis,lwd = blwd)
    if(watermark==TRUE){text(0,0,paste(" created by biobits\n",Sys.Date()),col="grey",pos=4,cex=0.6)}

    # median einzeichnen
    segments(0, 0.5, survmedian,0.5, col = "#a6202a", lty = 2, lwd = 1.5)
    segments(survmedian, 0.5, survmedian,0, col = "#a6202a", lty = 2, lwd = 1.5)
    text(survmedian,-0.02,paste(c("Median=\n\n",survmedian_text)),cex=0.6,col = "#a6202a")
    if (is.null(survtime)==FALSE)
    {
      survtimepoint<-formatC((summary(fitg, times=survtime)$surv),digits=2,format="f")
      title(main = list(c(titel,"\n",paste("n =",aucount)),cex= 0.8, col="#60759B"), sub=paste(survtimetext," = ",survtimepoint))
    }
    else{
      title(main = list(c(titel,"\n",paste("n =",aucount)),cex= 0.8, col="#60759B"))
    }

    # ausgabe ende
    if(file.out==TRUE)
    {dev.off()}
  }
  else
  {
    daten[,gruppe]<-as.factor(daten[[gruppe]])
    Acount <- length(daten$time)
    fits<-survival::survfit(survival::Surv(time,status) ~daten[[gruppe]], data =daten)
    gruppen_list  <- levels(daten[[gruppe]])
    MEDIANE <-list()
    ANZAHL<-list()
    SURVTIMEPOINT<-list()
    {if (length(gruppen_list)==1)
    {
      MEDIANE[gruppen_list[c]]<-formatC((summary(fits)$table["median"][[c]]),digits=1,format="f")
      ANZAHL[gruppen_list[c]]<-formatC((summary(fits)$table["records"][[c]]),digits=0,format="f")
      if (is.null(survtime)==FALSE)
      {
        SURVTIMEPOINT[gruppen_list[c]]<-formatC((summary(fits, times=survtime)$surv[c]),digits=2,format="f")
      }
    }
      else{
        for (c in 1:length(gruppen_list))
        {
          MEDIANE[gruppen_list[c]]<-formatC((summary(fits)$table[,"median"][[c]]),digits=1,format="f")
          ANZAHL[gruppen_list[c]]<-formatC((summary(fits)$table[,"records"][[c]]),digits=0,format="f")
          if (is.null(survtime)==FALSE)
          {
            #MIT Errorhandling, da bei zu gro?em Zeitraum ein fehler geworfen wird
            res<-try((SURVTIMEPOINT[gruppen_list[c]]<-formatC((summary(fits, times=survtime)$surv[c]),digits=2,format="f")),silent=TRUE)
            if(class(res) == "try-error"){SURVTIMEPOINT[gruppen_list[c]]<-NA}
          }
        }
      }}
    mypalette<-bbhelper::getBBColors(length(MEDIANE))

    # Legende ausserhalb des Plotbereichs

    if(is.null(legendeout)==FALSE)
    {
      if(file.out==TRUE)
      {png(file=paste(c(filename,".png"),collapse = ""),bg="white",res=300,width=1800,height=1600)}
      par(mar=c(5, 4, 4, 10)+.1,xpd=TRUE)#oma=c(0,0,0,2)plt=c(1,1,1,1.5),
      plot(fits,xlab=xlab,ylab=ylab,xmax=xmax,col=mypalette,cex.lab=cex.lab,cex.axis=cex.axis,lwd = blwd)
      tmp.u <- par('usr')
      leg.pos<-list(x=tmp.u[2], y=tmp.u[4], xjust=0, yjust=0,outer=TRUE)
    }
    else
    {
      leg.pos<-"topright"

      if(file.out==TRUE)
      {png(file=paste(c(filename,".png"),collapse = ""),bg="white",res=300,width=1600,height=1600)}
      plot(fits,xlab=xlab,ylab=ylab,xmax=xmax,col=mypalette,cex.lab=cex.lab,cex.axis=cex.axis,lwd = blwd)
    }



    if(watermark==TRUE){text(0,0,paste(" created by biobits\n",Sys.Date()),col="grey",pos=4,cex=0.6)}
    title(main = list(paste(titel,sep="" ),cex= 0.8, col="#60759B"),sub=list(paste(subtext,sep="" ),cex= 0.6))


    if (is.null(survtime)==FALSE)
    {
      legend(leg.pos,paste(names(MEDIANE),"\nMedian =",MEDIANE," n =",ANZAHL,"\n",survtimetext,"=",SURVTIMEPOINT,"\n"),lty=1,col=mypalette, adj = c(0, .6),cex= leg.cex,bty="n")
    }
    else
    {
      legend(leg.pos,paste(names(MEDIANE),"\nMedian =",MEDIANE," n =",ANZAHL,"\n"),lty=1,col=mypalette, adj = c(0, .6),cex= leg.cex,bty="n")
    }
    #Logrank test
    lograng_txt<-""
    if (logrank==TRUE){
      cdiff<-survdiff(Surv(time,status) ~daten[,gruppe], data =daten)
      cp.diff<-pchisq(cdiff$chisq, df=1, lower=FALSE)
      if(cp.diff<0.00001)
      {lograng_txt="Logrank: p < 0.00001"}
      else
      {if (cp.diff>0.001)
      {lograng_txt<-paste("Logrank: p = ",formatC(cp.diff,digits = 3,format = "f"))}
        else
        {lograng_txt<-paste("Logrank: p = ",formatC(cp.diff,digits = 5,format = "f"))}
      }

      pypos<-(par('xaxp')[2]/par('xaxp')[3])/2
      pxpos<-0.12
      text(pypos,pxpos,lograng_txt,pos=4,cex=leg.cex)
    }

    if(file.out==TRUE)
    {dev.off()}
  }
}






