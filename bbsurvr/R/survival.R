
#########################################################################################################################################################################################################
##Funktion f?r die Berechnung der ?berlebenswahrscheinlichkeit f?r einen bestimmten Zeitpunkt
##########################################################################################################################################################################################################
#' Function to calculate the probability of survival for a point in time
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param S
#' @param totaltimes
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' kmsurv(S,M)
#'
#'@export
kmsurv <- function(S, totaltimes) {
  f <- survfit.km(factor(rep(1,nrow(S))), S)
  tt <- c(0, f$totaltimes)
  ss <- c(1, f$surv) # add first point to survival curve
  approx(tt, ss, xout=totaltimes, method='constant', f=0)$y
}
##########################################################################################################################################################################################################
##Funktion f?r das ektrahieren des p-Wertes aus der surfdiff funktion
##########################################################################################################################################################################################################
#' Function to exctract p-Value via surfdiff object
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param list of vectors of identical length
#'
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' coalesce(c(NA,NA,53))
#'
#'@export
p.value.survdiff<-function(daten,group=NULL, time="time",status="status")
{
  sdf<-survdiff(Surv(daten[[time]],daten[[status]])~daten[,group], data =daten)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  return(p.val)

}

##########################################################################################################################################################################################################
##Funktion zum erstellen eines einheitlichen KM-Plots
##
##########################################################################################################################################################################################################
#' R coalesce Function
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param list of vectors of identical length
#'
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' coalesce(c(NA,NA,53))
#'
#'@export
bbKMPlot <- function (daten,gruppe=NULL,titel,filename,survtime=NULL,survtimetext=NULL,legendeout=NULL,file.out=TRUE,legfontsize=NULL,subtext="",
                       xmax=61.25,xlab="Monate (Anzeige bis max. 5 Jahre)",cex.lab=1,cex.axis=1,watermark=TRUE,ylab="",logrank=FALSE)
{
  #status als numeric definieren
  daten$status<-as.numeric(daten$status)
  blwd<-2
  #Schriftgröße der Legende
  if (is.null(legfontsize)==TRUE){leg.cex<-0.6} else {leg.cex<-legfontsize}

  if ((is.null(gruppe)==TRUE) || ((length(levels(as.factor(daten[[gruppe]])))==1)))
  {
    #aucount <- length(daten$time)
    if(file.out==TRUE)
    {png(file=paste(c(filename,".png"),collapse = ""),bg="white",res=300,width=1600,height=1600)}
    fitg<-survfit(Surv(time, status)~1 ,data = daten,type="kaplan-meier")
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
    fits<-survfit(Surv(time,status) ~daten[[gruppe]], data =daten)
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
    mypalette<-getBBColors(length(MEDIANE))

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
##########################################################################################################################################################################################################
##Funktion f?r eine zum erstellen eines einheitlichen BoxPlots
##########################################################################################################################################################################################################
#' R coalesce Function
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param list of vectors of identical length
#'
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' coalesce(c(NA,NA,53))
#'
#'@export
bbBoxPlot<- function(daten,gruppe=NULL,filename=NULL,ylab,titel=NULL)
{

  mypalette<-getBBColors(length(levels(gruppe)))
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
##Funktion f?r das erzeugen eines angepassten RGB-GmbH Farbthemas
##########################################################################################################################################################################################################
#' R coalesce Function
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param list of vectors of identical length
#'
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' coalesce(c(NA,NA,53))
#'
#'@export
rgb.theme <- function() {
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
##Funktion für die Darstellung des relativen 1-5-Jahresüberlebn nach Gruppe
##########################################################################################################################################################################################################
#' R coalesce Function
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param list of vectors of identical length
#'
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' coalesce(c(NA,NA,53))
#'
#'@export
bbRelSurvPlot<- function(x,gruppe=NULL,ylab=NULL,xlab=NULL,titel=NULL,jahrstart=NULL,jahrende=NULL)
{


  j_start<-1997
  j_end<-2012
  if (is.null(jahrstart)==FALSE){j_start<-jahrstart}
  if (is.null(jahrende)==FALSE){j_end<-jahrende}

  probs.male<-read.csv2(paste(c(GetRDirPath(),"/Data/PeriodensterbetafelD_M.csv"),collapse=""),header = TRUE, sep = ";",dec=",", quote="\"")
  names.male<-gsub("X","",names(probs.male))
  names(probs.male)<-names.male
  probs.female<-read.csv2(paste(c(GetRDirPath(),"/Data/PeriodensterbetafelD_W.csv"),collapse=""),header = TRUE, sep = ";",dec=",", quote="\"")
  names.female<-gsub("X","",names(probs.female))
  names(probs.female)<-names.female
  survdata<-new.df(c("group","diagyear","fu1","fu2","fu3","fu4","fu5","obs"))
  for (gr in levels(x[,gruppe]))
  {

    subdata<-drop.levels(subset(x,x[,gruppe]==gr))
    data.year=levels(as.factor(subdata$dy))
    period.result <- period(subdata, 5,
                            probs.male, probs.female, j_start, j_end,
                            method="hakulinen")
    period.result
    rel<-period.result["rel.surv"]
    #fürs debugging
    # print(obs=obs$observations[1])
    #debugging ende
    obs<-period.result["observations"]
    rbind(survdata,data.frame(group=gr,diagyear=gr,fu1=rel$rel.surv[1],fu2=rel$rel.surv[2],
                              fu3=rel$rel.surv[3],fu4=rel$rel.surv[4],fu5=rel$rel.surv[5],obs=obs$observations[1]))->survdata
  }

  survdata[3:7][survdata[3:7]==0]<-NA
  survdata[3:7][survdata[3:7]>100]<-100
  survdata
  obs.dat<-t(survdata[c(8)])
  yhights<-max(obs.dat)*2


  cust.col<-c("gray75","darkblue","darkred","orange","darkgreen","black")


  op<-par(mar=c(5, 4, 4, 12) + 0.1,xpd=TRUE)
  mp<-barplot(obs.dat,yaxt="n",ylim=c(0,yhights),names.arg=survdata[[1]],col=cust.col[1])
  text(mp, obs.dat, labels = obs.dat, pos = 3,cex=0.75)

  axis(4)

  mtext("no. of observations", side=4, line=2, cex=0.8,las=3)
  op<-par(new=TRUE)
  #plot.default(survdata$group, survdata$fu1,ylim=c(0,100),col=2,type="n")
  plot.default(survdata$group, survdata$fu1,ylim=c(0,100),col="white",xaxt="n",xlab=xlab,ylab=ylab,main=titel,cex.main=0.8)
  lines(survdata$group, survdata$fu1,type="l",col=cust.col[2])
  lines(survdata$group, survdata$fu2,type="l",col=cust.col[3])
  lines(survdata$group, survdata$fu3,type="l",col=cust.col[4])
  lines(survdata$group, survdata$fu4,type="l",col=cust.col[5])
  lines(survdata$group, survdata$fu5,type="l",col=cust.col[6])
  tmp.u <- par('usr')
  #leg.pos<-list(x=tmp.u[2]+2, y=tmp.u[4], xjust=0, yjust=0,outer=TRUE)
  leg.pos<-list(x=tmp.u[2]*1.15, y=tmp.u[4], outer=TRUE)
  legend(leg.pos, c(paste("Observations\n( n =",sum(obs.dat),")"),"1 year","2 years","3 years", "4 years", "5 years"), col = cust.col,
         lty = c(1),lwd=c(10,1,1,1,1,1),bty="n",cex=0.8,merge = TRUE)
  par(op)

}

##########################################################################################################################################################################################################
##Funktion für die Darstellung des relativen 1-5-Jahresüberlebn als Kohortenanalyse (Kohorte = Gruppen der Diagnosehjahre)
##########################################################################################################################################################################################################
#' R coalesce Function
#'
#' accepts a list of vectors of identical length and returns one vector with the first non-NA value
#'
#' @param list of vectors of identical length
#'
#' @return one vector with the first non-NA value
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' coalesce(c(NA,NA,53))
#'
#'@export
bbRelSurvCohortPlot<- function(x,ylab=NULL,xlab=NULL,titel=NULL,jahrstart=NULL,jahrende=NULL,jahrintervall=1,method='edererII', lang='de')
{


  j_start<-2003
  j_end<-2012
  if (is.null(jahrstart)==FALSE){j_start<-jahrstart}
  if (is.null(jahrende)==FALSE){j_end<-jahrende}



  probs.male<-read.csv2(paste(c(GetDataPath(),"/PeriodensterbetafelD_M.csv"),collapse=""),header = TRUE, sep = ";",dec=",", quote="\"")
  names.male<-gsub("X","",names(probs.male))
  names(probs.male)<-names.male
  probs.female<-read.csv2(paste(c(GetDataPath(),"/PeriodensterbetafelD_W.csv"),collapse=""),header = TRUE, sep = ";",dec=",", quote="\"")
  names.female<-gsub("X","",names(probs.female))
  names(probs.female)<-names.female
  survdata<-new.df(c("group","diagyear","fu1","fu2","fu3","fu4","fu5","obs"))
  jseq<-seq(j_start,j_end)
  counter=0
  for (j in 1:floor(length(jseq)/(jahrintervall+1)))
  {
    if(counter==0){counter=counter+1}
    lower_j<-jseq[counter]
    upper_j<-jseq[counter+jahrintervall]
    subdata<-drop.levels(subset(x,dy>=lower_j & dy<=upper_j))
    data.year=levels(as.factor(subdata$dy))
    period.result <- period(subdata, 5,
                            probs.male, probs.female, j_start, j_end,
                            method=method)
    period.result
    rel<-period.result["rel.surv"]

    obs<-period.result["observations"]
    rbind(survdata,data.frame(group=if(jahrintervall==0){lower_j}else{paste(lower_j,"-",upper_j)},
                              diagyear=paste(lower_j,"-",upper_j),
                              fu1=if(lower_j>=j_end){NA}else{rel$rel.surv[1]},
                              fu2=if(lower_j+1>=j_end){NA}else{rel$rel.surv[2]},#rel$rel.surv[2],
                              fu3=if(lower_j+2>=j_end){NA}else{rel$rel.surv[3]},
                              fu4=if(lower_j+3>=j_end){NA}else{rel$rel.surv[4]},
                              fu5=if(lower_j+4>=j_end){NA}else{rel$rel.surv[5]},
                              obs=obs$observations[1]))->survdata
    counter=counter+jahrintervall+1
    print(period.result)
  }

  survdata[3:7][survdata[3:7]==0]<-NA
  #survdata[3:7][survdata[3:7]>100]<-100
  print(survdata)
  obs.dat<-t(survdata[c(8)])
  yhights<-max(obs.dat)*2


  cust.col<-c("gray75","darkblue","darkred","orange","darkgreen","black")


  op<-par(mar=c(5, 4, 4, 12) + 0.1,xpd=TRUE)
  mp<-barplot(obs.dat,yaxt="n",ylim=c(0,yhights),names.arg=survdata[[1]],col=cust.col[1],cex.names=0.75)
  text(mp, obs.dat, labels = obs.dat, pos = 3,cex=0.75)

  axis(4)

  mtext("no. of observations", side=4, line=2, cex=0.8,las=3)
  op<-par(new=TRUE)
  #plot.default(survdata$group, survdata$fu1,ylim=c(0,100),col=2,type="n")
  plot.default(survdata$group, survdata$fu1,ylim=c(0,100),col="white",xaxt="n",xlab=xlab,ylab=ylab,main=titel,cex.main=0.8)
  lines(survdata$group, survdata$fu1,type="l",col=cust.col[2])
  lines(survdata$group, survdata$fu2,type="l",col=cust.col[3])
  lines(survdata$group, survdata$fu3,type="l",col=cust.col[4])
  lines(survdata$group, survdata$fu4,type="l",col=cust.col[5])
  lines(survdata$group, survdata$fu5,type="l",col=cust.col[6])
  tmp.u <- par('usr')
  #leg.pos<-list(x=tmp.u[2]+2, y=tmp.u[4], xjust=0, yjust=0,outer=TRUE)
  leg.pos<-list(x=tmp.u[2]*1.15, y=tmp.u[4], outer=TRUE)
  switch(lang,
         en=legend(leg.pos, c(paste("Observations\n( n =",sum(obs.dat),")"),"1 year","2 years","3 years", "4 years", "5 years"), col = cust.col,lty = c(1),lwd=c(10,1,1,1,1,1),bty="n",cex=0.8,merge = TRUE)
         ,de=legend(leg.pos, c(paste("Observationen\n( n =",sum(obs.dat),")"),"1 Jahr","2 Jahre","3 Jahre", "4 Jahre", "5 Jahre"), col = cust.col,lty = c(1),lwd=c(10,1,1,1,1,1),bty="n",cex=0.8,merge = TRUE)
  )
  par(op)

}

