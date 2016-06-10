

##########################################################################################################################################################################################################
## R Function to determine the rough UICC Stage
##########################################################################################################################################################################################################

#' Returns the rough UICC Stage from a detailed UICC Stage
#'
#'
#'
#' @param uiccfull the detailed uicc stage
#'
#' @return Vector of single character string containing the rough UICC Stage
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'

#'
#' @examples
#' roughUICC("IIIa")
#'
#' @export
roughUICC <- function(uiccfull) {

  res<-ifelse(substr(x,1,3)=="III","III",
              ifelse(substr(x,0,2)=="IV","IV",
                     ifelse(substr(x,0,2)=="II","II",
                            ifelse(substr(x,0,1)=="I","I",
                                   ifelse(substr(x,0,1)=="0","0","ND")))))
  return(res)
}

##########################################################################################################################################################################################################
## R Function to determine Age Group
##########################################################################################################################################################################################################

#' Returns an age group for a given age in years
#'
#'
#'
#' @param age the age in years
#' @param br the intervall in years to build groups on
#' @return Vector of single character string containing a string like '20-30'
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' getAgeGroup(53,10)
#'
#'@export
getAgeGroup<-function(age,br){

  res1<-as.character(cut(age,seq(0,120,br)))
  res1<- substr(res1,2,nchar(res1)-1)
  res1<-gsub(",","-",res1)
  return(res1)
}

##########################################################################################################################################################################################################
## R coalesce Function - accepts a list of vectors of identical length and returns one vector with the first non-NA value
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
coalesce <- function(...) {
  Reduce(function(x,y) {
    i<-which(is.null(x))
    x[i]<-y[i]
    x},
    list(...))
}

##########################################################################################################################################################################################################
## Function not in / without
##########################################################################################################################################################################################################
#' R not in / without operator
#'
#' %w/o% is a more intuitive interface as a binary operator, which returns a logical vector indicating if there is no match for its left operand
#'
#' @param x vector: the values not to be
#' @param y vector: the values to be searched against
#' @return
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' c(1:6,7:2) %w/o% c(3,7,12)
#'
#'@export
"%w/o%" <- function(x,y) x[!x %in% y]


##########################################################################################################################################################################################################
##Erstellt einen leeren DataFrame aus den Angaben im Header
##########################################################################################################################################################################################################
#' Function to create an empty named data.frame
#'
#'
#'
#' @param header vector of strings: The names for the columns of the data.frame
#'
#' @return A data.frame
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' ndf<-new.df(c("ID","NAME","AGE"))
#'
#'@export
new.df<- function(header){
  df<-data.frame(matrix(matrix(rep(1,length(header)),1),1))
  colnames(df)<-header
  return(df[NULL,])
}

##########################################################################################################################################################################################################
##Funktion für das erzeugen einer angepassten Farbpalette
##########################################################################################################################################################################################################
#' R Function für generating a custom colorpallette of length x
#'
#' @param x :numeric vector - size of colorpalette
#'
#' @return one character vector of length x containig the custom colors
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' cols<-getBBColors(7)
#'
#'@export
getBBColors<-function(x) {

  cols2<-c("#004992","#ED9B33","#B22229","#7296AF","#8ABD24","#AA9C8F","#BA9BC5","#575756")
  if(x<9){mypalette<-cols2[1:x]
  }else{
    mypalette<-brewer.pal(9,"Set1")
    mypalette<-colorRampPalette(mypalette, space = "Lab")
    mypalette<-mypalette(x)}
  return(mypalette)

}
##########################################################################################################################################################################################################
##Funktion f?r das wirklich l?schen von nicht vorhanden Levels in einem DataFrame
##########################################################################################################################################################################################################
#' Drop unused factor levels from all factors in a data.frame
#'
#' @param dat : data.frame
#'
#' @return dataframe
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' df<-drop.levels(df)
#'
#'@export
drop.levels <- function(dat){
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}

##########################################################################################################################################################################################################
##Funktion f?r Das Subsetting eines Dataframes "x" nach leveln "level", deren vorkommen gr??er/gleich  "mincounts"- mal ist
##########################################################################################################################################################################################################
#' Subset a data.frame by group y where each level has to be larger than size x
#'
#' @param x data.frame: The data.frame to subset
#' @param group character vector: the columnname of the data frame to build groups from
#' @param mincounts numeric vector: the minimal size allowed for a group zo stay in the dataframe
#' @return dataframe
#'
#' @author Stefan Bartels, \email{email@biobits.eu}
#'
#' @examples
#' ID1<-c(2,3,45,343,34,8,77,88,90)
#' ID2<-c("A","A","B","C","C","C","D","D","D")
#' df<-data.frame(ID1,ID2)
#' df<-subset.validlevels(df,"ID2",3)
#'
#'@export
subset.validlevels<-function(x,group,mincounts){

  #list  <- levels(x$group)
  list<-levels(x[,group])
  valid_levels <-c()
  for (t in list)
  {
    tsubdata <-subset(x,x[,group] == t)
    if (length(tsubdata[,group])>=mincounts)
    {
      valid_levels <- c(valid_levels,t)
    }
  }
  subdata<-subset(x,x[,group] %in% valid_levels)
  subdata<-drop.levels(subdata)
  return(subdata)


}

