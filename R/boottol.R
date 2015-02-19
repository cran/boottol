#' Creates 99 percent bootstrap percentile confidence intervals around credit scoring validation statistics
#' @import boot sqldf plyr tcltk
#' @param score - The score groupings.
#' @param target - The binary target variable.
#' @param bootsamp - How many bootstrap samples to be computed. When bootsamp is too low, a warning will be produced.
#' @param lth - If the score is rank ordered least likely to most likely then set equal to True, else False.
#' @param grp - An integer value of how the scores are grouped.
#' @examples
#' data("data")
#' \dontrun{tol.level=boottol(score=data$Score,target=data$Target,bootsamp=2000,lth="True",grp=10)}
#' @export
boottol <-function(score,target,bootsamp,lth, grp)
  {
    
    a="scre"
    b="trget"
    data1=setNames(data.frame(score,target), c(a,b))
    
    dataset <- data.frame(Statistic=character(),
                          Cutoff=integer(),
                          Lowerbound99=double(),
                          Upperbound99=double(),
                          ObservedStat=double(),
                          stringsAsFactors=FALSE) 
    
    for(j in seq(from=0, to=(max(score)-grp), by=grp))
    {
      ###############################################################################
      ################ SEGMENT DATA SET       #######################################
      ###############################################################################
      
      data2=data1[which(data1$scre>=j),]
      
      ###############################################################################
      ################ TOLERANCE FOR KS       #######################################
      ###############################################################################
      
      #KS function
      KSF <- function(d, i)
      {
        d2 <- d[i,]
        data=sqldf('
                   Select scre, sum(trget) as Bad, sum(1-trget) as Good  
                   from d2
                   group by scre
                   order by scre Desc')
        data$CumBad=cumsum(data$Bad)
        data$CumGood=cumsum(data$Good)
        data$CumBadPer=data$CumBad/max(data$CumBad)
        data$CumGoodPer=data$CumGood/max(data$CumGood)
        return(max(abs(data$CumGoodPer-data$CumBadPer)))
      }
      
      KSboot = boot(data2,KSF,bootsamp)
      KSCI=boot.ci(KSboot, conf = 0.99, type =c("norm","basic", "perc"))
      LCI_KS=KSCI$percent[4]
      UCI_KS=KSCI$percent[5]
      obsKS=KSboot$t0
      info=c('KS',j,LCI_KS,UCI_KS,obsKS)
      
      
      ###############################################################################
      ################ TOLERANCE FOR AUROC    #######################################
      ###############################################################################
      
      shift<-function(x,shift_by){
        stopifnot(is.numeric(shift_by))
        stopifnot(is.numeric(x))
        
        if (length(shift_by)>1)
          return(sapply(shift_by,shift, x=x))
        
        out<-NULL
        abs_shift_by=abs(shift_by)
        if (shift_by > 0 )
          out<-c(tail(x,-abs_shift_by),rep(NA,abs_shift_by))
        else if (shift_by < 0 )
          out<-c(rep(NA,abs_shift_by), head(x,-abs_shift_by))
        else
          out<-x
        out
      }
      
      
      
      if (substr(lth,1,1) == "T"){
        AUCF <- function(d, i)
        {
          d2 <- d[i,]
          data=sqldf('
                     Select scre, sum(trget) as Bad, sum(1-trget) as Good  
                     from d2
                     group by scre
                     order by scre')
          data$CumBad=cumsum(data$Bad)
          data$CumGood=cumsum(data$Good)
          data$CumBadPer=data$CumBad/max(data$CumBad)
          data$CumGoodPer=data$CumGood/max(data$CumGood)
          data$CumBadPerLag=shift(data$CumBadPer,-1)
          data$CumGoodPerLag=shift(data$CumGoodPer,-1)
          data[is.na(data)]=0
          data$auroc=.5*(data$CumBadPer+data$CumBadPerLag)*(data$CumGoodPer-data$CumGoodPerLag)
          return(sum(data$auroc))
        }
        
        AUCboot = boot(data2,AUCF,bootsamp)
        AUCCI=boot.ci(AUCboot, conf = 0.99, type =c("norm","basic", "perc"))
        LCI_AUROC=AUCCI$percent[4]
        UCI_AUROC=AUCCI$percent[5]
        obsAUC=AUCboot$t0
        info2=c('AUROC',j,LCI_AUROC,UCI_AUROC,obsAUC)
      } else {
        AUCF <- function(d, i)
        {
          d2 <- d[i,]
          data=sqldf('
                     Select scre, sum(trget) as Bad, sum(1-trget) as Good  
                     from d2
                     group by scre
                     order by scre desc')
          data$CumBad=cumsum(data$Bad)
          data$CumGood=cumsum(data$Good)
          data$CumBadPer=data$CumBad/max(data$CumBad)
          data$CumGoodPer=data$CumGood/max(data$CumGood)
          data$CumBadPerLag=shift(data$CumBadPer,-1)
          data$CumGoodPerLag=shift(data$CumGoodPer,-1)
          data[is.na(data)]=0
          data$auroc=.5*(data$CumBadPer+data$CumBadPerLag)*(data$CumGoodPer-data$CumGoodPerLag)
          return(sum(data$auroc))
        }
        
        AUCboot = boot(data2,AUCF,bootsamp)
        AUCCI=boot.ci(AUCboot, conf = 0.99, type =c("norm","basic", "perc"))
        LCI_AUROC=AUCCI$percent[4]
        UCI_AUROC=AUCCI$percent[5]
        obsAUC=AUCboot$t0
        info2=c('AUROC',j,LCI_AUROC,UCI_AUROC,obsAUC)
      }
      
      ###############################################################################
      ################ TOLERANCE FOR GINI     #######################################
      ###############################################################################
      
      
      LCI_GINI= (2*LCI_AUROC)-1
      UCI_GINI= (2*UCI_AUROC)-1
      obsGI= (2*obsAUC)-1
      info3=c('GINI',j,LCI_GINI,UCI_GINI,obsGI)
      
      
      ###############################################################################
      ################ PIECE DATASET TOGETHER #######################################
      ###############################################################################
      
      
      for(k in 1:ncol(dataset))
      {
        dataset[((j/10)+j+1),k]=info[k]
      }
      
      for(f in 1:ncol(dataset))
      {
        dataset[((j/10)+j+2),f]=info2[f]
      }
      
      for(x in 1:ncol(dataset))
      {
        dataset[((j/10)+j+3),x]=info3[x]
      }
      
    }
    na.omit(dataset)
  }
