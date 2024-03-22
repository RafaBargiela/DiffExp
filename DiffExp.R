#########################################################################################
###     DIFFERENTIAL EXPRESSION ANALYSIS    last update: 20th of March, 2024  #######
### Author: Rafael Bargiela, PhD
### Summary: Set of subroutines for the analysis of diferential expression of genes or transcriptomics data, or Taxa
### 		based on the same pipeline used for proteomics analysis 
#########################################################################################

DifExp<-function(M,groups.vector,reference.group=1,tune.sigma=1,q=0.01){
		# M: Abundance Matrix.
      	# groups.vector : Numerical or Character vector, defining the group assigned to each column
		# tune.sigma and q: these are values for the imputation by Minimum Probability Method of Na values after log2 calculation. Tune.sigma is constant and q are the probabilities to asses quantiles of the distribution
		# reference.group : Numeric. Position in the matrix or groups.vector where is located the reference group, the one which will be compared to the rest of them.
		DF<-{}
		# Data normalization over the median
		message("- Normalizing data over the median and getting log2")
		DF$Median.Norm<-sapply(1:ncol(M),function(x){
        log2<-log(as.matrix(M[,x]),base=2)
        log2[is.infinite(log2)]<-NA
        gMedian<-median(log2,na.rm=TRUE)
        log2-gMedian
      	})
      	rownames(DF$Median.Norm)<-rownames(M)
      	colnames(DF$Median.Norm)<-colnames(M)

      	# Data imputation (MinProb method)
      	message("- Imputation of missing values by Minimum probability method")
      	DF$Imputed<-DF$Median.Norm
      	MATGsd<-apply(DF$Median.Norm,1,sd,na.rm=TRUE) # Getting sd of each MATG, excluding missing values
      	MATGsdMed<-median(MATGsd,na.rm=TRUE)*tune.sigma # SD to get the imputed values, based on median proteins SD and tuned by constant set on tune.sigma (if 1, just ProtSD median)
      	for(c in 1:ncol(DF$Median.Norm)){
      		sample<-DF$Median.Norm[,c]
    		sample[is.infinite(sample)]<-NA # if inf haven't been switch to NA yet
    		min.quantile<-quantile(sample,probs = q, na.rm=TRUE) # Getting minimum value for the q quantile
    		dataset.to.impute.miss<-rnorm(length(sample),mean = min.quantile,sd=MATGsdMed) # Getting random distribution of values based on SD for imputation and minimum quantile value as mean
    		sample.NAs<-which(is.na(sample))  # Changing Missing values for values on the created dataset 
    		sample[sample.NAs]<-dataset.to.impute.miss[sample.NAs]
    		DF$Imputed[,c]<-sample
      	}

      	# Log2 Fold Change analysis
      	message("- Calculating Fold Change and Statistics")
      	message("NOTE: Equal variance among groups is assumed")
      	G<-unique(groups.vector)
      	Rname<-G[reference.group]
      	R<-DF$Imputed[,grep(G[reference.group],groups.vector,value=FALSE)] # Reference group Data
      	Rav<-sapply(1:nrow(R),function(x){mean(R[x,],na.rm = TRUE)})
      	names(Rav)<-rownames(R)
      	n1<-sapply(1:nrow(R),function(x){sum(!is.na(R[x,]))}) # nr of cases for each taxon on Reference group
      	# Summ of Squares for reference group
      	SSER<-sapply(1:nrow(R),function(x){sum((as.numeric(na.omit(R[x,]))-Rav[x])^2)})
      	C<-G[-reference.group] # Comparison groups
      	for(g in 1:length(C)){
      		Cgname<-C[g]
      		Cg<-DF$Imputed[,grep(C[g],groups.vector,value=FALSE)]
      		Cgav<-sapply(1:nrow(Cg),function(x){mean(Cg[x,],na.rm = TRUE)})
      		names(Cgav)<-rownames(Cg)
      		# Fold change
      		FC<-Rav-Cgav
      		  # T-test for each protein among two different groups ##
		  # NOTE: assuming equal variances and different groups size
		  n2<-sapply(1:nrow(Cg),function(x){sum(!is.na(Cg[x,]))})  # nr of cases for each taxon on each comparison group
		  # Summ of Squares for comparison group
		  SSECg<-sapply(1:nrow(Cg),function(x){sum(as.numeric((na.omit(Cg[x,]))-Cgav[x])^2)})
		  SSE<-SSER+SSECg # Total sum of squares
		  Sp2<-SSE/(n1+n2-2) # Pooled variance of the two groups (n1+n2-2 degrees of freedom)
		  ESE<-sqrt((Sp2/n1)+(Sp2/n2))
		  t.test<-FC/ESE
		  p.values<-sapply(1:length(t.test),function(x){
		    if(t.test[x]<0){
		      2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=TRUE)                
		    }else{
		      2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=FALSE)
		    }
		  })
		  p.values.adj<-p.adjust(p.values,method="fdr")

		  # Summary matrix
		  name<-paste(Rname," vs ",Cgname,sep="")
		  DF$FCstats[[name]]<-matrix(c(Rav,Cgav,FC,p.values,p.values.adj),nc=5,dimnames=list(rownames(DF$Imputed),c(paste(Rname,"log2 mean"),paste(Cgname,"log2 mean"),"Fold Change","t.test p-values","adjusted p-values")))
      	}

      	return(DF)
	}

volcano<-function(FC,logPv,FCthreshold=1,PVthreshold=2,vol.xlim=c(-5,5),vol.ylim=c(0,5),xlim.ticks=seq(vol.xlim[1],vol.xlim[2],2.5),vol.xlim.labs=TRUE,vol.ylim.labs=TRUE,volpch=21,points.labs=rownames(FC),highlight=NULL,points.cex=1.5,points.labs.cex=1,volcol=hcl.colors(sum(length(FCthreshold),length(PVthreshold))+1,"plasma"),legend.cols=TRUE,x.legcol=vol.xlim[1]+(abs(vol.xlim[1])*0.01),y.legcol=vol.ylim[2]-(abs(vol.ylim[2])*0.01),x.legcolabsep=((abs(vol.xlim[1])+abs(vol.xlim[2]))*0.02),y.legcolabbottom=(vol.ylim[2]-vol.ylim[1])/1.5,legcol.cex=1,legcol.points.cex=points.cex,volXlabCorr=rep(0.05,nrow(FC)),volYlabCorr=rep(0,nrow(FC)),legend.dots=FALSE,LegdotsLabs="",xlab=list("log\u2082(Fold change)",cex=1.2,font=2),ylab=list("-log\u2081\u2080(p-value)",cex=1.2,font=2),overRepLab="Overrepresented",x.overRepLab=vol.xlim[1],y.overRepLab=vol.ylim[1],x.underRepLab=vol.xlim[2],y.underRepLab=vol.ylim[1],underRepLab="Underrepresented",xlegcol=vol.xlim[1],ylegcol=vol.ylim[2],xlegcolabsep){

  # FC : log2 FC matrix.
  # logPV : -log 10 of the adjusted p-values.
  # legend.cols: Logical. If TRUE, it shows the legend corresponding to each color. Usually related with FC and p-value threshold
  # x.legcol and y.legcol: Initial x and y point where to draw the legend.cols
  # x.legcollabsep: Separation between dots and text in legend.cols
  # y.legcolabbottom: y bottom point for the legend, been y.legcol the top.
  # legcol.cex: size for legend.cols labels
  # legcol.points.cex: size for legend.cols dots
  # highlight: To mark the presence of a specific type of element, independently of been differentially expressed or not.

  # Chechings -----------------------------------------------------------------
    require(TeachingDemos)
    if(length(volcol)!=sum(length(FCthreshold),length(PVthreshold))+1){
      warning("- Number of colors for dots different than sum of thresholds divisions\n",paste("Colors:",length(volcol),"\n"),paste("Thresholds:",sum(length(FCthreshold),length(PVthreshold))+1,"\n"))
    }
    if(length(volXlabCorr)!=nrow(FC)){
      warning(paste("Nr of  X-axis corrections",length(volXlabCorr),"for labels must be equal to row names of FC",nrow(FC)))
    }
    if(length(volYlabCorr)!=nrow(FC)){
      warning(paste("Nr of  X-axis corrections",length(volYlabCorr),"for labels must be equal to row names of FC",nrow(FC)))
    }
    if(legend.dots==TRUE & length(volpch)==1){
      warning(paste("Only 1 dot type specified (volpch)\n","Samples to compare against the reference sample (first column): ",ncol(FC),"\n",sep=""))
    }       
  # End checkings -------------------------------------------------------------
    MostExpDF<-data.frame()
    MostExp<-vector("character")
    plot(NA,NA,xlim=vol.xlim,ylim=vol.ylim,ylab=ylab,xlab=xlab,axes=FALSE)
    axis(1,xlim.ticks,labels=vol.xlim.labs,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
    axis(2,seq(vol.ylim[1],vol.ylim[2],1),labels=vol.ylim.labs,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
    arrows(c(-FCthreshold,FCthreshold),vol.ylim[1],c(-FCthreshold,FCthreshold),vol.ylim[2],code=0,lty=3,lwd=2,col="grey50")
    arrows(vol.xlim[1],PVthreshold,vol.xlim[2],PVthreshold,code=0,lty=3,lwd=2,col="grey50")
    rect(0,0,vol.xlim[1],vol.ylim[2],border=NA,col=rgb(t(col2rgb("thistle2")),alpha=100,maxColorValue=255),xpd=TRUE)
    text(x.overRepLab,y.overRepLab,labels=overRepLab,adj=c(0,0),cex=1.2,font=2,xpd=TRUE,col="mediumorchid3")
    rect(0,0,vol.xlim[2],vol.ylim[2],border=NA,col=rgb(t(col2rgb("grey90")),alpha=100,maxColorValue=255),xpd=TRUE)
    text(x.underRepLab,y.underRepLab,labels=underRepLab,adj=c(1,0),cex=1.2,font=2,xpd=TRUE,col="grey50")            
    for(r in 1:nrow(FC)){
      if(length(highlight)>0){
          if(sum(grepl(rownames(FC)[r],highlight))>0){
            message(paste(rownames(FC)[r],"will be highlighted"))
            next
          }else{

          }
      }
      #
      for(c in 1:ncol(FC)){
        n<-1
        topFC<-sum(abs(FC[r,c])>=FCthreshold)
        topPV<-sum(abs(logPv[r,c])>=PVthreshold)
        if(topFC==0 & topPV==0){
          n<-1
        }else{
          if((topFC>0 & topPV==0) | (topFC==0 & topPV>0)){
            n<-n+1
          }else{
            n<-n+topFC+topPV
          }
        } 

        points(FC[r,c],logPv[r,c],pch=volpch[c],col=volcol[n],bg=rgb(t(col2rgb(volcol[n])),alpha=100,maxColorValue=255),cex=points.cex,lwd=2)
        if(n==3){
          # print(rownames(FC)[r])
          if(is.list(volXlabCorr)==TRUE){
              xc<-volXlabCorr[[c]][r]
            }else{
              xc<-volXlabCorr[r]
            }
          if(is.list(volYlabCorr)==TRUE){
              yc<-volYlabCorr[[c]][r]
            }else{
              yc<-volYlabCorr[r]
            }
            if(FC[r,c]<0){
              adj<-c(1,0.5)
              xc<-xc*(-1)
            }else{
              adj<-c(0,0.5)
            }
        
          # print(c(c,rownames(FC)[r],yc))
            text(FC[r,c]+xc,logPv[r,c]+yc,labels=points.labs[r],adj=adj,cex=points.labs.cex,font=4,xpd=TRUE)
            MostExp<-append(MostExp,rownames(FC)[r],length(MostExp))      
        }
      }
    }

    # Highlighting specific proteins
    if(length(highlight)>0){
      for(e in 1:length(highlight)){
        for(c in 1:ncol(FC)){
          n<-1
          topFC<-sum(abs(FC[highlight[e],c])>=FCthreshold)
          topPV<-sum(abs(logPv[highlight[e],c])>=PVthreshold)
          if(topFC==0 & topPV==0){
            n<-1
          }else{
            if((topFC>0 & topPV==0) | (topFC==0 & topPV>0)){
              n<-n+1
            }else{
              n<-n+topFC+topPV
            }
          } 
          points(FC[highlight[e],c],logPv[highlight[e],c],pch=volpch[c],col="white",bg="transparent",cex=points.cex+0.2,lwd=2)
          points(FC[highlight[e],c],logPv[highlight[e],c],pch=volpch[c],col=volcol[n],bg=rgb(t(col2rgb(volcol[n])),alpha=100,maxColorValue=255),cex=points.cex,lwd=2)         
          p<-grep(highlight[e],rownames(FC))
          if(is.list(volXlabCorr)==TRUE){
              xc<-volXlabCorr[[c]][p]
            }else{
              xc<-volXlabCorr[p]
            }
              if(is.list(volYlabCorr)==TRUE){
              yc<-volYlabCorr[[c]][p]
            }else{
              yc<-volYlabCorr[p]
            } 
          if(FC[highlight[e],c]<0){
            adj<-c(1,0.5)
            xc<-xc*(-1)
          }else{
            adj<-c(0,0.5)
          }
          if(length(points.labs)==length(highlight)){
            shadowtext(FC[highlight[e],c]+xc,logPv[highlight[e],c]+yc,labels=points.labs[highlight[e]],adj=adj,cex=1,font=2,xpd=TRUE,col="black",bg="white",r=0.05) 

          }else{          
            shadowtext(FC[highlight[e],c]+xc,logPv[highlight[e],c]+yc,labels=points.labs[p],adj=adj,cex=1,font=2,xpd=TRUE,col="black",bg="white",r=0.05)              
          }
        }
      }
    }
    box(lwd=2)
    # Sample type Legend, when using more than one dot shape
      if(legend.dots==TRUE){
        xi<-vol.xlim[1]
        y<-vol.ylim[2]
        yw<-(y-(y-(y*0.17)))/length(LegdotsLabs)
        for(p in 1:length(volpch)){
          points(xi,y,pch=volpch[p],cex=2,lwd=2,bg="white",col="black",xpd=TRUE)
          text(xi+0.3,y,labels=LegdotsLabs[p],font=2,cex=1.2,xpd=TRUE,adj=c(0,0.5))
          y<-y-yw
        } 
      }
    # Thresholds Legend related to colors
      if(legend.cols==TRUE){
        legcolsLabs<-c(paste("log\u2082(FC)<",FCthreshold[1],", p>",10^(-PVthreshold[1]),sep=""),paste("log\u2082(FC)>",FCthreshold[1],"OR p<",10^(-PVthreshold[1]),sep=""),paste("log\u2082(FC)>",FCthreshold[1],"& p<",10^(-PVthreshold[1]),sep=""))
        FCleg<-vector("character",length(FCthreshold))
        PVleg<-vector("character",length(PVthreshold))
        FCleg[1]<-as.character(FCthreshold[1])
        PVleg[1]<-as.character(10^(-PVthreshold[1]))
        if(length(FCthreshold)>1){
          for(th in 2:length(FCthreshold)){
            FCleg[th]<-paste("log\u2082(FC)>",FCthreshold[th],sep="")
          }
        }
        if(length(PVthreshold)>1){
          for(th in 2:length(PVthreshold)){
            PVleg[th]<-paste("p<",10^(-PVthreshold[th]),sep="")
          }
        }
        FCPVleg<-paste(FCleg,"&",PVleg)
        if(length(FCPVleg)>1){
          legcolsLabs<-append(legcolsLabs,FCPVleg[2:length(FCPVleg)],length(legcolsLabs))
        }
        y.legcolabsep<-(y.legcol-y.legcolabbottom)/length(volcol)
        for(col in 1:length(volcol)){
          points(x.legcol,y.legcol,pch=volpch[1],cex=legcol.points.cex,col=volcol[col],bg=rgb(t(col2rgb(volcol[col])),alpha=100,maxColorValue=255),lwd=2)
          if(x.legcol<0){
            adj<-c(0,0.5)
            sep<-x.legcolabsep
          }else{
            adj<-c(1,0.5)
            sep<-x.legcolabsep*(-1)
          }
          text(x.legcol+sep,y.legcol,labels=legcolsLabs[col],adj=adj,font=2,cex=legcol.cex)
          y.legcol<-y.legcol-y.legcolabsep
        }
      }
    MostExM<-FC[unique(MostExp),]
    return(MostExM)
}

QQPlot<-function(Alog,FCsd,SD,sdrange=c(1*SD,2*SD,3*SD),q.xlim=c(-5,5),q.ylim=c(-5,5),x.tick=seq(q.xlim[1],q.xlim[2],2.5),y.tick=seq(q.ylim[1],q.ylim[2],2.5),QXlabCorr=rep(0.05,nrow(FCsd)),QYlabCorr=rep(0,nrow(FCsd)),highlight=NULL,points.labs=rownames(Alog),points.cex=1.5,Qpch=21,Qcol=hcl.colors(length(sdrange)+1,"plasma"),legend=FALSE,leglabs=NULL,y.legsep=4,LegdotsLabs="",labcex=1,over.x.corr=0,over.y.corr=0,under.x.corr=0,under.y.corr=0,MostExpLabs=FALSE,MostExpThreshold=2,yl.corr=0,xl.corr=0,overRepLab="Overrepresented",underRepLab="Underrepresented",xlab="",ylab=""){
  # Alog: Abundances matrix, usually in logarithmns. NOTE: Ensure first column corresponds to control/reference
  # FCsd:  Matrix with the log2(fold-change)/SD per element (protein,gene,transcript)
  # SD:  Standard deviation of the log2(FC) for all elements
  # MostExpLabs : Show labs for elements most expressed above a threshold
  # MostExpThreshold : Threshold to show Most Expressed elements (genes, transcripts, proteins)
  # y.legsep: Percentage of separataion on y-axis for Color legend

  # CHECKINGS ----------------------------------------------------------------------------
  if(ncol(Alog)!=(ncol(FCsd)+1)){
    stop(paste("ERROR: number of columns in Abundance Matrix must be equal to number of columns in Fold Change matrix +1: ",ncol(Alog),"!=",ncol(FCsd),sep=""))
   }
  # END CHECKINGS ------------------------------------------------------------------------

  require(TeachingDemos)
  plot(NA,NA,xlim=q.xlim,ylim=q.ylim,axes=FALSE,xlab=xlab,ylab="")
    text(q.xlim[1]-1.5,0,labels=ylab,cex=1.2,font=2,srt=90,xpd=TRUE)
    axis(1,x.tick,labels=TRUE,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
    axis(2,y.tick,labels=TRUE,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
    polygon(c(q.xlim[1],q.xlim[2],q.xlim[1]),c(q.ylim[1],q.ylim[2],q.ylim[2]),border=NA,col=rgb(t(col2rgb("thistle2")),alpha=100,maxColorValue=255),xpd=TRUE)
  polygon(c(q.xlim[1],q.xlim[2],q.xlim[2]),c(q.ylim[1],q.ylim[1],q.ylim[2]),border=NA,col=rgb(t(col2rgb("grey90")),alpha=100,maxColorValue=255),xpd=TRUE)
  text(q.xlim[2]-0.75+over.x.corr,q.ylim[2]+over.y.corr,labels=overRepLab,adj=c(1,1),cex=1.2,font=2,xpd=TRUE,col="mediumorchid3")
  text(q.xlim[2]+under.x.corr,q.ylim[1]-0.1+under.y.corr,labels=underRepLab,adj=c(1,0),cex=1.2,font=2,xpd=TRUE,col="grey50")                              
    A<-Alog[,1]
    if(ncol(Alog)>2){
      B<-sapply(2:ncol(Alog),function(x){
        v<-vector("numeric")
        v<-append(v,Alog[,x],length(v))
        v
      })
    }else{
      B<-Alog[,2]
    }

    # Calculating quantiles ----------------------------------------------------
    # Q1<-Alog2[round(seq(0.1,0.9,0.1)*(length(Alog2[,1])-1)),1]
    Q1<-quantile(A,seq(0.1,0.9,0.1))
    # Q2<-B[round(seq(0.1,0.9,0.1)*(length(B)-1))]
    Q2<-quantile(B,seq(0.1,0.9,0.1))
    arrows(q.xlim[1]-0.5,Q2,q.xlim[2]+0.5,Q2,lwd=1,lty=3,code=0,col="grey50")
    arrows(Q1,q.ylim[1]-0.5,Q1,q.ylim[2]+0.5,lwd=1,lty=3,code=0,col="grey50")
    text(quantile(A,.10),q.ylim[1]-0.25,labels="10%",font=2,cex=0.9,adj=c(1,0.5))
    text(quantile(A,.9),q.ylim[1]-0.25,labels="90%",font=2,cex=0.9,adj=c(0,0.5))
    text(q.xlim[2],quantile(B,.10),labels="10%",font=2,cex=0.9,adj=c(0.5,1))
    text(q.xlim[2],quantile(B,.9),labels="90%",font=2,cex=0.9,adj=c(0.5,0))
    # Printing dots-------- ----------------------------------------------------
    MostExpDF<-data.frame()
    MostExp<-vector("character")
    for(r in 1:nrow(Alog)){
      for(c in 2:ncol(Alog)){ # Regards that it starts from the second column
        n<-sum(abs(FCsd[r,c-1])>=sdrange)
        col<-n+1
        points(Alog[r,1],Alog[r,c],pch=Qpch[c-1],cex=points.cex,lwd=2,col=Qcol[col],bg=rgb(t(col2rgb(Qcol[col])),alpha=100,maxColorValue=255))
          if(is.list(QXlabCorr)==TRUE){
            xc<-QXlabCorr[[c-1]][r]
          }else{
            xc<-QXlabCorr[r]
          }
            if(is.list(QYlabCorr)==TRUE){
            yc<-QYlabCorr[[c-1]][r]
          }else{
            yc<-QYlabCorr[r]
          }
        if(MostExpLabs==TRUE){
            if(abs(FCsd[r,c-1])>=MostExpThreshold){   # Printing labels for most expressed proteins
            adj<-c(0,0.5)
            if(Alog[r,1]<Alog[r,c]){      # ... but most expressed proteins on B are adjusted to the left
              adj<-c(1,0.5)
              xc<-xc*(-1)
            }
            text(Alog[r,1]+xc,Alog[r,c]+yc,labels=points.labs[r],font=4,cex=labcex,adj=adj,xpd=TRUE)
            MostExp<-append(MostExp,rownames(Alog)[r],length(MostExp))
          } 
        }         
      }
    }
    # Highlighting specific proteins
    if(length(highlight)>0){
      for(e in 1:length(highlight)){
        for(c in 2:ncol(Alog)){
          if(abs(FCsd[highlight[e],c-1])<SD){
              col<-1
            }else{
              col<-(sum(sdrange<=abs(FCsd[highlight[e],c-1])))+1
            }
          points(Alog[highlight[e],1],Alog[highlight[e],c],pch=Qpch[c-1],col="white",bg="transparent",cex=points.cex+0.2,lwd=2)         
          points(Alog2[highlight[e],1],Alog2[highlight[e],c],pch=Qpch[c-1],col=Qcol[col],bg=rgb(t(col2rgb(Qcol[col])),alpha=100,maxColorValue=255),cex=points.cex,lwd=2)

          p<-grep(highlight[e],rownames(Alog))
          if(is.list(QXlabCorr)==TRUE){
              xc<-QXlabCorr[[c-1]][p]
            }else{
              xc<-QXlabCorr[p]
            }
              if(is.list(QYlabCorr)==TRUE){
              yc<-QYlabCorr[[c-1]][p]
            }else{
              yc<-QYlabCorr[p]
            } 
          if(FCsd[highlight[e],c-1]<0){
            adj<-c(1,0.5)
            xc<-xc*(-1)
          }else{
            adj<-c(0,0.5)
          }
          if(length(points.labs)==length(highlight)){
            # print(highlight[e])
            shadowtext(Alog2[highlight[e],1]+xc,Alog2[highlight[e],c]+yc,labels=points.labs[highlight[e]],adj=adj,cex=1,font=4,xpd=TRUE,col="black",bg="white",r=0.05)      
          }else{
            # print(points.labs[highlight[e]])
            shadowtext(Alog[p,1]+xc,Alog[p,c]+yc,labels=points.labs[highlight[e]],adj=adj,cex=1,font=4,xpd=TRUE,col="black",bg="white",r=0.05)              
          }
        }
      }
    }


    # Legend
      if(legend==TRUE){
        yl<-q.ylim[2]+yl.corr
        yw<-(abs(q.ylim[1])+abs(q.ylim[2]))*y.legsep/100
        xw<-((abs(q.xlim[1])-0.2)+abs(q.xlim[2]))*2/100
        if(is.null(leglabs)==TRUE){
          leglabs<-vector("character",length(sdrange)+1)
          for(r in 1:(length(sdrange)+1)){
            if(r==1){
              leglabs[r]<-"FC < \u03c3"
            }else{
              if(r==2){
                leglabs[r]<-"FC > \u03c3"               
              }else{
                leglabs[r]<-paste("FC >",r-1,"\u03c3")
              }
            }
          }
          # print(leglabs)
        }
        for(c in 1:length(Qcol)){
          points(q.xlim[1]-0.2,yl,pch=19,col=Qcol[c],cex=1.2)
          text(q.xlim[1]+xw-0.2,yl,labels=leglabs[c],font=2,cex=1.2,xpd=TRUE,adj=c(0,0.5))
          yl<-yl-yw
        }
        if(length(Qpch)>1){
        xi<-q.xlim[1]+xl.corr
        y<-yl-0.1
        yw<-(y-2)/length(LegdotsLabs)
        for(p in 1:length(Qpch)){
          points(xi,y,pch=Qpch[p],cex=2,lwd=2,bg="white",col="black",xpd=TRUE)
          text(xi+0.3,y,labels=LegdotsLabs[p],font=2,cex=1.2,xpd=TRUE,adj=c(0,0.5))
          y<-y-yw
        }   
        }     
      }
    box(lwd=2)
    MostExM<-FCsd[unique(MostExp),]
    return(MostExM)
}

