
forOptions<-function(data,Sm,Km,R,Psh,selSets,PbSsh,Bopt,colsNopt,rowsNopt){
library(numbers)
library(tidyr)
library(dplyr)
library(stringr)

data$Group=as.factor(data$Group)
data$Role=as.factor(data$Role)
data$GID=as.factor(data$GID)
data=data[order(data$Group),]

dataCc=data[which(data$Role=="Check"),]
dataCc$SetNum=rep("All",dim(dataCc)[1])
dataNCc=data[-which(data$Role=="Check"),]
dataNCc$Group=as.factor(as.character(dataNCc$Group))
dataNCc$Role=as.factor(as.character(dataNCc$Role))
dataNCc$GID=as.factor(as.character(dataNCc$GID))

ng=nlevels(dataNCc$Group)
minGen=unlist(lapply(split(dataNCc$GID,dataNCc$Group),length))
#minGen=min(unlist(lapply(split(dataNCc$GID,dataNCc$Group),length)))
#H=minGen*ng
H=length(unique(dataNCc$GID))
Hfile=H
Ncheck=length(unique(dataCc$GID))
#newG=unlist(lapply(split(dataNCc$GID,dataNCc$Group),function(x){x[1:minGen]}))
#dataNCc=dataNCc[which(dataNCc$GID%in%newG==TRUE),]
dataNCc$IDnum=1:dim(dataNCc)[1]
dataNCc$SetNum=rep("Set",dim(dataNCc)[1])
rm(data)  

######################################################
#Km=Sm-1
#By managment
#Q:Number of sets
#P:Percentage of saving, lower than 30%
#Z1:site 1 size
#Zq:site 2-q size

if(is.null(Sm) && is.null(Km)){
  Sm=list()
  Km=list()
  cont=1
  for (i in 10:3){
    for ( j in (i-1):2){
      Sm[[cont]]=i
      Km[[cont]]=j
      cont=cont+1
    }
  }
  Sm=unlist(Sm);Km=unlist(Km)
}

if(is.null(Sm) && !is.null(Km)){
  Sm=list()
  cont=1
  for (i in 10:3){
      Sm[[cont]]=i
      cont=cont+1
    }
  Sm=unlist(Sm);Km=rep(Km,length(Sm))
}

if(!is.null(Sm) && is.null(Km)){
  Km=list()
  cont=1
    for ( j in (Sm-1):2){
      Km[[cont]]=j
      cont=cont+1
    }
  Sm=rep(Sm,length(Km));Km=unlist(Km)
}
SmKm=paste(Sm,Km,sep="_")

#Values for H and R
#if(!is.null(PbS) && is.null(R)){
#  R=2:3
#  H=round(PbS/R)
#  PbS=H*R
#}
#if(!is.null(PbS) && !is.null(R)){
#  H=round(PbS/R)
#  PbS=H*R
#}
if(is.null(R)){
  R=2:3
  PbS=H*R
}else{
	PbS=H*R
}

#Values for Percent of genotypes in SET-1
#if(is.null(P)){
  P=seq(1,H-1,by=1)
#}

  RHPbS=paste(R,H,PbS,sep="_")
  allini=crossing(RHPbS,SmKm,P)
  allini=as.data.frame(cbind(str_split_fixed(as.character(allini$RHPbS), "_", 3),str_split_fixed(as.character(allini$SmKm), "_", 2),allini$P))
  names(allini)=c("R","H","PbS","Sm","Km","P")
  allini$R=as.numeric(as.character(allini$R))
  allini$H=as.numeric(as.character(allini$H))
  allini$PbS=as.numeric(as.character(allini$PbS))
  allini$Sm=as.numeric(as.character(allini$Sm))
  allini$Km=as.numeric(as.character(allini$Km))
  allini$P=as.numeric(as.character(allini$P))
  
  Q=choose(allini$Sm,allini$Km)+1
  
  #Z1=round(allini$H*allini$P/100)
  Z1=P
  Z2_q=(allini$H-Z1)/(Q-1)
  
#C:Number of checks
#N:Number of genotypes
tot=Z1+(Q-1)*Z2_q

allini=cbind(allini,Q,Z1,Z2_q,tot)

allini=allini[which((allini$H-allini$Z1)%%(allini$Q-1)==0),]

allinilist=split(allini,1:dim(allini)[1])

#datatmp=allinilist[[1]]
adjparms<-function(datatmp){
  Padj=datatmp$P
  Z1adj=datatmp$Z1
  Z2_qadj=datatmp$Z2_q
  Q_adj=datatmp$Q
  
  #if(datatmp$tot>datatmp$H){
  #  a=seq(datatmp$P+1,datatmp$P+5)
  #  b=seq(datatmp$P-1,datatmp$P-5)
  #  Ptest=1:10
  #  Ptest[seq(1,10,by=2)]=a
  #  Ptest[seq(2,10,by=2)]=b
  #  deltest=which(Ptest<=0 & Ptest>30)
  #  if(length(deltest)!=0){Ptest=Ptest[-deltest]}
  #  Z1tmp=round((Ptest*datatmp$H)/100)
  #	Z2_qtmp=round((datatmp$H-Z1tmp)/(datatmp$Q-1))
  #  totmp=Z1tmp+(datatmp$Q-1)*Z2_qtmp
  #  if(length(which(totmp<=datatmp$H))!=0){
  #    seltmp=min(which(totmp<=datatmp$H))
  #    Padj=Ptest[seltmp]
  #    Z1adj=Z1tmp[seltmp]
  #    Z2_qadj=Z2_qtmp[seltmp]
  # }
  #}

setsinenv=length(table(combn(datatmp$Sm,datatmp$Km)))
PbS=2
B=3
#C=2
C=Ncheck-1
while(all(PbS%%B!=0) & C<50){
      C=C+1
      N=Z1adj+((length(which(combn(datatmp$Sm,datatmp$Km)==1)))*Z2_qadj)+C
      dn=divisors(N)
      dn=dn[which(dn%in%5:20==TRUE)]
      if(length(dn)!=0){
        B=dn
        PbS=N*datatmp$R
      }
}
onlyuse=which((PbS%%B==0)==TRUE)
B=B[onlyuse]

#colsN must be multiple of block sizes
colsN=list()
rowsN=list()
if(length(B)!=0){
for (k in 1:length(B)){
  dPbS=divisors(PbS)
  colsN[[k]]=dPbS[which(dPbS%%B[k]==0)]
  rowsN[[k]]=PbS/colsN[[k]]
  
  if(any(colsN[[k]]==1)){
    idc1=which(colsN[[k]]%in%1==TRUE)
    colsN[[k]]=colsN[[k]][-idc1]
    rowsN[[k]]=rowsN[[k]][-idc1]
    }
  if(any(rowsN[[k]]==1)){
    idr1=which(rowsN[[k]]%in%1==TRUE)
    colsN[[k]]=colsN[[k]][-idr1]
    rowsN[[k]]=rowsN[[k]][-idr1]
  }
}
}

noset=dim(combn(datatmp$Sm,datatmp$Km))[2]-length(which(combn(datatmp$Sm,datatmp$Km)==1))
return(list(C,Z1adj,Z2_qadj,PbS,N,B,colsN,rowsN,noset))
}
alladj1=lapply(allinilist,adjparms)
alliniA=cbind(allini,C_adj=sapply(alladj1, "[[", 1),Z1_adj=sapply(alladj1, "[[", 2),Z2_qadj=sapply(alladj1, "[[", 3),PbS_adj=sapply(alladj1, "[[", 4),N_adj=sapply(alladj1, "[[", 5))

alliniB=cbind(allini,noset=sapply(alladj1, "[[", 9))

#TP:Total plots for complete experiment
TP_adj=alliniA$N_adj*alliniA$R*alliniA$Sm
tot_adj=alliniA$Z1_adj+(alliniA$Q-1)*alliniA$Z2_qadj
tot_adj1=alliniA$Z1_adj+(alliniA$Q-1)*alliniA$Z2_qadj+alliniA$C_adj
P_adj=round(((tot_adj1-alliniA$N_adj)/(tot_adj1))*100,2)
#P_adj=round(((alliniA$Z2_qadj*alliniB$noset)*100)/(tot_adj))#corregir (572-140)/572

alliniA=cbind(alliniA,TP_adj,tot_adj,P_adj)

completa=as.data.frame(matrix(0,1,21))
names(completa)=c("R","H","PbS","Sm","Km","P","Q","Z1","Z2_q","tot","C_adj","Z1_adj","Z2_qadj","PbS_adj","N_adj","TP_adj","tot_adj","P_adj","B_adj","cols_adj","rows_adj")
for ( k in 1:length(alladj1)){
ddtmp=list()
for(i in 1:length(alladj1[[k]][6][[1]])){
rows_adj=alladj1[[k]][8][[1]][[i]]
cols_adj=alladj1[[k]][7][[1]][[i]]
B_adj=rep(alladj1[[k]][6][[1]][i],times=length(cols_adj))
ddtmp[[i]]=cbind(B_adj,cols_adj,rows_adj)
}
ddtmp=do.call(rbind,ddtmp)
completa=rbind(completa,cbind(do.call("rbind", replicate(dim(ddtmp)[1], alliniA[k,], simplify = FALSE)),ddtmp))
}
completa=completa[-1,]
#Block size optional
if(!is.null(Bopt)){
  if(length(which(completa$B_adj==Bopt))!=0){
    completa=subset(completa, B_adj == Bopt)
  }else{
    completa=completa
    print("Don't found options availables with the block size desired")
  }
}

#Rows and columns optional
if(!is.null(colsNopt) && is.null(rowsNopt)){
  if(length(which(completa$cols_adj==colsNopt))!=0){
    completa=subset(completa, cols_adj == colsNopt)
  }else{
    completa=completa
    print("Don't found options availables with the columns number desired")
  }
}
if(is.null(colsNopt) && !is.null(rowsNopt)){
  if(length(which(completa$rows_adj==rowsNopt))!=0){
    completa=subset(completa, rows_adj == rowsNopt)
  }else{
    completa=completa
    print("Don't found options availables with the rows number desired")
  }
}
if(!is.null(colsNopt) && !is.null(rowsNopt)){
  if(length(which(completa$cols_adj == colsNopt & completa$rows_adj==rowsNopt))!=0){
    completa=completa[which(completa$cols_adj == colsNopt & completa$rows_adj==rowsNopt),]
  }else{
    completa=completa
    print("Don't found otions availables with the rows and columns number desired")
  }
}
if(length(which(completa$Z1_adj<ng))!=0){
	completa=completa[-which(completa$Z1_adj<ng),]
}

if(!is.null(Psh)){
	if(length(which(completa$P_adj==Psh))!=0){
		completa=completa[which(completa$P_adj==Psh),]
	}else{
		if(length(which(completa$P_adj>Psh))<dim(completa)[1] && length(which(completa$P_adj>Psh))!=0){	
			completa=completa[-which(completa$P_adj>Psh),]
		}
	}
}

if(!is.null(PbSsh)){
	if(length(which(completa$PbS_adj==PbSsh))!=0){
		completa=completa[which(completa$PbS_adj==PbSsh),]
	}else{
		if(length(which(completa$PbS_adj>PbSsh))<dim(completa)[1] && length(which(completa$PbS_adj>PbSsh))!=0){	
			completa=completa[-which(completa$PbS_adj>PbSsh),]
		}
	}
}

if(length(which(completa$tot_adj==Hfile))!=0){
	completa=completa[which(completa$tot_adj==Hfile),]
}else{
	if(length(which(completa$tot_adj<Hfile))!=0){
	completa=completa[-which(completa$tot_adj<Hfile),]
	}
}


completa=completa[,-c(2,3,6,8:10)]
rownames(completa)=as.character(1:dim(completa)[1])
names(completa)=c("R","Sm","Km","Q","C","Z1","Z2_q","PbS","N","TP","tot","P","B","cols","rows")
completa$BbS=(completa$PbS/completa$R)/completa$B

return(list(completa,dataCc,dataNCc,ng,minGen,Hfile))
#return(list(ng,minGen,P,Q,Z1,Z2_q,N,C,B,PbS,colsN,rowsN,dataCc,dataNCc,H))
}

