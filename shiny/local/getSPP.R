genSPP=function(datos,M,Q,Z1,Z2_q,B,N,Sm,R,RowCols,name)
{
library(xlsx)
#Once the calculus was defined, provided the file with GID,group,role
datos$Group=as.factor(datos$Group)
datos$Role=as.factor(datos$Role)
datos$GID=as.factor(datos$GID)
datos=datos[order(datos$Group),]

dataCc=datos[which(datos$Role=="Check"),]
dataCc$SetNum=rep("All",dim(dataCc)[1])
dataNCc=datos[-which(datos$Role=="Check"),]
dataNCc$Group=as.factor(as.character(dataNCc$Group))
dataNCc$Role=as.factor(as.character(dataNCc$Role))
dataNCc$GID=as.factor(as.character(dataNCc$GID))
dataNCc$IDnum=1:dim(dataNCc)[1]
dataNCc$SetNum=rep("Set",dim(dataNCc)[1])
rm(datos)  
ng=nlevels(dataNCc$Group)
aldesall=list()
for(i in 1:M){
  dataNC=dataNCc
  dataC=dataCc
  data_splitZ=list()
  data_splitZ[[1]]=unlist(lapply(split(dataNC$IDnum,dataNC$Group),function(x) sample(x,Z1[i]/ng)))
  dataNC$SetNum[data_splitZ[[1]]]="Set1"
  datatmp=dataNC[-which(dataNC$IDnum%in%data_splitZ[[1]]==TRUE),]
  setname="Set1"
  for( j in 1:(Q[i]-2)){
    data_splitZ[[j+1]]=unlist(lapply(split(datatmp$IDnum,datatmp$Group),function(x) sample(x,Z2_q[i]/ng)))
    dataNC$SetNum[data_splitZ[[j+1]]]=paste("Set",j+1,sep="")
    datatmp=datatmp[-which(datatmp$IDnum%in%data_splitZ[[j+1]]==TRUE),]
    setname=c(paste("Set",j+1,sep=""),setname)
  }
  dataNC$SetNum[datatmp$IDnum]=paste("Set",Q[i],sep="")
  setname=c(paste("Set",Q[i],sep=""),setname)
  setname=setname[-Q[i]]
  #This is because maybe we dont have the same number of hybrids in each group and is not enough for include in each set
  if(dim(datatmp)[1]<Z2_q[i]){
    how=Z2_q[i]-dim(datatmp)[1]
    for(k in 1:how){
      addf=as.data.frame(t(c(rep(paste("Zaux_",k,sep=""),(dim(dataNC)[2]-1)),paste("Set",Q[i],sep=""))))
      names(addf)=names(dataNC)
      addf$Role="Hybrid"
      dataNC=rbind(dataNC,addf)
    }
  }
  rm(datatmp)
  #generate alpha(0,1) for each ENV
  aldesfin=list()
  for(j in 1:Sm[i]){
    aldes=alfaDRA(k=B[i],v=N[i],r=R, pp=1)$design
    aldes=aldes[order(aldes$Entry),]
    
    datatmpenv=rbind(dataNC[-which(dataNC$SetNum==setname[j]),!names(dataNC)%in%c("IDnum")],dataC)
    datatmpenv=datatmpenv[order(datatmpenv$GID),]
    datatmpenv$IDNum=1:dim(datatmpenv)[1]
    datatmpenv=datatmpenv[rep(seq_len(nrow(datatmpenv)), each=R),]
    aldes2=cbind(aldes,datatmpenv)
    aldesfin[[j]]=cbind(Managment=rep(i,N[i]*R),Site=rep(j,N[i]*R),aldes2[,!names(aldes2)%in%c("Entry","IDNum")])
    rm(aldes,datatmpenv,aldes2)
  }
  
  aldesall[[i]]=do.call("rbind", aldesfin)    
  rm(aldesfin)  
  cat("Managment:",i,"of",M,"\n")
  rm(setname)  
}

aldesall=do.call("rbind", aldesall)  
aldesall=aldesall[order(aldesall$Managment,aldesall$Site,aldesall$Plot),]
aldesall=cbind(aldesall,RowCols)

#### Create workbook ####
wb <- createWorkbook(type="xlsx")
sheet <- createSheet(wb, sheetName = "FieldBook")
xlsx.addTable(wb, sheet, data= aldesall, col.names=TRUE, row.names=FALSE, columnWidth=10,
              fontColor="darkblue", fontSize=12, rowFill="white")

##Create fieldMap by Managment
nlevMag=unique(aldesall$Managment)
for(i in 1:M){
  MapMag<-aldesall[which(aldesall$Managment==nlevMag[i]),]
  nameSheet <- paste("FieldMap_Manag_",nlevMag[i],sep="")
  sheet <- createSheet(wb, sheetName = nameSheet)
  
  nlevSite=unique(MapMag$Site)
  for (j in 1:Sm[i]){
    MapMagEnv<-MapMag[which(MapMag$Site==nlevSite[j]),]
    MapField <- reshape(MapMagEnv[,c("GID","Row","Col")],v.names = "GID", idvar = "Row", timevar = "Col", direction = "wide")
    row.names(MapField) <- MapField[,1]
    MapField <- MapField[,-1]
    names(MapField) <- gsub("GID.","",names(MapField))
    
    xlsx.addHeader(wb, sheet, value=paste("Field Map Entries in Site ",nlevSite[j],sep=""),level=6, startCol=3, color="darkblue", underline=1, backGroundColor="white")
    xlsx.addLineBreak(sheet, 1)
    xlsx.addTableBlocks(wb, sheet, data = MapField[nrow(MapField):1,], col.names=TRUE, row.names=TRUE, columnWidth=7,
                        fontColor="darkblue", fontSize=12, rowFill="white", nband=1, nstack=B[i])
    xlsx.addLineBreak(sheet, 2)
  }
}

#setwd(savefile)
saveWorkbook(wb, name)
return("yes")
}