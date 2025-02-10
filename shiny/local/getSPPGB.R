genSPPGB<-function(rowsN,colsN,PbS,Sm,Km,dataCc,C,dataNCc,Z1,ng,Z2_q,Q,B,N,R,name,sumdat){
  library(xlsx)
  library(rfunctions)
  library(Matrix)
  library(MASS)
  
  if(any(rowsN*colsN!=PbS)){
    print("The arrange row*col is not possible")
  }else{
    RowCols=list()
    rc=1
    mit=rowsN[rc]%%2
    if(mit!=0){
      Ccol=rep(c(rep(c(1:colsN[rc],colsN[rc]:1),floor(rowsN[rc]/2)),c(1:colsN[rc])),Sm[rc])
    }else{
      Ccol=rep(rep(c(1:colsN[rc],colsN[rc]:1),(rowsN[rc]/2)),Sm[rc])
    }
    Crow=rep(rep(1:rowsN[rc],each=colsN[rc]),Sm[rc])
    RowCols[[rc]]=as.data.frame(cbind(Row=Crow,Col=Ccol))
  }
  RowCols=do.call("rbind",RowCols)
  
  if(dim(dataCc)[1]<C){
	print("Not enought Checks names in the input file, we will generated in automatically")
	Cadd=C-dim(dataCc)[1]
	forrep=dataCc[dim(dataCc)[1],]
	for (i in 1:Cadd){
		forreptmp=apply(forrep,2,as.character)
		forreptmp=as.data.frame(t(paste(forreptmp,"_Add",i,sep="")))
		names(forreptmp)=names(dataCc)
		forreptmp$Role="Check"
		forreptmp$Group="All"		
		forreptmp$SetNum="All"
		dataCc=rbind(dataCc,forreptmp)
	}
  
  }
  dataCc=dataCc[1:C,]
  
  aldesall=list()
  #M=1
  #rm(setname)
  dataNC=dataNCc
  dataC=dataCc
  data_splitZ=list()
  
  
  if(Z1%%ng==0){
    data_splitZ[[1]]=unlist(lapply(split(dataNC$IDnum,dataNC$Group),function(x) sample(x,Z1/ng)))
    dataNC$SetNum[data_splitZ[[1]]]="Set1"
    datatmp=dataNC[-which(dataNC$IDnum%in%data_splitZ[[1]]==TRUE),]
  }else{
    splitmp=unlist(lapply(split(dataNC$IDnum,dataNC$Group),function(x) sample(x,floor(Z1/ng))))
    datatmp=dataNC[-which(dataNC$IDnum%in%splitmp==TRUE),]
    rest1=Z1-(ng*floor(Z1/ng))
    splitmp2=sample(datatmp$IDnum,rest1)
    datatmp=datatmp[-which(datatmp$IDnum%in%splitmp2==TRUE),]
    data_splitZ[[1]]=c(splitmp,splitmp2)
    dataNC$SetNum[data_splitZ[[1]]]="Set1"
  }
  setname="Set1"
  
  for( j in 1:(Q-2)){
    data_splitZ[[j+1]]=datatmp$IDnum[1:Z2_q]
    dataNC$SetNum[data_splitZ[[j+1]]]=paste("Set",j+1,sep="")
    datatmp=datatmp[-which(datatmp$IDnum%in%data_splitZ[[j+1]]==TRUE),]
    setname=c(paste("Set",j+1,sep=""),setname)
  } 
  
  dataNC$SetNum[datatmp$IDnum[1:Z2_q]]=paste("Set",Q,sep="")
  setname=c(paste("Set",Q,sep=""),setname)
  setname=setname[-Q]
  
  rm(datatmp)
  setnameg=setname[length(setname):1]
  setsinenv=combn(Sm,Km)
  setsbyenv=list()
  for (i in 1:Sm){
    setsbyenv[[i]]=c("Set1",setnameg[which(apply(setsinenv,2,function(x){i%in%x})==TRUE)])
  }
  #generate alpha(0,1) for each ENV
  aldesfin=list()
  EffEnv=list()
  for(j in 1:Sm){
    aldes=alfaDRA(k=B,v=N,r=R, pp=1)
    EffEnv[[j]]=c(round(aldes$eteo,5),round(aldes$E,5))
    aldes=aldes$design
    aldes=aldes[order(aldes$Entry),]
    
    datatmpenv=rbind(dataNC[which(dataNC$SetNum%in%setsbyenv[[j]]),!names(dataNC)%in%c("IDnum")],dataC)
    #datatmpenv=rbind(dataNC[-which(dataNC$SetNum==setname[j]),!names(dataNC)%in%c("IDnum")],dataC)
    datatmpenv=datatmpenv[order(datatmpenv$GID),]
    datatmpenv$IDNum=1:dim(datatmpenv)[1]
    datatmpenv=datatmpenv[rep(seq_len(nrow(datatmpenv)), each=R),]
    aldes2=cbind(aldes,datatmpenv)
    aldesfin[[j]]=cbind(Managment=rep(1,N*R),Site=rep(j,N*R),aldes2[,!names(aldes2)%in%c("Entry","IDNum")])
    rm(aldes,datatmpenv,aldes2)
  }
  
  aldesall[[1]]=do.call("rbind", aldesfin)    
  rm(aldesfin)  
  
  #For calculate efficiency
  aldesall=do.call("rbind", aldesall)  
  aldesall=aldesall[order(aldesall$Managment,aldesall$Site,aldesall$Plot),]
  aldesall=cbind(aldesall,RowCols)
  #aldesall$Site=as.factor(aldesall$Site)
  #aldesall$GID=as.factor(aldesall$GID)
  #aldesall$Block=as.factor(aldesall$Block)
  #aldesall$Rep=as.factor(aldesall$Rep)
  #Xtmp=sparse.model.matrix(~ GID, aldesall)
  #Ztmp=sparse.model.matrix(~ Site:Rep + Site + Site:GID + Site:Rep:Block, aldesall)
  #X=matrix(0,dim(Xtmp)[1],dim(Xtmp)[2])
  #Xx=as.data.frame(summary(Xtmp))
  #X[as.matrix(Xx[,1:2])]<-Xx[,3]
  #Z=matrix(0,dim(Ztmp)[1],dim(Ztmp)[2])
  #Zz=as.data.frame(summary(Ztmp))
  #Z[as.matrix(Zz[,1:2])]<-Zz[,3]
  #V=crossprodcpp(t(Z),diag(dim(Z)[1]))+diag(dim(X)[1])
  #invV=geninv(V)
  #A1=t(X)%*%invV%*%X
  #EffAll=round(sum(diag(A1))/(dim(Z)[2]+dim(X)[1]),5)

  
  EffEnvM=do.call("rbind",EffEnv)
  #EffEnvM=rbind(EffEnvM,c("",EffAll))
  colnames(EffEnvM)=c("TheoreticalEff","ReachedEff")
  rownames(EffEnvM)=c(apply(as.data.frame(1:Sm),1,function(x){paste("Env_",x,sep="")}))
  #rownames(EffEnvM)=c(apply(as.data.frame(1:Sm),1,function(x){paste("Env_",x,sep="")}),"AllSparse")
  EffEnvM=as.data.frame(EffEnvM)
  
  #### Create workbook ####
  wb <- createWorkbook(type="xlsx")
  sheet <- createSheet(wb, sheetName = "Summary")
  xlsx.addTable(wb, sheet, data= sumdat, col.names=TRUE, row.names=FALSE, columnWidth=30,
                fontColor="darkblue", fontSize=12, rowFill="white")
  
  sheet <- createSheet(wb, sheetName = "Efficiency")
  xlsx.addTable(wb, sheet, data= EffEnvM, col.names=TRUE, row.names=TRUE, columnWidth=15,
                fontColor="darkblue", fontSize=12, rowFill="white")
  
  sheet <- createSheet(wb, sheetName = "FieldBook")
  xlsx.addTable(wb, sheet, data= aldesall, col.names=TRUE, row.names=FALSE, columnWidth=10,
                fontColor="darkblue", fontSize=12, rowFill="white")
  
  ##Create fieldMap by Managment
  nlevMag=unique(aldesall$Managment)
  MapMag<-aldesall[which(aldesall$Managment==nlevMag),]
  nameSheet <- paste("FieldMap_Manag_",nlevMag,sep="")
  sheet <- createSheet(wb, sheetName = nameSheet)
  
  nlevSite=unique(MapMag$Site)
  for (j in 1:Sm){
    MapMagEnv<-MapMag[which(MapMag$Site==nlevSite[j]),]
    MapField <- reshape(MapMagEnv[,c("GID","Row","Col")],v.names = "GID", idvar = "Row", timevar = "Col", direction = "wide")
    row.names(MapField) <- MapField[,1]
    MapField <- MapField[,-1]
    names(MapField) <- gsub("GID.","",names(MapField))
    
    xlsx.addHeader(wb, sheet, value=paste("Field Map Entries in Site ",nlevSite[j],sep=""),level=6, startCol=3, color="darkblue", underline=1, backGroundColor="white")
    xlsx.addLineBreak(sheet, 1)
    xlsx.addTableBlocks(wb, sheet, data = MapField[nrow(MapField):1,], col.names=TRUE, row.names=TRUE, columnWidth=7,
                        fontColor="darkblue", fontSize=12, rowFill="white", nband=1, nstack=B)
    xlsx.addLineBreak(sheet, 2)
  }
  
  saveWorkbook(wb, name)
  return("yes")
}
