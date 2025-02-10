#rm(list=ls())
#ls()
# Alfa design 
#
# INPUTS
# k  - length of longer block (6)
# v  - number of varieties (24)
# r  - number of replications (4)
# pp - number of rows per replication (1)
#
# OUTPUTS

# E_teo - theoretical eficiency of the design
# E     - reached eficiency of the design
# design - Field Book - Experimental design
# Note that both Xpuv and Xpl include NaNs to keep structure of the design, 
# where all blocks have length k. This full structure allows identification
# of individual blocks and replications when evaluating the experiment.

alfaDRA=function(k,v,r,pp)
{
	pocetGenerovani=50; # number of created designs in searching design with sufficient efficiency
	delka_d_pasu=v/pp; # length of row
	ad=((delka_d_pasu%%1)==0)&((ceiling(delka_d_pasu)%%k)==0); # alpha design indicator
	
	# generation of alpha design
	if(ad)
	{
		s=v/k;
		H=matrix(0,k,r)
		G_ext=matrix(0,k,r*s)
		E=0
		E_teo=(v-1)*(r-1)/((v-1)*(r-1)+r*(s-1)) # theoretical efficiency
		for(l in 1:pocetGenerovani)
		{
			rand = matrix(runif((k-1)*(r-1)),k-1,r-1)
			H[2:k,2:r]=floor(s*rand)
			# calculation of alpha design efficiency
			C=matrix(0,v,v)
			for(a in 1:(k-1))
			{
				for(b in (a+1):k)
				{
					R = (H[b,]-H[a,])%%s
					eye = matrix(0,s,s)
					diag(eye) = 1
					Blok=sum(R==0)*eye
					for(dd in 1:(s-1))
					{
						ones = matrix(1,s-dd,s-dd)
						parte1=matrix(0,s,s)
						parte2=matrix(0,s,s)
						if(dd == s-1)
						{
						  parte1[1,s]=diag(sum(dd==R)*ones)
						  parte2[s,1]=diag(sum(s-dd==R)*ones)
						}
						else
						{
						  parte1[1:(s-dd),(dd+1):s]=diag(diag(sum(dd==R)*ones))
						  parte2[(dd+1):s,1:(s-dd)]=diag(diag(sum(s-dd==R)*ones))
						}
						Blok=Blok+parte1+parte2
					}
					C[(((a-1)*s)+1):(a*s),(((b-1)*s)+1):(b*s)]=Blok
				}
			}
			I=diag(1,v)
			CM=r*I+C+t(C)
			A=I-1/r/k*CM
			vl_cisla=eigen(A)$values
			nen_vl_cisla=vl_cisla[1:(v-1)]
			h=1/nen_vl_cisla
			y=(v-1)/sum(h)
			# calculation of alpha design efficiency
			if (y > E) # searching for designs with higher efficiency
			{
				E=y
				G=H
			}
		}
		# rotation of blocks
		Q=matrix(1,k,1)%*%matrix(0:(r-1),1)
		G=G+Q%%s
		vec=matrix(1,k,r)
		for(i in 1:s)
		{
		G_ext[,seq(i,(r*s),s)]=G+(i-1)*vec
		}
		G_ext=G_ext%%s
		vec=matrix(1,1,(r*s))
		for(i in 1:k-1)
		{
			G_ext[i+1,]=G_ext[i+1,]+i*s*vec
		}
		P=array(0,c(k,s,r))
		for(i in 1:r)
		{
			P[,,i]=G_ext[,(1+((i-1)*s)):(i*s)]
		}
		P=P+1
	X=t(matrix(P,nrow=v,byrow=F))[r:1,]
	X2=X[r:1,]
	X3=t(matrix (t(X2),nrow=v/pp,byrow=F))
	Xpl=X3[r:1,]
	Xpuv=Xpl
	}

# Aleatoriza Genotipos (xgen) a números en Xpl
xgen=sample(seq(1:v), v, replace = FALSE, prob = NULL)
XXX=matrix(0,nrow(Xpl),ncol(Xpl))
for(i in 1:nrow(Xpl)){
  XXX[i,]=xgen[Xpl[i,]]
}

Entry=as.vector(t(XXX))

Rep=rep(sample(seq(1:r), r, replace = FALSE, prob = NULL),each=v)

BB=list()
for(qq in 1:r){
  BB[[qq]]=rep(sample(seq(1:(v/k)), v/k, replace = FALSE, prob = NULL),each=k)
}
Block=unlist(BB, use.names=FALSE)

# Aleatoriza tratamientos dentro de bloques
rando=runif(v*r)
design1=cbind(Rep,Block,Entry,rando)
design2=design1[order(Rep,Block,rando),-4]

# Agrega dato de plot
Plot=seq(1,v*r)
design3=as.data.frame(cbind(Plot,design2))

	lista = list(eteo=E_teo,E=E,design=design3)
	return(lista)
}

#a=alfaDRA(k=4,v=8,r=4, pp=1)

