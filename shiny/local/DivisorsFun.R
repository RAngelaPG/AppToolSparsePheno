onlyDiv=function(M,N){
  library(numbers)
  dn=list()
  for (i in 1:M){
    dn[[i]]=divisors(N[i])
    dn[[i]]=dn[[i]][which(dn[[i]]%in%5:20==TRUE)]
  }
return(dn)
}