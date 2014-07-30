#test approach for sampling species totals directly

lambda=c(16,64)

N1=16
N2=64

pi11=0.5
pi21=0.25
pi1=c(pi11,pi21)

c11=0.75
c21=0.5
c1=c(c11,c21)

#expected value data
X111=N1*pi11*c11
X121=N1*(1-pi11)*c11
X112=N1*pi11*(1-c11)
X122=N1*(1-pi11)*(1-c11)
X211=N2*pi21*c21
X221=N2*(1-pi21)*c21
X212=N2*pi21*(1-c21)
X222=N2*(1-pi21)*(1-c21)

O=matrix(0,2,2)
O[1,1]=sum(X111,X211)
O[1,2]=sum(X112,X212)
O[2,1]=sum(X121,X221)
O[2,2]=sum(X122,X222)

X=array(0,dim=c(2,2,2,1000))
for(i in 1:1000){
  for(io in 1:2){
    for(ic in 1:2){
      if(ic==1)curc=c1
      else curc=1-c1
      if(io==1)curpi=pi1
      else curpi=1-pi1
      pars=c(lambda*curc*curpi)
      X[,io,ic,i]=rmultinom(1,O[io,ic],pars)
    }
  }
}

crap=apply(X,4,'mean')

