######Function 1:-Computes the kth moment of r order statistic using Monte carlo Integration. #############

ordermean<-function(r,k,beta,n){
  A<-factorial(n)/((factorial(r-1))*(factorial(n-r)))
  integrand<- function(x) A*(x^k)*(((1-exp(beta-x))^(r-1))*((1+exp(beta))*exp(-x))^(n-r+1))/(1+exp(-x))^(n+1)
  v<- integrate(integrand,beta,Inf)
  v$value
}

###############################################################################################################

##########Function 2:- Computes cross product means of ordered statistics using Monte Carlo Integration #######

prodmeans<-function(r,s,beta,n) { 
  a<-factorial(n)/((factorial(r-1))*(factorial(s-r-1))*(factorial(n-s)))
  n2<-integrate(function(y) { 
    sapply(y, function(y) {
      integrate(function(x) {a*x*y*(((1-exp(-(x-beta)))*((1+exp(-x))^(-1)))^(r-1))*((((1-exp(-(y-beta)))*((1+exp(-y))^(-1)))-((1-exp(-(x-beta)))*((1+exp(-x))^(-1))))^(s-r-1))*((1-((1-exp(-(y-beta)))*((1+exp(-y))^(-1))))^(n-s))*((1+exp(beta))*exp(-x)*((1+exp(-x))^(-2)))*((1+exp(beta))*exp(-y)*((1+exp(-y))^(-2)))}, beta, y)$value
    })
  }, beta, Inf)
  
  n2$value
}

##################################################################################################################

############ Function 3:- Computes Mean Vector for actual set of n ordered statistics ###########################

muvector<-function(m,n,R,beta){
  mu<-matrix(0,n,1)
  for(i in 1:n){
    value<-ordermean(i,1,beta,n)
    mu[i,1]<-value
  }
  mu
}

##################################################################################################################

############## Function 4:- Computes the Variance of the actual set of n ordered statistics ######################

variance_vector<- function(m,n,R,beta){
  sig <- matrix(0,n,n)
  for (i in 1:n){
    for(j in i:n){
      if(i==j){
        value1<-ordermean(i,2,beta,n)
        value2<-ordermean(i,1,beta,n)
        value<- value1-(value2*value2)
        sig[i,j]<-round(value,4)
      }
      else{
        value1<-prodmeans(i,j,beta,n)
        value2<-ordermean(i,1,beta,n)
        value3<-ordermean(j,1,beta,n)
        value<-value1-(value2*value3)
        sig[i,j]<-round(value,4)
      }
    }
  }
  for(i in 1:n){
    for(j in 1:n){
      if(j<i){
        sig[i,j]<-sig[j,i]+10^(-10)
      }
    }
  }
  sig
}

####################################################################################################################


################## Function 5:- Computes K1,K2,....KM from the given R vector 
##################              for Thomas Wilson Probability computation ############################################

kvalues<-function(m,R){
  count<-0
  d<-data.frame()
  while(count<=10^4){
    count<-count+1
    for(i in 1:m){
      if(i==1){
        val<-i
        b<-val
      }
      else{
        s<-0
        for(k in 1:(i-1)){
          s<-s+R[k]
        }
        s<-s+i
        val<-sample(i:s,1)
        b<-c(b,val)
      }
    }
    z<-0
    for(i in 1:(length(b)-1)){
      if(b[i]<b[i+1]){
        z<-z+1
      }
    }
    if(z==m-1){
      d<-rbind(d,b)
    }
  }
  d<-unique(d)
  d
}
####################################################################################################################

############### Function 6:- Computes Matrix D which is of order mxn ###############################################
computematrix<-function(K,m,n){
  D<-matrix(0,m,n)
  for(r in 1:m){
    for(s in 1:n){
      if(s==K[r]){
        D[r,s]<-1
      }
    }
  }
  D
}

############### Function 7:- Computes probability associated with each D matrix using Thomas and Wilson Formula ########
computeprob<-function(K,R,n,m){
  p<-1
  for(i in 2:m){
    S<-0
    for(j in 1:(i-1)){
      S<-S+R[j]+1
    }
    value1<-choose(n-K[i],S-K[i]+1)
    value2<-choose(n-K[i-1],S-K[i-1])
    value<-value1/value2
    p<-p*value
  }
  p
}

############ Function 8:- Computes the mean vector to be used for computation of order mx1 ##########################
muvalue<-function(m,n,R,beta){
  mu<-muvector(m,n,R,beta)
  b<-kvalues(m,R)
  M<-length(b[,1])
  Smat<-matrix(0,m,n)
  for(i in 1:M){
    K<-as.numeric(b[i,])
    D<-computematrix(K,m,n)
    p<-computeprob(K,R,n,m)
    Smat<-Smat+D*p
  }
  meanvec<-Smat%*%mu
  meanvec
}

############# Function 9:- Computes the final Variance Covariance matrix to be used for computations of order mxm #######

sigvalues<-function(m,n,R,beta){
  b<-kvalues(m,R)
  M<-length(b[,1])
  Smat<-matrix(0,m,m)
  sig<-variance_vector(m,n,R,beta)
  mu<-muvector(m,n,R,beta)
  value<-sig+mu%*%t(mu)
  for(i in 1:M){
    K<-as.numeric(b[i,])
    D<-computematrix(K,m,n)
    p<-computeprob(K,R,n,m)
    Smat<-Smat+(D%*%value %*% t(D))*p
  }
  v1<-muvalue(m,n,R,beta)
  varvect<-Smat-(v1)%*%t(v1)
  varvect
}

#########################################################################################################

########## Function 10:- Computing Ci's #################################################################
Cvalues<-function(m,n,R,beta){
  sig<-sigvalues(m,n,R,beta)
  mu<-muvalue(m,n,R,beta)
  value1<-t(mu)%*%solve(sig)
  value2<-value1%*%mu
  s<-"The value of Var(sigma^2)/sigma^2 is:"
  print(paste(s,round(1/as.numeric(value2),5)))
  c<-value1/as.numeric(value2)
  print("The C_i vector is given as:-")
  c
}
#########################################################################################################
#########################################################################################################

########## Function 11:- Computing PHI #################################################################
phi<-function(A,B,S){
  v1<-t(A)%*%solve(S)
  v2<-v1%*%B
  v2
     }
##################################################################################################################################################################################################################
#########################################################################################################

########## Function 12:- Computing PHY #################################################################
phy<-function(A,B,C,S){
  v1<-t(A)%*%solve(S)
  v2<-v1%*%B
  v3<-t(C)%*%solve(S)
  v4<-v2%*%v3
  v4
}
##################################################################################################################################################################################################################
#########################################################################################################

########## Function 13:- Computing ai's #################################################################
avalues<-function(m,n,R,beta){
  sig<-sigvalues(m,n,R,beta)
  mu<-muvalue(m,n,R,beta)
  one<-matrix(1,m,1)
  v1<-phy(mu,mu,one,sig)
  v2<-phy(mu,one,mu,sig)
  v3<-v1-v2
  v4<-phi(mu,mu,sig)
  v5<-phi(one,one,sig)
  v6<-phi(mu,one,sig)
  v7<-v4%*%v5-v6%*%v6
  a<-v3/as.numeric(v7)
  
  s<-"The value of Var. of estimate is:"
  print(paste(s,round(as.numeric(v4)/as.numeric(v7),5)))
  
  print("The a_i vector is given as:-")
  a
  }
#########################################################################################################
#########################################################################################################

########## Function 14:- Computing bi's #################################################################
bvalues<-function(m,n,R,beta){
  sig<-sigvalues(m,n,R,beta)
  mu<-muvalue(m,n,R,beta)
  one<-matrix(1,m,1)
  v1<-phy(one,one,mu,sig)
  v2<-phy(one,mu,one,sig)
  v3<-v1-v2
  v4<-phi(mu,mu,sig)
  v5<-phi(one,one,sig)
  v6<-phi(mu,one,sig)
  v7<-v4%*%v5-v6%*%v6
  b<-v3/as.numeric(v7)
  s<-"The value of Var. of estimate is:"
  print(paste(s,round(as.numeric(v5)/as.numeric(v7),5)))
  s1<-"The value of Cov. of estimate is:"
  print(paste(s,round(as.numeric(-v6)/as.numeric(v7),5)))
  print("The b_i vector is given as:-")
  b
  
}
#########################################################################################################

