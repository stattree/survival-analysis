#library(NLRoot)

#### Finding root by bisection method #######
BFfzero1 = function (f, a, b, num = 10, eps = 1e-08) 
{
  h = abs(b - a)/num
  i = 0
  j = 0
  a1 = b1 = 0
  while (i <= num) {
    a1 = a + i * h
    b1 = a1 + h
    if (f(a1) == 0) {
      print(a1)
      print(f(a1))
    }
    else if (f(b1) == 0) {
      print(b1)
      print(f(b1))
    }
    else if (f(a1) * f(b1) < 0) {
      repeat {
        if (abs(b1 - a1) < eps) 
          break
        x <- (a1 + b1)/2
        if (f(a1) * f(x) < 0) 
          b1 <- x
        else a1 <- x
      }
#      print(j + 1)
      j = j + 1
 #     print((a1 + b1)/2)
 #     print(f((a1 + b1)/2))
      return((a1+b1)/2)
    }
    i = i + 1
  }
  if (j == 0) 
  return(-1000)    # return -1000 in failing case
    #  print(-1000)
}


##### parameter setup #############

alpha = 0.3

beta = function(s)
{
  exp(-s)
  #2
  #0.01*s+1
  #sin(pi*s)
}

##### f ######
f = function(s)
{
  #0.2*s*I()
  #2*dnorm(s)*I(s>0)
  exp(-s)
  #1*(abs(s)<1)
  #exp(-s)*(abs(s)<3)/(1-exp(-3))
  #2*s*(abs(s)<1)
}

# integrand
int <- function(x0,s)
{
  exp(x0*beta(s))*f(s)
}


# conditional survival function
Stx = function(x,t)
{
   exp(- exp(alpha)*integrate(int,lower=0,upper=t,x0=x)$value)
#  inc = 0.01
#  g = seq(0,t,inc)
#  exp( - sum(inc*exp(alpha+x*beta(g))*f(g)))
}



est_beta <- function(X_b,Z_b,D_b,alpha_b,f_b,h_b)
{
  
  s_T = sort(Z_b[D_b==1])
  
  kerf <- function(x,h)
  {
    dnorm(x/h)
  }
  
  LL_f <- function(t1,b)
  {
    ll1=ll2=c()
    
    for (i in 1:length(X_b))
    {
      ll1[i] = (D_b[i]!=2)*exp(alpha_b+X_b[i]*b)*sum( (s_T<Z_b[i])*f_b*kerf(s_T-t1,h_b) )
      ll2[i] = (D_b[i]==2)*exp(alpha_b+X_b[i]*b)*sum( f_b*kerf(s_T-t1,h_b) )
    }
    
    ll =  - ( sum((D_b==1)*b*X_b*kerf(Z_b-t1,h_b)) - sum(ll1+ll2) )
    return(ll)
  }
  
  hatb = c()
  
  for (j in 1:length(s_T))
  {
    hatb[j] = optimize(LL_f,t1=s_T[j],lower=0,upper=5)$minimum
  }
  return(hatb)
}
######################################
###### estimation of F (p) #########
######################################




est_f <- function(X_f,Z_f,D_f,alpha_f,beta_f)
{
  int_2 <- function(x0,s)
  {
    exp(x0*beta(s))
  }
  
  
  int_w <- function(X0,T0,D0)
  {
    s_T = sort(T0[D0==1])
    part = c(0,s_T)
    int0 = matrix(NA,length(s_T),length(X0))
    for (k in 1:length(s_T))
    {
      for (j in 1:length(X0))
      {
        int0[k,j] = integrate(int_2,lower=part[k],upper=part[k+1],x0=X0[j])$value
      }
    }
    return(int0)
  }
  
  
  s_T = sort(T[D_f==1])      ############# event times #############
  oldeta = eta0 = 2/n0
  #eP = Dis_P1(eta0,T,D,X)
  part = c(0,s_T)
  
  c_m = int_w(X_f,T,D_f)/diff(part)  ############ numerical integration ##############
  
  
  Dis_P2 <- function(eta,T0,D0,X0)
  {
    p=c()
    s_T = sort(T0[D0==1])
    for (k in seq(1,length(s_T)-1))
    {
      temp =  
        (D0!=2)*(s_T[k]<=T0)*exp(X0*beta_f[k])
      + (D0==2)*exp(X0*beta_f[k])
      - (D0!=2)*(eta<=T0)*exp(X0*beta_f[length(s_T)])
      - (D0==2)*exp(X0*beta_f[length(s_T)])
      p[k] = 1/(1/eta + exp(alpha_f)*sum(temp))
    }
    return(append(p,eta))
  }
  
  
  ###### Newton raphson ###########
  for (k in 1:100)
  {
    eP = Dis_P2(oldeta,T,D_f,X_f)
    neweta = oldeta - ( sum(eP)-1 )/( (sum(eP[-length(eP)]^2))/oldeta^2 + 1  )
    if (abs(neweta-oldeta)<1e-6) break
    oldeta = neweta
    #sum(Dis_P1(oldeta,T,D,X))
  }
  
  eP = Dis_P2(neweta,T,D_f,X_f)
  return(eP)
  #temp_f = approx(s_T,cumsum(eP), method="constant",yleft=0,yright=1,xout=eg)$y
  #return(temp_f)
}


est_alpha <- function(X_a,Z_a,D_a,beta_a,f_a)
{
  s_T = sort(Z_a[D_a==1])

  s0 = c()  
  
  ll1=ll2=c()
  
  for (i in 1:length(X_a))
  {
    ll1[i] = (D_a[i]!=2)*sum( exp(X_a[i]*beta_a)*(s_T<Z_a[i])*f_a )
    ll2[i] = (D_a[i]==2)*sum( exp(X_a[i]*beta_a)*f_a )
    s0[i] = ll1[i]+ll2[i]
  }
  
  return(log(N/sum(s0)))
}


################ integration ##########################



#s_T = sort(T[D==1])
#eP = Dis_P1(0.2,T,D,X)
#plot(s_T,cumsum(eP))
#lines(s_T,1-exp(-s_T))

#eg = seq(0,max(s_T),length=200)
#est = approx(s_T,cumsum(eP), yleft=0,xout=eg)
#plot(est$x,est$y)
#lines(eg,1-exp(-eg))
  
##### random number generation ###########
set.seed(3)
n0 = 200   # sample size


Nrep = 2 # number of replication
tg = seq(0,5,0.05)

h=0.35

hata=rep(NA,Nrep)   # estimated alpha
hat_beta = matrix(0,Nrep,length(tg))  
#eg = seq(0,3,length=200)
#ef = matrix(NA,Nrep,length(eg)) 
ef = list(NA)

a=proc.time()
  
for (rep in 1:Nrep)
{

  X = runif(n0)    # covariate
  U = runif(n0)    # for inverse method
  T = c()          # event time 
  px = c()         # cure probability
  
  # generate event times
  for (i in 1:n0)
  {
      px[i] = Stx(X[i],100)   # cure probability for i subject
      # define a function to apply inverse method
      BS = function(t){( (1 - Stx(X[i],t))/(1-(px[i])))-U[i]}
      T[i] = BFfzero1(BS, 0, 50)
  }
  
  T[runif(n0)<px] <- 999   # cured case

  T = T[T!=-1000]          # exclude the failure case of finding root
  X = X[T!=-1000]
  
  n = length(T)            # effective sample size
  
  C = runif(n,0,3)      # censoring variable
  Z = pmin(C,T)         # observed time
  D = 1*(T<C)           # censoring indicator
  D[T==999] = 2         # cure indicator (D=2)
  
  N0 = sum(D==0)    # number of censoring
  N = sum(D==1)     # number of event
  N2 = sum(D==2)    # number of cure
  

####### estimation of f when other parameters are known #############

alpha_ini = alpha_o = 0.2
beta_ini = beta_o = rep(1,N)
#f_ini = f_o = f

for (iter in 1:10)
{

f_n = est_f(X,Z,D,alpha_o,beta_o)
#f_n = stepfun(eg,c(f_n[1],f_n))

beta_n = est_beta(X,Z,D,alpha_o,f_n,h)
#beta_n = stepfun(tg,c(beta_n[1],beta_n))

alpha_n = est_alpha(X,Z,D,beta_n,f_n)

f_o = f_n 
beta_o = beta_n 
alpha_o = alpha_n 



print(c(rep,iter))
print(alpha_n)
print(sum(f_n^2))
}

ef[[rep]] = f_n
#hat_beta[rep,] = beta_n
hata[rep] = alpha_n
#print(rep)

}

proc.time()-a


plot(eg,colMeans(ef),ylim=c(0,1.5))
lines(eg,1-exp(-eg),col=2)


