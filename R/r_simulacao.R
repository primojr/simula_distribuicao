# Gerando valores

rm(list=ls(all=TRUE))

n = 25
B = 5

#lambda=0.90, alpha=0.5,theta=2.5, a=1.5,b=0.50; b0=-0.5;b1=0.7

lambda = 0.9
alpha = 0.5
theta = 2.5
a = 1.50
b = 0.5
beta0 = -0.5 # Intercepto
beta1 = 0.7 # Coef. da regressÃ£o

cova = phi = z = N = u = v = te = ce = tempo = delta = matrix(nrow=n,ncol=B)

cova = sapply(1:B, function(x) rbinom(n,1,.5))
phi = sapply(1:B, function(x) exp(beta0+beta1*cova))
z = sapply(1:B, function(x) rpois(n, lambda=lambda))
N = sapply(1:B, function(x) rpois(n, phi[,B]*z[,B]))

# Variavel condicional
for (j in 1:B){
  for(i in 1:n){
    if (N[i,j]>0) {u[i,j]=max(runif(N[i,j],0,1))} else {u[i,j] = Inf}
  }}
#
v = sapply(1:B, function(x) ((1/(((1/u[,B])-1)^(1/alpha)+1))^(1/(1/theta))))
te = sapply(1:B, function(x) qweibull(v[,B], shape=a, scale=b))
ce = sapply(1:B, function(x) runif(n,0,15))  # tau=15, este valor controla a % de censura
#tempo = sapply(1:B, function(x) min(te[,B], ce[,B]))
delta = sapply(1:B, function(x) ifelse(te[,B] < ce[,B], 1, 0)) # indicador de censura
for(j in 1:B){
  for(i in 1:n){
    tempo[i,j] = min(te[i,j], ce[i,j])
  }}


#phi.real=exp(beta0+beta1*cova)
#p.real=exp(lambda*exp(-phi.real)-lambda)
#p.real[cova==1]
#mean(p.real[cova==1])
#p.real[cova==0]
#mean(p.real[cova==0])


##############################
######## ESTIMAÇÃO OPTIM #####
##############################

b0.max=numeric()
b0.max=matrix(nrow=1,ncol=B)

b1.max=numeric()
b1.max=matrix(nrow=1,ncol=B)


phi2=numeric()
phi2=matrix(nrow=n,ncol=B)

lambda.max=numeric()
lambda.max=matrix(nrow=1,ncol=B)

alpha.max=numeric()
alpha.max=matrix(nrow=1,ncol=B)

theta.max=numeric()
theta.max=matrix(nrow=1,ncol=B)

a.max=numeric()
a.max=matrix(nrow=1,ncol=B)

b.max=numeric()
b.max=matrix(nrow=1,ncol=B)

gt2=numeric()
gt2=matrix(nrow=n,ncol=B)

Gt2=numeric()
Gt2=matrix(nrow=n,ncol=B)

ft2=numeric()
ft2=matrix(nrow=n,ncol=B)

Ft2=numeric()
Ft2=matrix(nrow=n,ncol=B)

St2=numeric()
St2=matrix(nrow=n,ncol=B)



dens=numeric()
dens=matrix(nrow=n,ncol=B)

sob=numeric()
sob=matrix(nrow=n,ncol=B)

max_vero=numeric()
max_vero=matrix(nrow=n,ncol=B)

f_vero=numeric()
f_vero=matrix(nrow=n,ncol=B)

log_vero=numeric()
log_vero=matrix(nrow=n,ncol=B)

soma_log_vero=numeric()
soma_log_vero=matrix(nrow=1,ncol=B)

resu_max_vero_maximo=numeric()
resu_max_vero_maximo=matrix(nrow=1,ncol=B)

convergencia_maximo=numeric()
convergencia_maximo=matrix(nrow=1,ncol=B)

for (i in 1:B){
  maximo=function(x){

    b0=x[1]
    b1=x[2]
    lambda2=x[3]
    alpha2=x[4]
    theta2=x[5]
    a2=x[6]
    b2=x[7]

    gt2[,i]=(a2*tempo[,i]^(a2-1)/(b2^a2))*(exp(-(tempo[,i]/b2)^a2))
    Gt2[,i]=(1-exp(-(tempo[,i]/b2)^a2))
    ft2[,i]=(alpha2*(1/theta2)*gt2[,i]*(Gt2[,i]^(alpha2*(1/theta2)-1))*(1-Gt2[,i]^(1/theta2))^(alpha2-1))/(((Gt2[,i]^(alpha2*(1/theta2)))+(1-Gt2[,i]^(1/theta2))^(alpha2))^2)
    Ft2[,i]=(Gt2[,i]^(alpha2*(1/theta2)))/(((Gt2[,i]^(alpha2*(1/theta2)))+(1-Gt2[,i]^(1/theta2))^alpha2))
    St2[,i]=(1-Ft2[,i])
    phi2[,i]=exp(b0+b1*cova[,i])
    dens[,i]=lambda2*phi2[,i]*(exp(lambda2*(exp(-phi2[,i]*St2[,i])-1)))*exp(-phi2[,i]*St2[,i])*ft2[,i]
    sob[,i]=(1+exp(lambda2*(exp(-phi2[,i])-1))-exp(lambda2*(exp(-phi2[,i]*St2[,i])-1)))
    max_vero[,i]=(dens[,i]^(delta[,i]))*(sob[,i]^(1-delta[,i]))
    log_vero[,i]=log(max_vero[,i])
    soma_log_vero[,i]=-sum(log_vero[,i])
    return(soma_log_vero[,i])
  }

  chute=c(beta0,beta1,lambda,alpha,theta,a,b)
  resultados=optim(chute,fn=maximo,method="BFGS",hessian=FALSE)

  b0.max[i]=resultados$par[1]
  b1.max[i]=resultados$par[2]
  lambda.max[i]=resultados$par[3]
  alpha.max[i]=resultados$par[4]
  theta.max[i]=resultados$par[5]
  a.max[i]=resultados$par[6]
  b.max[i]=resultados$par[7]

  resu_max_vero_maximo[i]=maximo(c(resultados$par[1],resultados$par[2],resultados$par[3],resultados$par[4],resultados$par[5],resultados$par[6],resultados$par[7]))
  convergencia_maximo[i]=resultados$convergence

}


#__________RESULTADOS #__________________
resultados_optim_max=matrix(c(lambda.max,alpha.max,theta.max,a.max,b.max,b0.max,b1.max,resu_max_vero_maximo,convergencia_maximo),nrow=B,ncol=9)
dados_max_optim=data.frame(resultados_optim_max)
names(dados_max_optim)=c("lambda","alpha","theta", "a","b","b0.max","b1.max","vero","convergencia")

#________Valores estimados
dados_max_optim

#________Média das estimativas
media_optim_max=c(mean(dados_max_optim$lambda),
                  mean(dados_max_optim$alpha),
                  mean(dados_max_optim$theta),
                  mean(dados_max_optim$a),
                  mean(dados_max_optim$b),
                  mean(dados_max_optim$b0.max),
                  mean(dados_max_optim$b1.max)
)
media_optim_max

