library('randtoolbox')
library('doSMP')
library('Brobdingnag')
library('matrixStats')

#############################################################################
## smp code

rmSessions(qnames=NULL, all.names=TRUE)
w = startWorkers(workerCount = 4)
registerDoSMP(w)

#############################################################################
## Generate data - simulate mixture of .5*N(10,2)+.5*N(30,5)

y.n=100
y=c(rep(0,y.n))
for(i in 1:y.n){
	temp=runif(1)
	if(temp < .5){
		y[i]=rnorm(1,10,2)
	}
	else {
		y[i]=rnorm(1,30,5)
	}
}
ybar=sum(y)/y.n
ss=sum((y-ybar)^2)

##
###################################################################################
## calculate variance using more efficient estimate

my.var=function(x){return(sum((x-mean(x))^2/(length(x)+1)))} #for Dr. bratcher ;-)

##
###################################################################################
## Density function specification

f.3=function(para,data){
	eval=matrix(0,length(para[,1]),length(data))
	
	for(j in 1:length(data)){
		eval[,j]=log(para[,5]*dnorm(data[j],mean=para[,1],sd=para[,3])+(1-para[,5])*dnorm(data[j],mean=para[,2],sd=para[,4]))
	}
	return(rowSums(eval))		
}

##
###################################################################################

###################################################################################
## Get.seq -> returns a new random sequence in the quasi-MC-chain
## int -> start a new chain?
## eval-> return the new sequence for all parameters

get.seq=function(limits,n,p,int=FALSE){
	eval=NULL
	new.seq=halton(n,p,init=int)
	for(j in 1:p){
		eval=cbind(eval,qunif(new.seq[,j],limits$lower[j],limits$upper[j]))
	}
	return(eval)
}
#
###################################################################################
## Enforces order for location parameters -> VERY IMPORTANT

enforce.loc=function(grid,loc.par=NULL){
	n=length(loc.par)
	
	grid.loc=foreach(i =1:n,.combine='cbind') %do%{
		grid[,loc.par[i]]
	}
	for(i in 1:n){
		col.max=max.col(grid.loc)
		grid[,loc.par[n-i+1]]=rowCollapse(grid.loc,col.max)
		grid.loc[col(grid.loc)==col.max]=-Inf
	}
	return (grid)
}

##
###################################################################################
## calculates the mean/var of the estimates to deterime the new bounds
##

cal.stats=function(
		theta)
{
	eval = foreach(i = 1:dim(theta)[1],.combine='rbind') %do% {
				cbind(mean(theta[i,]),sqrt(my.var(theta[i,])))
			}
	return(eval)
}

##
###################################################################################
## Integration processing of pincus theorem

pincus.integrate=function(
		grid,
		lam,
		fn,
		data,
		enforce.order=NULL,
		loc.par=NULL)
{
	p 			 = dim(grid)[2]
	numer 		 = c(rep(0,p))
	estimate 	 = c(rep(0,p))
	den 		 = 0
	est.brob	 = NULL
	
	t			 = enforce.order(grid,loc.par=loc.par)
	f.nag		 = exp(lam*as.brob(fn(t,data)))
	den			 = den+sum(f.nag[f.nag!=as.brob(+exp(-Inf))])
	
	for(k in 1:p){			
		temp	 = sum(t[f.nag!=as.brob(+exp(-Inf)),k]*f.nag[f.nag!=as.brob(+exp(-Inf))])
		temp2	 = numer[k]
		est.brob = cbind(est.brob,(temp+temp2)/den)
	}
	for(k in 1:p){
		estimate[k] = as.numeric(est.brob[[k]])
	}
	return (estimate)
}
##
###################################################################################
##

get.bounds=function(para.stats,limits,loc.par=NULL,scale.par=NULL,p.par=NULL){
	ll = limits$lower
	ul = limits$upper
	
	if(!is.null(loc.par)){
		ll[loc.par] = para.stats[loc.par,1]-3*para.stats[loc.par,2]
		ul[loc.par] = para.stats[loc.par,1]+3*para.stats[loc.par,2]
	}
	if(!is.null(scale.par)){
		ll[scale.par] = foreach(i = 1:length(scale.par),.combine='cbind') %do% 
				{
					max(0,para.stats[scale.par[i],1]-3*para.stats[scale.par[i],2])
				}
		ul[scale.par] = para.stats[scale.par,1]+3*para.stats[scale.par,2]
	}
	if(!is.null(p.par)){
		ll[p.par] = foreach(i = 1:length(p.par),.combine='cbind') %do%
				{
					max(para.stats[p.par,1]-3*para.stats[p.par,2],0)
				}
		ul[p.par] = foreach(i =1:length(p.par),.combine='cbind') %do%
				{
					min(para.stats[p.par,1]+3*para.stats[p.par,2],1)
				}
	}
	return(list(lower=ll,upper=ul))
}

###################################################################################
##  Begin integration processing (need to do some optimization here)
pincus=function(
	data		= NULL,			#model data
	lambda		= 10000,		#lambda to be used in pincus theorem
	seq.n		= 10,			#number of values to sample per parameter
	seq.p		= 5,			#number of parameters
	bounds		= NULL,			#boundry limits for parameters
	chunk_size	= 10000,		#size to split up integration (depends on memory)
	fn			= NULL,			#function to be evaluated
	.precision	= .01,			#convergence precision
	.loc.par 	= NULL,			#column of parameter list of any location parameters
	.scale.par 	= NULL,
	.p.par 		= NULL
	)
{
	grid.size 	= seq.n^seq.p
	Theta.est	= matrix(0,seq.p,seq.n)
	loop		= TRUE
	loop.num	= 0
	while(loop){ 
		for(i in 1:10){
			next.seq 	= get.seq(bounds,seq.n,seq.p,int=if(i==1){TRUE}else{FALSE}) 
			tmp 		= apply(next.seq,2,list) #these 2 lines will dynamically build the grid on entered parameters
			Theta 		= lapply(tmp,unlist)
			iter_chunk 	= iter(data.matrix(expand.grid(Theta)), chunksize = chunk_size,by = 'row')
			out 		= foreach(chunk = iter_chunk,
									.packages 		= c('matrixStats','Brobdingnag'),
									.export 		= c('enforce.loc','pincus.integrate'),
									.combine 		= 'rbind',
									.multicombine 	= TRUE) %dopar% 
							{  
								pincus.integrate(chunk, data,fn = fn, lam = lambda, enforce.order = enforce.loc,loc.par = .loc.par)
							}
			if(grid.size > chunk_size){
				Theta.est[,i] = colMeans(out)
			}else{
				Theta.est[,i] = out
			}
		}
		parameter.stats = cal.stats(Theta.est)
		bounds 			= get.bounds(parameter.stats,bounds,loc.par=.loc.par,scale.par=.scale.par,p.par=.p.par)
		loop 			= max(bounds$upper-bounds$lower)>.precision
		loop.num 		= loop.num+1
		print(bounds)
		if(loop.num > 10){
			loop = FALSE
		}
	}
	return (parameter.stats)
}

pincus(data=y,fn=f.3,bounds=list(lower=c(-100,-100,0,0,0),upper=c(100,100,30,30,1)),
		.precision=.1,.loc.par=c(1,2),.scale.par=c(3,4),.p.par=c(5))
