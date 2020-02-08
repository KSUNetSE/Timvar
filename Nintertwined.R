

# here we load the functions required for the simulation
source("NeighborhoodDataWD.R");
source("Net_Import.R");
source("NetCmbn.R");
source("Para_switching.R");
source("GEMF_ODE.R");


# here we load two layer networks in a format that the main solver function uses
N=500;
Net1=Net_Import("L1er.txt",N);
Net2=Net1;
Net3=Net_Import("L2er.txt",N);
NetSet=list(Net1,Net2,Net3);
Net=NetCmbn(NetSet,N);
rm(Net1,Net2,Net3,NetSet);

#here we set the model parameters(rates) and arranged them in a format that is used by the main solver
p2=0.5;
p0=0.5;
bet=0.24;
gam1=5;
gam2=gam1*((1/p2)-1);
Para=Para_switching(delta=1,gamma_1=gam1,gamma_2=gam2,bet,betprim=bet*p0)
M<-Para[[1]];


#we set initial sates of the nodes
x0=matrix(as.integer(4),1,N); 

P0=matrix(0,M,N);
for(i in 1:N ){
  P0[x0[i],i]=1;
}

#here we solve the N-intertwined equation
TT=10;# final time
sol=GEMF_ODE(Para,Net,P0,TT); # solving ODE


#here we plot the results
t=sol[[1]];  # time points that we have mean field solution for them.
P_t=sol[[2]]; #P_t is a list of M (number of state) elements. each element is a 
# n_t*N matrix (number of time points* number of nodes). for example P_t[[m]] 
#is a matrix that stores probablity of finding nodes in state m at 
#different time points.  Probablity of finding node n in state m at 
#time point j, is  p_t[[m]][j,n]

#ploting population of each comportmant through time 
Comp_P_t=matrix(0,size(t)[2],M)
leg=character(M)

for (i in 1:M){
  Comp_P_t[,i]=rowSums(P_t[[i]]);
  
  leg[i]= paste("state",toString(i));
}


matplot(t, Comp_P_t, type = "l", lty = 1,lwd=2,  xlab="time",ylab="states' population",col= rainbow(M))
legend("topright",legend=leg,lty = 1,lwd=2,col= rainbow(M))  


