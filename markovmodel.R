

# here we load the functions required for the simulation


source("NeighborhoodDataWD.R");
source("Net_Import.R");
source("NetCmbn.R");
source("Post_NodesState.R");
source("GEMF_SIM.R");
source("Para_switching.R");


# here we load two layer networks in a format that the main simulation function uses
N=500;#number of nodes
Net1=Net_Import("L1er.txt",N);
Net2=Net1;
Net3=Net_Import("L2er.txt",N);
NetSet=list(Net1,Net2,Net3);
Net=NetCmbn(NetSet,N);
rm(Net1,Net2,Net3,NetSet);


#here we set the model parameters(rates) and arranged them in a format that is used by the main simulation
p2=0.5;
p0=0.5;
bet=0.24;
gam1=5;
gam2=gam1*((1/p2)-1);
Para=Para_switching(delta=1,gamma_1=gam1,gamma_2=gam2,bet,betprim=bet*p0)
M<-Para[[1]];


#we set initial sates of the nodes
x0=matrix(as.integer(4),1,N);


#simulation terminator 
maxNumevent=1000000;Runtime=150;


# here we run the simulation 2 times
numrun=2;
dd={};
for (i in 1:numrun){
  lst<-GEMF_SIM(Para,Net,x0,maxNumevent,Runtime,N);
  dd[[i]]=lst;
  print(i);
}





timstp=0.5;
Tr<-seq(0,Runtime,timstp);
EI=matrix(numeric(0), numrun,length(Tr)) 

  
  for(runn in 1:numrun){
    
    ts<-dd[[runn]][[1]];
    n_index<-dd[[runn]][[2]];
    i_index<-dd[[runn]][[3]];
    j_index<-dd[[runn]][[4]];
    Tf<-dd[[runn]][[5]];
    lasteventnumber<-dd[[runn]][[6]]; 
    lst3<-Post_NodesState(x0,M,N,ts,n_index,j_index,lasteventnumber,timstp,Runtime)
    nodstt=lst3[[2]];
    regI=rowSums(nodstt>=3);
    EI[runn,]=regI;
    
  };  



matplot(Tr, colMeans(EI)/N, type = "l", lty = 1,lwd=2,  xlab="Time",ylab="Infection Prevalence",col= rainbow(M))
legend("topright",legend="",lty = 1,lwd=2,col= rainbow(M))  



