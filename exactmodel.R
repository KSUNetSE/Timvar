
# here we load the functions required for the simulation

source("NeighborhoodDataWD.R");
source("Net_Import.R");
source("GEMF_SIM_several_switching_exact_rev.R");
source("Para_switching_exact.R");



# here we load two layer networks in a format that the main simulation function uses

N=500;#N is number of nodes in the network
Net1=Net_Import("L1er.txt",N);
Net2=Net_Import("L2er.txt",N);
v1=Net1[[2]]; v2=Net1[[3]];
IN1=rbind(v1,v2); Nei1=Net1[[1]][[1]];

v1=Net2[[2]]; v2=Net2[[3]];
IN2=rbind(v1,v2); Nei2a=Net2[[1]][[1]];

Nei2b=Nei2a*0;
for (n in 1:N){
  if(IN2[1,n]!=0){
    for(ind in IN2[1,n]:IN2[2,n]){
      neigbor<-Nei2a[1,ind];
      
      neisofneighbor<-Nei2a[1,IN2[1,neigbor]:IN2[2,neigbor]];
      corind=(which(neisofneighbor == n)[[1]])-1+IN2[1,neigbor];
      Nei2b[1,ind]<-corind;
    }
    
  }
  
}
Nei2=rbind(Nei2a,Nei2b)
Nei2[4,]<-1;
rm(v1,v2,Net1,Net2,Nei2a,Nei2b)



#we set initial sates of the nodes
M=4;
x0=matrix(as.integer(4),1,N);



#simulation terminator 
maxNumevent=1000000;Runtime=20;


#here we set the model parameters(rates) and arranged them in a format that is used by the main simulation
p2=0.5;
p0=0.5;
bet=0.24;
gam1=5;
gam2=gam1*((1/p2)-1);
Para=Para_switching_exact(delta=1,gamma_1=gam1,gamma_2=gam2,beta1=bet,beta2=bet,p0);



# here we run the simulation 200 times and save the result as an R file
lst<- GEMF_SIM_several_switching_exact_rev(Para,IN1,IN2,Nei1,Nei2,x0,maxNumevent,Runtime,N,numrun=2,res=1);  


Tr=lst[[1]]; 
Ipop=lst[[5]]+lst[[4]];


matplot(Tr, rowMeans(Ipop)/N, type = "l", lty = 1,lwd=2,  xlab="Time",ylab="Infection Prevalence",col= rainbow(M))
legend("topright",legend="",lty = 1,lwd=2,col= rainbow(M))  



