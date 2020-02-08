Para_switching_exact<-function(delta,gamma_1,gamma_2,beta1,beta2,p0){
  M=4;#number of compartments
  

  
  
  A_d=matrix(0,M,M);# node base transition matrix
  A_d[1,2]<-gamma_1; A_d[3,4]<-gamma_1;
  A_d[2,1]<-gamma_2; A_d[4,3]<-gamma_2;
  A_d[3,1]<-delta; A_d[4,2]<-delta;
  

  
 
  list(A_d,beta1,beta2,p0);
}