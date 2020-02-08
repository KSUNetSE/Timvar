GEMF_ODE<-function(Para,Net,P0,Tf){
  
  library(pracma)
  
  
  M=Para[[1]]; q=Para[[2]]; L=Para[[3]]; A_d=Para[[4]]; A_b=Para[[5]];
  
  Neigh=Net[[1]]; I1=Net[[2]]; I2=Net[[3]]; N=dim(I1)[2];
  
  Q_d=diag(rowSums(A_d))-A_d;
  
  
  Q_b=list();
         for(i in 1:L){
             Q_b[[i]]=diag(rowSums(A_b[,,i]))-A_b[,,i];
            } ;
  
  P0_vec=matrix(P0,M*N,1);
 
  
#------definition of function 
  
  Meanfield<-function(t,P_vec){
  
    P=matrix(P_vec,M,N);
    y=matrix(0,N,L);
    dP=matrix(0,M,N);
    
    for (i in 1:N){
      for (l in 1:L){
        if(I1[l,i]!=0){     
          
            Nei=Neigh[[l]][1,I1[l,i]:I2[l,i]];
            Wei=Neigh[[l]][2,I1[l,i]:I2[l,i]];
            y[Nei,l]=y[Nei,l]+Wei*P[q[l],i];
          
        };
      };
    };
    
    
    for (i in 1:N){
      RN=-t(Q_d);
      for(l in 1:L){
           RN=RN-y[i,l]*t(Q_b[[l]]);
      };
      dP[,i]=RN %*% P[,i]
    
    };
    
    dP_vec=matrix(dP,M*N,1)
    
  }
    
#end of function------    
    
    
  sol = ode45(Meanfield, 0, Tf, P0_vec, hmax = 0.1);
  
  t=sol$t;
  y=sol$y;
 
   P_t=list();
  
  for(i in 1:M){
    
    P_t[[i]]=y[,seq(i,M*N,M)];
  };
  
  list(t,P_t);
  
   }

