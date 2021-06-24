function [A,b] = treat_Dirichlet_2d_v2(A,b,boundarynodesU,boundarynodesP,Dirichlet_funU,Dirichlet_funV,Dirichlet_funP,PbU,PbP)
                                       
%-1:Dirichlet
%-2:Neumann
%-3:Robin

nbn1 = size(boundarynodesU,2);
nbn2 = size(boundarynodesP,2);

u1 = size(PbU,2);

%considering u at boundaries
for k = 1:nbn1
    
 if boundarynodesU(1,k) == -1
     
     i = boundarynodesU(2,k);
     
     A(i,:) = 0;
     A(i,i) = 1;
     
     b(i) = feval(Dirichlet_funU,PbU(1,i),PbU(2,i));
    
 end
end

%considering v at boundaries
for k = 1:nbn1
    
 if boundarynodesU(1,k) == -1
     
     i = boundarynodesU(2,k);
     
     A(i+u1,:) = 0;
     A(i+u1,i+u1) = 1;
     
     b(i+u1) = feval(Dirichlet_funV,PbU(1,i),PbU(2,i));
    
 end
end

%considering p at boundaries
for k = 1:nbn2
    
 if boundarynodesP(1,k) == -1
     
     i = boundarynodesP(2,k);
     
     A(i+u1+u1,:) = 0;
     A(i+u1+u1,i+u1+u1) = 1;
     
     if PbP(1,i) == 0 && PbP(2,i) == 0
        
         b(i+u1+u1) = feval(Dirichlet_funP,PbP(1,i),PbP(2,i));
        
     end
    
 end
end