function boundaryedges = generate_boundary_edge_info(left,right,bottom,top,h,basis_type,T,bc)

N1 = (right-left)/h(1);
N2 = (top-bottom)/h(2);
boundaryedges = zeros(4,N1+N2+N1+N2);

%bottom edge

for k = 1:N1
    
    if k == 1
        
        boundaryedges(1,k) = -1;
        boundaryedges(2,k) = 2*N2*(k-1)+1;
        boundaryedges(3,k) = T(1,boundaryedges(2,k));
        boundaryedges(4,k) = T(2,boundaryedges(2,k));
        
    end
    
    boundaryedges(1,k) = bc;
    boundaryedges(2,k) = 2*N2*(k-1)+1;
    boundaryedges(3,k) = T(1,boundaryedges(2,k));
    boundaryedges(4,k) = T(2,boundaryedges(2,k));
    
end

%right edge

for k = 1:N2
    
     boundaryedges(1,N1+k) = -1;
     boundaryedges(2,N1+k) = 2*N2*(N1-1)+2*k;
     boundaryedges(3,N1+k) = T(2,boundaryedges(2,N1+k));
     boundaryedges(4,N1+k) = T(3,boundaryedges(2,N1+k));

end
%top edge

for k = N1:-1:1
    
     boundaryedges(1,N1+N2+N1-k+1) = -1;
     boundaryedges(2,N1+N2+N1-k+1) = 2*N2*k;
     boundaryedges(3,N1+N2+N1-k+1) = T(3,boundaryedges(2,N1+N2+N1-k+1));
     boundaryedges(4,N1+N2+N1-k+1) = T(1,boundaryedges(2,N1+N2+N1-k+1));
end
%left edge
for k = N2:-1:1
    
     boundaryedges(1,N1+N2+N1+N2-k+1) = -1;
     boundaryedges(2,N1+N2+N1+N2-k+1) = 2*k-1;
     boundaryedges(3,N1+N2+N1+N2-k+1) = T(3,boundaryedges(2,N1+N2+N1+N2-k+1));
     boundaryedges(4,N1+N2+N1+N2-k+1) = T(1,boundaryedges(2,N1+N2+N1+N2-k+1));
end