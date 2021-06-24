function boundarynodes = generate_boundary_nodes_info(left,right,bottom,top,h,basis_type,T,bc)

N1 = (right-left)/h(1);
N2 = (top-bottom)/h(2);


if basis_type == 101
    boundarynodes = zeros(2,N1+N2+N1+N2);
    %bottom edge

    for i = 1:N1
    
        if i == 1
        
            boundarynodes(1,i) = -1;
            boundarynodes(2,i) = (i-1)*(N2+1)+1;
        
        end
        
        boundarynodes(1,i) = bc;
        boundarynodes(2,i) = (i-1)*(N2+1)+1;
       
    end

    %right edge

    for i = 1:N2
    
        boundarynodes(1,N1+i) = -1;
        boundarynodes(2,N1+i) = (N2+1)*N1+i;
     
    end
    %top edge

    for i = N1:-1:1
    
        boundarynodes(1,N1+N2+N1-i+1) = -1;
        boundarynodes(2,N1+N2+N1-i+1) = (N2+1)*(i+1);
     
    end
    %left edge
    for i = N2:-1:1
    
        boundarynodes(1,N1+N2+N1+N2-i+1) = -1;
        boundarynodes(2,N1+N2+N1+N2-i+1) = i+1;
     
    end
    
elseif basis_type == 102
    boundarynodes = zeros(2,2*(N1+N2+N1+N2));
    %bottom edge

    for i = 1:2*N1
    
        if i == 1
            
            boundarynodes(1,i) = -1;
            boundarynodes(2,i) = (i-1)*(2*N2+1)+1;
            
        end
        
        boundarynodes(1,i) = bc;
        boundarynodes(2,i) = (i-1)*(2*N2+1)+1;
       
    end

    %right edge

    for i = 1:2*N2
    
        boundarynodes(1,2*N1+i) = -1;
        boundarynodes(2,2*N1+i) = (2*N2+1)*2*N1+i;
     
    end
    %top edge

    for i = 2*N1:-1:1
    
        boundarynodes(1,2*N1+2*N2+2*N1-i+1) = -1;
        boundarynodes(2,2*N1+2*N2+2*N1-i+1) = (2*N2+1)*(i+1);
     
    end
    %left edge
    for i = 2*N2:-1:1
    
        boundarynodes(1,2*N1+2*N2+2*N1+2*N2-i+1) = -1;
        boundarynodes(2,2*N1+2*N2+2*N1+2*N2-i+1) = i+1;
     
    end
    
end