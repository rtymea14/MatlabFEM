function [P,T] = generate_P_T_2D_v2(left,right,bottom,top,h,basis_type)

N1 = (right-left)/h(1);
N2 = (top-bottom)/h(2);

if basis_type == 101
    
    i=0;
    for column = 1:(N1+1)
        for row = 1:(N2+1)
            
            i = i+1;
            
            P (1,i) = left + (column-1)*h(1);
            P (2,i) = bottom + (row-1)*h(2);
    
        end
    end
    
    for column = 1:N1
        for row = 1:N2
            
        element_index_1 = (column-1)*2*N2 + (row-1)*2 + 1;
        element_index_2 = element_index_1 + 1;
       
        node1 = (column-1)*(N2+1)+row;
        node2 = column*(N2+1)+row;
        node3 = node1 + 1;
        
        T(:,element_index_1) = [node1 node2 node3];
       
        node1_2 = node3;
        node2_2 = node2;
        node3_2 = node2_2 + 1;
       
        T(:,element_index_2) = [node1_2 node2_2 node3_2];
        
        end
    end
    
elseif basis_type == 102
    
    i=0;
    for column = 1:(2*N1+1)
        for row = 1:(2*N2+1)
            
            i = i+1;
            
            P (1,i) = left + (column-1)*h(1)/2;
            P (2,i) = bottom + (row-1)*h(2)/2;
    
        end
    end
    
    for column = 1:N1
        j = 0;
        for row = 1:N2
            
        element_index_1 = (column-1)*2*N2 + (row-1)*2 + 1;
        element_index_2 = element_index_1 + 1;
       
        node1 = 2*(column-1)*(2*N2+1)+row+j;
        node2 = 2*column*(2*N2+1)+row+j;
        node3 = node1 + 2;
        node4 = node2 - (2*N2+1);
        node5 = node4 + 1;
        node6 = node1 + 1;
        
        T(:,element_index_1) = [node1 node2 node3 node4 node5 node6];
       
        node1_2 = node3;
        node2_2 = node2;
        node3_2 = node2_2 + 2;
        node4_2 = node5;
        node5_2 = node2_2 + 1;
        node6_2 = node4_2 + 1;
       
        T(:,element_index_2) = [node1_2 node2_2 node3_2 node4_2 node5_2 node6_2];
        j = j+1;
        end
    end
    
end