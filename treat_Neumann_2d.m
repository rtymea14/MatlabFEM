function v = treat_Neumann_2d(cfun,pfun,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,boundary_edges,local_basis_test,P,T,Tb,vector_size,basis_type_test,der_order_test_x,der_order_test_y)

nbe = size(boundary_edges,2);

v = sparse(vector_size,1);

for k = 1:nbe
    
    if boundary_edges(1,k) == -2
              
        nk = boundary_edges(2,k);
        
        end_point_1 = P(:,boundary_edges(3,k));
        
        end_point_2 = P(:,boundary_edges(4,k));
        
        [Gauss_nodes_2d_line,Gauss_weights_2d_line] = generate_2d_line_Gauss(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2);
        
        vertices = P(:,T(:,nk));
        
        for beta = 1:local_basis_test
            
            temp = Gauss_quadrature_2D_line_test(cfun,pfun,Gauss_weights_2d_line,Gauss_nodes_2d_line,vertices,basis_type_test,der_order_test_x,der_order_test_y,beta);
                
            v(Tb(beta,nk),1) = v(Tb(beta,nk),1) + temp; 
            
        end
        
    end
    
end