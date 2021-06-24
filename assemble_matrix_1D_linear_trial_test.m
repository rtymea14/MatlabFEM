function A = assemble_matrix_1D_linear_trial_test(ceo,Pb,T_basis,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type_trial,der_order_trial,basis_type_test,der_order_test)

Nb = size(Pb,2);

A = sparse(Nb,Nb);

for n = 1:number_of_elements
    
    vertices = Pb(:,T_basis(:,n));
  
    lower_bound=min(vertices(1),vertices(2));
    upper_bound=max(vertices(1),vertices(2));
         
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    
    for alpha = 1:number_of_local_basis
        
        for beta = 1:number_of_local_basis
            
            r = Gauss_quadrature_1D_trial_test(ceo,Gauss_coefficient_local_1D,Gauss_point_local_1D,vertices,basis_type_trial,der_order_trial,alpha,vertices,basis_type_test,der_order_test,beta);
            
            A(T_basis(beta,n),T_basis(alpha,n)) = A(T_basis(beta,n),T_basis(alpha,n)) + r;
            
        end
        
    end
    
end