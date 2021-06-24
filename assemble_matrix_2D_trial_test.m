function A = assemble_matrix_2D_trial_test(ceo,Pb,T_basis,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type_trial,r,s,basis_type_test,p,q)

Nb = size(Pb,2);

A = sparse(Nb,Nb);

for n = 1:number_of_elements
    
    vertices = P(:,T(:,n));
         
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    for alpha = 1:number_of_local_basis
        
        for beta = 1:number_of_local_basis
            
            result = Gauss_quadrature_2D_trial_test(ceo,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,basis_type_trial,r,s,alpha,vertices,basis_type_test,p,q,beta);
                                               
            A(T_basis(beta,n),T_basis(alpha,n)) = A(T_basis(beta,n),T_basis(alpha,n)) + result;
            
        end
        
    end
    
end