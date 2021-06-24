function A = assemble_matrix_2D_trial_test_v2(ceo,PbU,T_basisU,PbP,T_basisP,P,T,number_of_elements,number_of_local_basis1,number_of_local_basis2,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type_trial,r,s,basis_type_test,p,q)

Nb1 = size(PbU,2);
Nb2 = size(PbP,2);

A = sparse(Nb1,Nb2);

for n = 1:number_of_elements
    
    vertices = P(:,T(:,n));
         
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    for alpha = 1:number_of_local_basis1
        
        for beta = 1:number_of_local_basis2
            
            result = Gauss_quadrature_2D_trial_test(ceo,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,basis_type_trial,r,s,alpha,vertices,basis_type_test,p,q,beta);
                                               
            A(T_basisU(beta,n),T_basisP(alpha,n)) = A(T_basisU(beta,n),T_basisP(alpha,n)) + result;
            
        end
        
    end
    
end