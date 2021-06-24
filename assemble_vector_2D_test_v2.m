function b = assemble_vector_2D_test_v2(ceo1,ceo2,Pb,T_basis,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type_test,p,q)

Nb = size(Pb,2);

b = zeros(Nb,1);

for n = 1:number_of_elements
    
    vertices = P(:,T(:,n));
  
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    for beta = 1:number_of_local_basis
            
        result = Gauss_quadrature_2D_test_v2(ceo1,ceo2,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices,basis_type_test,p,q,beta);
            
        b(T_basis(beta,n),1) = b(T_basis(beta,n),1) + result; 
            
    end
    
end