function result = Gauss_quadrature_2D_test(ceo,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices_test,basis_type_test,der_order_test_x,der_order_test_y,basis_index_test)

Gpn = length(Gauss_coefficient_local_triangle);

result = 0;

for i = 1:Gpn
    
    result = result + Gauss_coefficient_local_triangle(i)*feval(ceo,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2))*local_basis_2D(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),vertices_test,basis_type_test,der_order_test_x,der_order_test_y,basis_index_test);
    
end