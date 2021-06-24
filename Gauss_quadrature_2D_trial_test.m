function result = Gauss_quadrature_2D_trial_test(ceo_fun,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices_trial,basis_type_trial,der_order_trial_x,der_order_trial_y,basis_index_trial,vertices_test,basis_type_test,der_order_test_x,der_order_test_y,basis_index_test)

Gpn = length(Gauss_coefficient_local_triangle);

result = 0;

for i = 1:Gpn
    
    result = result + Gauss_coefficient_local_triangle(i)*feval(ceo_fun,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2))*local_basis_2D(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),vertices_trial,basis_type_trial,der_order_trial_x,der_order_trial_y,basis_index_trial)*local_basis_2D(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),vertices_test,basis_type_test,der_order_test_x,der_order_test_y,basis_index_test);
    
end