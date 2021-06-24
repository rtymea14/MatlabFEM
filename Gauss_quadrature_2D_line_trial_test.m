function result = Gauss_quadrature_2D_line_trial_test(cfun,qfun,Gauss_coefficient_local_triangle,Gauss_point_local_triangle,vertices_trial,basis_type_trial,der_order_trial_x,der_order_trial_y,basis_index_trial,vertices_test,basis_type_test,der_order_test_x,der_order_test_y,basis_index_test)

Gpn = length(Gauss_coefficient_local_triangle);

result = 0;

for i = 1:Gpn
    
    result = result + Gauss_coefficient_local_triangle(i)*feval(cfun,Gauss_point_local_triangle(1,i),Gauss_point_local_triangle(2,i))*feval(qfun,Gauss_point_local_triangle(1,i),Gauss_point_local_triangle(2,i))*local_basis_2D(Gauss_point_local_triangle(1,i),Gauss_point_local_triangle(2,i),vertices_trial,basis_type_trial,der_order_trial_x,der_order_trial_y,basis_index_trial)*local_basis_2D(Gauss_point_local_triangle(1,i),Gauss_point_local_triangle(2,i),vertices_test,basis_type_test,der_order_test_x,der_order_test_y,basis_index_test);
    
end