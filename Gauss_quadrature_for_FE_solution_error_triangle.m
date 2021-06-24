function r = Gauss_quadrature_for_FE_solution_error_triangle(uh_local,exact_function,vertices,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_order_x,derivative_order_y)

Gpn=length(Gauss_coefficient_reference_triangle);

r=0;

[Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);

for i=1:Gpn
    
    r=r+Gauss_coefficient_local_triangle(i)*(feval(exact_function,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2))-FE_solution_triangle(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),uh_local,vertices,basis_type,derivative_order_x,derivative_order_y))^2;

end