function r=FE_solution_triangle(x,y,uh_local,vertices,basis_type,derivative_order_x,derivative_order_y)

r=0;

number_of_local_basis=length(uh_local);

for i=1:number_of_local_basis
    
    r=r+uh_local(i)*local_basis_2D(x,y,vertices,basis_type,derivative_order_x,derivative_order_y,i);
    
end