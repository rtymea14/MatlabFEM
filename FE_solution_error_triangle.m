function r=FE_solution_error_triangle(uh,exact_function,left,right,bottom,top,h,basis_type,derivative_order_x,derivative_order_y,Gauss_point_number)

N1=(right-left)/h(1);
N2=(top-bottom)/h(2);

number_of_elements=2*N1*N2;

[P,T]=generate_P_T_2D_v2(left,right,bottom,top,h,101);

if basis_type==101
    
    Pb = P;
    Tb = T;
    
elseif basis_type==102
    
    [Pb,Tb]=generate_P_T_2D_v2(left,right,bottom,top,h,102);
    
end

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

r=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    
    vertices=P(:,T(:,n));
    
    uh_local=uh(Tb(:,n));
    
    r = r +Gauss_quadrature_for_FE_solution_error_triangle(uh_local,exact_function,vertices,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,derivative_order_x,derivative_order_y);

end

r=sqrt(r);