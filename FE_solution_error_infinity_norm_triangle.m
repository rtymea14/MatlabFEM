function r=FE_solution_error_infinity_norm_triangle(uh,exact_function,left,right,bottom,top,h,basis_type,derivative_order_x,derivative_order_y,Gpn)

N1=(right-left)/h(1);
N2=(top-bottom)/h(2);

number_of_elements=2*N1*N2;

[P,T]=generate_P_T_2D_v2(left,right,bottom,top,h,101);

if basis_type == 101
    
    Pb = P;
    
    Tb = T;
    
elseif basis_type == 102
    
    [Pb,Tb]=generate_P_T_2D_v2(left,right,bottom,top,h,102);
    
end

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle] = generate_Gauss_reference_triangle(Gpn);

r=0;

for n=1:number_of_elements
    
    vertices = P(:,T(:,n));
    
    uh_local = uh(Tb(:,n));
    
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    s = size(Gauss_point_local_triangle,1);
    
    for i = 1:s
    
        temp(i) = max(abs(feval(exact_function,Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2))-FE_solution_triangle(Gauss_point_local_triangle(i,1),Gauss_point_local_triangle(i,2),uh_local,vertices,basis_type,derivative_order_x,derivative_order_y)));    
    
    end
    
    temp = max(temp);
    
    if temp>r
        
        r=temp;
        
     end
       
 end