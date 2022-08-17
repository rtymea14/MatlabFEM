%function [error_on_nodes,infinity_norm_error,L2_error,H1_error] = Steady_2D_elliptic(left,right,bottom,top,h,basis_type,Gauss_point_number,bc)

function [error_on_nodes,infinity_norm_error,L2_error,H1_error] = Steady_2D_elliptic(left,right,bottom,top,h,Gauss_point_number,bc)

%obtain necessary parameter

[P,T] = generate_P_T_2D_v2(left,right,bottom,top,h,101);
%boundaryedges = generate_boundary_edge_info(left,right,bottom,top,h,101,T,bc);

  
PbP = P;
TbP = T;
boundarynodesP = generate_boundary_nodes_info(left,right,bottom,top,h,101,T,bc);
    
    
[PbU,TbU] = generate_P_T_2D_v2(left,right,bottom,top,h,102);
boundarynodesU = generate_boundary_nodes_info(left,right,bottom,top,h,102,T,bc);
  
number_of_elements = size(T,2);

number_of_local_basisP = 3;

number_of_local_basisU = 6;

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number);

a1 = assemble_matrix_2D_trial_test('function_2mu',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,1,0,102,1,0);

a2 = assemble_matrix_2D_trial_test('function_mu',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,0,1,102,0,1);

A11 = a1 + a2;

a3 = assemble_matrix_2D_trial_test('function_mu',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,0,1,102,1,0);

A12 = a3;

a4 = assemble_matrix_2D_trial_test_v2('function_minusone',PbU,TbU,PbP,TbP,P,T,number_of_elements,number_of_local_basisP,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,101,0,0,102,1,0);

A13 = a4;

A21 = A12.';

a5 = assemble_matrix_2D_trial_test('function_mu',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,1,0,102,1,0);

a6 = assemble_matrix_2D_trial_test('function_2mu',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,0,1,102,0,1);

A22 = a5 + a6;

a7 = assemble_matrix_2D_trial_test_v2('function_minusone',PbU,TbU,PbP,TbP,P,T,number_of_elements,number_of_local_basisP,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,101,0,0,102,0,1);

A23 = a7;

A31 = A13.';

A32 = A23.';

sizeP = size(PbP,2);

A33 = zeros(sizeP,sizeP);

A = [A11 A12 A13;A21 A22 A23;A31 A32 A33];

sizeU = size(PbU,2);

b1 = assemble_vector_2D_test_v2('function_rho','function_fx',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,0,0);

b2 = assemble_vector_2D_test_v2('function_rho','function_fy',PbU,TbU,P,T,number_of_elements,number_of_local_basisU,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,102,0,0);

b3 = zeros(sizeP,1);

b = [b1;b2;b3];

%treat Dirichlet

[A,b] = treat_Dirichlet_2d_v2(A,b,boundarynodesU,boundarynodesP,'function_gU','function_gV','function_gP',PbU,PbP);

solution = A\b;

for i = 1:sizeU
    
    U(i) = solution(i);
    V(i) = solution(i+sizeU);
    
end

for i = 1:sizeP
    
    P(i) = solution(i+sizeU+sizeU);
    
end

hhx = h(1,1);
hhy = h(1,2);
hh2x = hhx/2;
hh2y = hhy/2;
N1 = 1/hh2x + 1;
N2 = 1/hh2y + 1;
N11 = 1/hhx + 1;
N22 = 1/hhy + 1;

x = 0:hh2x:1;
y = 1:-hh2y:0;
x1 = 0:hhx:1;
y1 = 1:-hhy:0;
[X,Y] = meshgrid(x,y);
[X1,Y1] = meshgrid(x1,y1);
Z1 = sparse(N1,N2);
Z2 = sparse(N1,N2);
Z3 = sparse(N11,N22);

k = 1;
for j = 1:1:N1
    for i = N2:-1:1
        Z1(i,j) = U(k);
        Z2(i,j) = V(k);
        k = k+1;
    end
end 

k = 1;
for j = 1:1:N11
    for i = N22:-1:1
        Z3(i,j) = P(k);
        k = k+1;
    end
end
%disp(Z);
%disp(U);
%figure
contourf(X,Y,Z1,'ShowText','on')
title('x-Velocity (u) Profile for h = 1/32')

contourf(X,Y,Z2,'ShowText','on')
title('y-Velocity (v) Profile for h = 1/32')

plot(Y(:,(N2+1)/2),Z1(:,(N2+1)/2));
title('Centerline x-Velocity Profile for h = 1/32')
ylabel('u at x/L = 0.5')
xlabel('y')
%contourf(X,Y,Z)
surf(X1,Y1,Z3)
title('Pressure Profile for h = 1/32')

error_on_nodes = generate_error_on_nodes_2d('function_exact', PbU, U);

infinity_norm_error=FE_solution_error_infinity_norm_triangle(U,'function_exact',left,right,bottom,top,h,102,0,0,Gauss_point_number);

L2_error=FE_solution_error_triangle(U,'function_exact',left,right,bottom,top,h,102,0,0,Gauss_point_number);

H1_error_x=FE_solution_error_triangle(U,'exact_solution_x_derivative',left,right,bottom,top,h,102,1,0,Gauss_point_number);

H1_error_y=FE_solution_error_triangle(U,'exact_solution_y_derivative',left,right,bottom,top,h,102,0,1,Gauss_point_number);

H1_error=sqrt(H1_error_x^2+H1_error_y^2);

%solution = error_on_nodes;
