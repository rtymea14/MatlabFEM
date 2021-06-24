function [error_on_nodes,infinity_norm_error,L2_error,H1_error] = Transient_2D_parabolic(left,right,bottom,top,h,basis_type,Gauss_point_number,bc,start_t,end_t,dt,theta)

%obtain necessary parameter

[P,T] = generate_P_T_2D_v2(left,right,bottom,top,h,101);
boundaryedges = generate_boundary_edge_info(left,right,bottom,top,h,101,T,bc);

if basis_type == 101
    
    Pb = P;
    Tb = T;
    boundarynodes = generate_boundary_nodes_info(left,right,bottom,top,h,basis_type,T,bc);
    
elseif basis_type == 102
    
  [Pb,Tb] = generate_P_T_2D_v2(left,right,bottom,top,h,basis_type);
  boundarynodes = generate_boundary_nodes_info(left,right,bottom,top,h,basis_type,T,bc);
  
else
end

number_of_elements = size(T,2);

if basis_type == 101
    number_of_local_basis = 3;
elseif basis_type == 102
    number_of_local_basis = 6;
end

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_point_number);

[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number);

A1 = assemble_matrix_2D_trial_test('function_c',Pb,Tb,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,1,0,basis_type,1,0);

A2 = assemble_matrix_2D_trial_test('function_c',Pb,Tb,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,0,1,basis_type,0,1);
   
A = A1 + A2;

%Assemble Mass matrix

Mass = assemble_matrix_2D_trial_test('function_one',Pb,Tb,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,0,0,basis_type,0,0);

A_tilde = Mass/dt + theta*A;
A_fixed = Mass/dt - (1-theta)*A;

%time iteration
M = (end_t - start_t)/dt;

X_old = generate_initial_vector('function_initial',Pb);

for m = 0:M-1
    
    tm = m*dt;
    tmp1 = (m+1)*dt;

    bm = assemble_time_vector_2D_test('function_f',tm,Pb,Tb,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,0,0);
    bmp1 = assemble_time_vector_2D_test('function_f',tmp1,Pb,Tb,P,T,number_of_elements,number_of_local_basis,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,basis_type,0,0);
    
    b_tilde = theta*bmp1 + (1-theta)*bm + A_fixed*X_old;

    %treat Neumann
    vector_size = size(b_tilde,1);

    v = treat_Neumann_2d('function_c','function_p',Gauss_coefficient_reference_1D,Gauss_point_reference_1D,boundaryedges,number_of_local_basis,P,T,Tb,vector_size,basis_type,0,0);

    b_tilde = b_tilde + v;

    %treat Robin

    matrix_size = size(A_tilde,1);

    [w,R] = treat_Robin_2d('function_c','function_r','function_q',Gauss_coefficient_reference_1D,Gauss_point_reference_1D,boundaryedges,number_of_local_basis,number_of_local_basis,P,T,Tb,vector_size,matrix_size,basis_type,basis_type,0,0,0,0);

    b_tilde = b_tilde + w;

    A_tilde = A_tilde + R;

    %treat Dirichlet

    [A,b] = treat_Dirichlet_2d_time(A_tilde,b_tilde,boundarynodes,'function_g',tmp1,Pb);
    
    %complete numerical solution

    solution = A\b;
    
    X_old = solution;
    
end

error_on_nodes = generate_error_on_nodes_2d('function_exact', end_t, Pb, solution);

infinity_norm_error=FE_solution_error_infinity_norm_triangle(solution,'function_exact',end_t,left,right,bottom,top,h,basis_type,0,0,Gauss_point_number);

L2_error=FE_solution_error_triangle(solution,'function_exact',end_t,left,right,bottom,top,h,basis_type,0,0,Gauss_point_number);

H1_error_x=FE_solution_error_triangle(solution,'exact_solution_x_derivative',end_t,left,right,bottom,top,h,basis_type,1,0,Gauss_point_number);

H1_error_y=FE_solution_error_triangle(solution,'exact_solution_y_derivative',end_t,left,right,bottom,top,h,basis_type,0,1,Gauss_point_number);

H1_error=sqrt(H1_error_x^2+H1_error_y^2);

%solution = error_on_nodes;