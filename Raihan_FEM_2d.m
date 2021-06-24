left = 0; right = 1;bottom = 0; top = 1;
h = [1/4 1/4];

disp('P matrix generated using version 1');
P1 = generate_P_T_2D(left,right,bottom,top,h,101);
disp(P1);

disp('P matrix generated using version 2');
[P,T] = generate_P_T_2D_v2(left,right,bottom,top,h,101);
disp(P);
disp('T matrix');
disp(T);
disp('boundaryedges matrix');
boundaryedges = generate_boundary_edge_info(left,right,bottom,top,h,101,T);
disp(boundaryedges);
disp('boundarynodes for linear 2d elements');
boundarynodes1 = generate_boundary_nodes_info(left,right,bottom,top,h,101,T);
disp(boundarynodes1);

basis_type = 102;
if basis_type == 101
    
    Pb = P;
    Tb = T;
    
elseif basis_type == 102
    
  [Pb,Tb] = generate_P_T_2D_v2(left,right,bottom,top,h,basis_type);
  
else
end

disp('Pb matrix for quadratic 2d elements');
disp(Pb);
disp('Tb matrix for quadratic 2d elements');
disp(Tb);
disp('boundarynodes for quadratic 2d elements');
boundarynodes = generate_boundary_nodes_info(left,right,bottom,top,h,basis_type,T);
disp(boundarynodes);