format LONG
left = 0; right = 2;bottom = 0; top = 1;
h = [1/4 1/4;1/8 1/8;1/16 1/16;1/32 1/32];
%h = [1/8 1/8;1/16 1/16;1/32 1/32;1/64 1/64];
basis_type = 101;
Gauss_point_number = 4;

%function f should be -10e^(x+y)

% Using Dirichlet Boundary Condition
%function g should be
%u = e^y when x = 0
%u = e^(2+y) when x = 2
%u = e^x when y = 0
%u = e^(x+1) when y = 1

for i = 1:1:4
[max_error(i),infinity_norm_error(i),L2_error(i),H1_error(i)] = Steady_2D_elliptic(left,right,bottom,top,h(i,:),basis_type,Gauss_point_number,-1);
end

disp('max error on nodes for linear element = ');
disp(max_error');
disp('inifinity norm error on nodes for linear element = ');
disp(infinity_norm_error');
disp('L2_error on nodes for linear element = ');
disp(L2_error');
disp('H1_error on nodes for linear element = ');
disp(H1_error');

basis_type = 102;

for i = 1:1:4
[max_error(i),infinity_norm_error(i),L2_error(i),H1_error(i)] = Steady_2D_elliptic(left,right,bottom,top,h(i,:),basis_type,Gauss_point_number,-1);
end

disp('max error on nodes for quadratic element = ');
disp(max_error');
disp('inifinity norm error on nodes for quadratic element = ');
disp(infinity_norm_error');
disp('L2_error on nodes for quadratic element = ');
disp(L2_error');
disp('H1_error on nodes for quadratic element = ');
disp(H1_error');