format LONG
left = 0; right = 1;bottom = 0; top = 1;
%h = [1/32 1/32];
h = [1/4 1/4;1/8 1/8;1/16 1/16;1/32 1/32];

Gauss_point_number = 4;


for i = 1:1:4
[max_error(i),infinity_norm_error(i),L2_error(i),H1_error(i)] = Steady_2D_elliptic(left,right,bottom,top,h(i,:),Gauss_point_number,-1);
end

disp('max error on nodes for U = ');
disp(max_error');
disp('inifinity norm error on nodes for U = ');
disp(infinity_norm_error');
disp('L2_error on nodes for U = ');
disp(L2_error');
disp('H1_error on nodes for U = ');
disp(H1_error');