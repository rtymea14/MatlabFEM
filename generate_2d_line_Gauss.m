function [Gauss_nodes_2d_line,Gauss_weights_2d_line] = generate_2d_line_Gauss(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,end_point_1,end_point_2)

if end_point_1(1) == end_point_2(1) 
    
    lower_bound=min(end_point_1(2),end_point_2(2));
    upper_bound=max(end_point_1(2),end_point_2(2));
    
    [Gauss_coefficient_local_1D,Gauss_point_local_1D_y]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    
    Gpn = size(Gauss_point_local_1D_y,2);
    
    Gauss_point_x = ones(1,Gpn);
    
    Gauss_point_x = end_point_1(1)*Gauss_point_x;
    
    Gauss_nodes_2d_line = [Gauss_point_x; Gauss_point_local_1D_y];
    
    Gauss_weights_2d_line = Gauss_coefficient_local_1D;  
    
elseif end_point_1(2) == end_point_2(2)
    
    lower_bound=min(end_point_1(1),end_point_2(1));
    upper_bound=max(end_point_1(1),end_point_2(1));
    
    [Gauss_coefficient_local_1D,Gauss_point_local_1D_x]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    
    Gpn = size(Gauss_point_local_1D_x,2);
    
    Gauss_point_y = ones(1,Gpn);
    
    Gauss_point_y = end_point_1(2)*Gauss_point_y;
       
    Gauss_nodes_2d_line = [Gauss_point_local_1D_x; Gauss_point_y];
    
    Gauss_weights_2d_line = Gauss_coefficient_local_1D;  
    
else
    
    lower_bound=min(end_point_1(1),end_point_2(1));
    upper_bound=max(end_point_1(1),end_point_2(1));
    
    [Gauss_coefficient_local_1D,Gauss_point_local_1D_x]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    
    Gpn = size(Gauss_point_local_1D_x,2);
    
    for i = 1:Gpn
    
        Gauss_point_y(i) = a * Gauss_point_local_1D_x(i) + b; 
    
    end
       
    Gauss_nodes_2d_line = [Gauss_point_local_1D_x; Gauss_point_y];
    
    Gauss_weights_2d_line = Gauss_coefficient_local_1D;  
    
end