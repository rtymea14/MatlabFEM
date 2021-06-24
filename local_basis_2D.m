function result = local_basis_2D(x,y,vertices_triangle,basis_type,r,s,basis_index)

%using reference functions
%derivative w. r. t. x: r = some value && s = 0
%derivative w. r. t. y: r = 0 && s = some value
%derivative mix: r = some value && s = some value


x1 = vertices_triangle(1,1);
y1 = vertices_triangle(2,1);

x2 = vertices_triangle(1,2);
y2 = vertices_triangle(2,2);

x3 = vertices_triangle(1,3);
y3 = vertices_triangle(2,3);

J = [(x2-x1) (x3-x1);(y2-y1) (y3-y1)];

%Jacobian = det(J);
Jacobian=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

xhat = ((y3-y1)*(x-x1) - (x3-x1)*(y-y1))/Jacobian;

yhat = (-(y2-y1)*(x-x1)+(x2-x1)*(y-y1))/Jacobian;


if basis_type == 101
   
    if r == 0 && s == 0
        
        if basis_index == 1
            
            result = -xhat - yhat + 1;
                       
        elseif basis_index == 2
            
            result = xhat ;
                                    
        elseif basis_index == 3
            
            result = yhat;
            
        else
            
             warning = 'I am wrong!'
            
        end
       
    elseif r == 1 || s == 1
        
        if basis_index == 1
            
            dpsihat_by_dxhat1 = -1;
            dpsihat_by_dyhat1 = -1;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat1*(y3-y1)/Jacobian + dpsihat_by_dyhat1*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat1*(x1-x3)/Jacobian + dpsihat_by_dyhat1*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 2
            
            dpsihat_by_dxhat2 = 1;
            dpsihat_by_dyhat2 = 0;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat2*(y3-y1)/Jacobian + dpsihat_by_dyhat2*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat2*(x1-x3)/Jacobian + dpsihat_by_dyhat2*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 3
            
            dpsihat_by_dxhat3 = 0;
            dpsihat_by_dyhat3 = 1;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat3*(y3-y1)/Jacobian + dpsihat_by_dyhat3*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat3*(x1-x3)/Jacobian + dpsihat_by_dyhat3*(x2-x1)/Jacobian;
            end
            
        else
            
             warning = 'I am wrong!'
            
        end
        
    elseif r >= 2 || s >= 2
        
        result = 0;
        
    else
        
        warning = 'I am wrong!'
        
    end
    
elseif basis_type == 102
    
    if r == 0 && s == 0
        
        if basis_index == 1
            
            result = 2*xhat^2 + 2*yhat^2 + 4*xhat*yhat - 3*yhat - 3*xhat + 1;
            
        elseif basis_index == 2
            
            result = 2*xhat^2 - xhat;
            
        elseif basis_index == 3
            
            result = 2*yhat^2 - yhat;
            
        elseif basis_index == 4
            
            result = -4*xhat^2 - 4*xhat*yhat + 4*xhat;
            
        elseif basis_index == 5
            
            result = 4*xhat*yhat;
            
        elseif basis_index == 6
            
            result = -4*yhat^2 - 4*xhat*yhat + 4*yhat;
            
        else
            
             warning = 'I am wrong!'
            
        end
       
    elseif r == 1 || s == 1
        
        if basis_index == 1
            
            dpsihat_by_dxhat1 = 4*xhat + 4*yhat - 3;
            dpsihat_by_dyhat1 = 4*yhat + 4*xhat - 3;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat1*(y3-y1)/Jacobian + dpsihat_by_dyhat1*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat1*(x1-x3)/Jacobian + dpsihat_by_dyhat1*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 2
            
            dpsihat_by_dxhat2 = 4*xhat - 1;
            dpsihat_by_dyhat2 = 0;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat2*(y3-y1)/Jacobian + dpsihat_by_dyhat2*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat2*(x1-x3)/Jacobian + dpsihat_by_dyhat2*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 3
            
            dpsihat_by_dxhat3 = 0;
            dpsihat_by_dyhat3 = 4*yhat - 1;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat3*(y3-y1)/Jacobian + dpsihat_by_dyhat3*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat3*(x1-x3)/Jacobian + dpsihat_by_dyhat3*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 4
            
            dpsihat_by_dxhat4 = -8*xhat - 4*yhat + 4;
            dpsihat_by_dyhat4 = -4*xhat;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat4*(y3-y1)/Jacobian + dpsihat_by_dyhat4*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat4*(x1-x3)/Jacobian + dpsihat_by_dyhat4*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 5
            
            dpsihat_by_dxhat5 = 4*yhat;
            dpsihat_by_dyhat5 = 4*xhat;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat5*(y3-y1)/Jacobian + dpsihat_by_dyhat5*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat5*(x1-x3)/Jacobian + dpsihat_by_dyhat5*(x2-x1)/Jacobian;
            end
            
        elseif basis_index == 6
            
            dpsihat_by_dxhat6 = -4*yhat;
            dpsihat_by_dyhat6 = -8*yhat - 4*xhat + 4;
            
            if r == 1 && s == 0
                result = dpsihat_by_dxhat6*(y3-y1)/Jacobian + dpsihat_by_dyhat6*(y1-y2)/Jacobian;
            elseif r == 0 && s == 1
                result = dpsihat_by_dxhat6*(x1-x3)/Jacobian + dpsihat_by_dyhat6*(x2-x1)/Jacobian;
            end
            
        else
            
             warning = 'I am wrong!'
            
        end
        
     elseif r == 2 || s == 2 || (r+s) == 2
        
        if basis_index == 1
            
            d_dpsihat_by_dxhat1 = 4;
            d_dpsihat_by_dyhat1 = 4;
            d_dpsihat_by_dxhat_dyhat1 = 4;
            
            if r == 2 && s == 0
                result = d_dpsihat_by_dxhat1*((y3-y1)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat1*((y3-y1)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dyhat1*((y1-y2)^2)/(Jacobian^2);
            elseif r == 0 && s == 2
                result = d_dpsihat_by_dxhat1*((x1-x3)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat1*((x1-x3)*(x2-x1))/(Jacobian^2) + d_dpsihat_by_dyhat1*((x2-x1)^2)/(Jacobian^2);
            elseif r == 1 && s == 1
                result = d_dpsihat_by_dxhat1*((x1-x3)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat1*((x1-x3)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat1*((x2-x1)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dyhat1*((x2-x1)*(y1-y2))/(Jacobian^2);
            end
            
        elseif basis_index == 2
            
            d_dpsihat_by_dxhat2 = 4;
            d_dpsihat_by_dyhat2 = 0;
            d_dpsihat_by_dxhat_dyhat2 = 0;
            
            if r == 2 && s == 0
                result = d_dpsihat_by_dxhat2*((y3-y1)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat2*((y3-y1)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dyhat2*((y1-y2)^2)/(Jacobian^2);
            elseif r == 0 && s == 2
                result = d_dpsihat_by_dxhat2*((x1-x3)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat2*((x1-x3)*(x2-x1))/(Jacobian^2) + d_dpsihat_by_dyhat2*((x2-x1)^2)/(Jacobian^2);
           elseif r == 1 && s == 1
                result = d_dpsihat_by_dxhat2*((x1-x3)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat2*((x1-x3)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat2*((x2-x1)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dyhat2*((x2-x1)*(y1-y2))/(Jacobian^2);
            end
            
        elseif basis_index == 3
            
            d_dpsihat_by_dxhat3 = 0;
            d_dpsihat_by_dyhat3 = 4;
            d_dpsihat_by_dxhat_dyhat3 = 0;
            
            if r == 2 && s == 0
                result = d_dpsihat_by_dxhat3*((y3-y1)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat3*((y3-y1)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dyhat3*((y1-y2)^2)/(Jacobian^2);
            elseif r == 0 && s == 2
                result = d_dpsihat_by_dxhat3*((x1-x3)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat3*((x1-x3)*(x2-x1))/(Jacobian^2) + d_dpsihat_by_dyhat3*((x2-x1)^2)/(Jacobian^2);
            elseif r == 1 && s == 1
                result = d_dpsihat_by_dxhat3*((x1-x3)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat3*((x1-x3)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat3*((x2-x1)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dyhat3*((x2-x1)*(y1-y2))/(Jacobian^2);
            end
            
        elseif basis_index == 4
            
            d_dpsihat_by_dxhat4 = -8;
            d_dpsihat_by_dyhat4 = 0;
            d_dpsihat_by_dxhat_dyhat4 = -4;
            
            if r == 2 && s == 0
                result = d_dpsihat_by_dxhat4*((y3-y1)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat4*((y3-y1)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dyhat4*((y1-y2)^2)/(Jacobian^2);
            elseif r == 0 && s == 2
                result = d_dpsihat_by_dxhat4*((x1-x3)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat4*((x1-x3)*(x2-x1))/(Jacobian^2) + d_dpsihat_by_dyhat4*((x2-x1)^2)/(Jacobian^2);
           elseif r == 1 && s == 1
                result = d_dpsihat_by_dxhat4*((x1-x3)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat4*((x1-x3)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat4*((x2-x1)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dyhat4*((x2-x1)*(y1-y2))/(Jacobian^2);
            end
            
        elseif basis_index == 5
            
            d_dpsihat_by_dxhat5 = 0;
            d_dpsihat_by_dyhat5 = 0;
            d_dpsihat_by_dxhat_dyhat5 = 4;
            
            if r == 2 && s == 0
                result = d_dpsihat_by_dxhat5*((y3-y1)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat5*((y3-y1)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dyhat5*((y1-y2)^2)/(Jacobian^2);
            elseif r == 0 && s == 2
                result = d_dpsihat_by_dxhat5*((x1-x3)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat5*((x1-x3)*(x2-x1))/(Jacobian^2) + d_dpsihat_by_dyhat5*((x2-x1)^2)/(Jacobian^2);
            elseif r == 1 && s == 1
                result = d_dpsihat_by_dxhat5*((x1-x3)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat5*((x1-x3)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat5*((x2-x1)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dyhat5*((x2-x1)*(y1-y2))/(Jacobian^2);
            end
            
        elseif basis_index == 6
            
            d_dpsihat_by_dxhat6 = 0;
            d_dpsihat_by_dyhat6 = -8;
            d_dpsihat_by_dxhat_dyhat6 = -4;
            
            if r == 2 && s == 0
                result = d_dpsihat_by_dxhat6*((y3-y1)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat6*((y3-y1)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dyhat6*((y1-y2)^2)/(Jacobian^2);
            elseif r == 0 && s == 2
                result = d_dpsihat_by_dxhat6*((x1-x3)^2)/(Jacobian^2) + 2*d_dpsihat_by_dxhat_dyhat6*((x1-x3)*(x2-x1))/(Jacobian^2) + d_dpsihat_by_dyhat6*((x2-x1)^2)/(Jacobian^2);
            elseif r == 1 && s == 1
                result = d_dpsihat_by_dxhat6*((x1-x3)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat6*((x1-x3)*(y1-y2))/(Jacobian^2) + d_dpsihat_by_dxhat_dyhat6*((x2-x1)*(y3-y1))/(Jacobian^2) + d_dpsihat_by_dyhat6*((x2-x1)*(y1-y2))/(Jacobian^2);
            end
            
        else
            
             warning = 'I am wrong!'
            
        end
        
    elseif r >= 3 || r >= 3 || (r+s) >= 3
        
        result = 0;
        
    else
        
        warning = 'I am wrong!'
        
    end
    
end
    