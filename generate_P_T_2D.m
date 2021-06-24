function P = generate_P_T_2D(left,right,bottom,top,h,basis_type)

N1 = (right-left)/h(1);
N2 = (top-bottom)/h(2);

Nm = (N1+1)*(N2+1);

for i = 1:Nm
    
    r = mod(i,(N2+1));
    if r ~= 0
       column = fix(i/(N2+1)) + 1;
       row    = r;
    else 
        column = fix(i/(N2+1));
        row    = N2+1;
    end
    P (1,i) = left + (column-1)*h(1);
    P (2,i) = bottom + (row-1)*h(2);
    
end