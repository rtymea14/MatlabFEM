function r = function_g(x,y)

if x == 0
    
    r = exp(y);
    
elseif x == 2
    
    r = exp(2+y);
    
elseif y == 0
    
    r = exp(x);
    
elseif y == 1
    
    r = exp(x+1);

end