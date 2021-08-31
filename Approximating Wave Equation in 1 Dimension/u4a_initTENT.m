function [u] = u4a_initTENT(x)
u=x;
if(x<=1)    
    u = x;
 
elseif(x>1)
    u = (-1/3).*x;
u;
end 