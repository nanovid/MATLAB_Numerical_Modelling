function [u] = u4a_bndryTENT(t)
    u=x;
    if(x<=1)
        u=x;
    elseif(x>1)
        u = (-1/3).*x +4/3
end