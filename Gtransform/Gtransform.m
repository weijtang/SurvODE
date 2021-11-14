function [G, dG]=Gtransform(temp, rho,r)
% specifying the G-transform function and its gradient given (rho1, r1) 

if (rho>0)
    G=(power(1+temp, rho)-1)/rho;
    dG=power(1+temp, rho-1);
end

if (rho==0)
    if r==0 %cox as a special case
        G = temp;
        dG = ones(size(temp));
    else
        newr=1/r;
        G=newr*log(1+temp/newr);
        dG=1./(1+temp/newr);
    end
end


