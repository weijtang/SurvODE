function [G, dG, ddG, d3G]=Gtransform(temp, rho,r)


if (rho>0)
    G=(power(1+temp, rho)-1)/rho;
    dG=power(1+temp, rho-1);
    ddG=(rho-1)*power(1+temp, rho-2);
    d3G=(rho-1)*(rho-2)*power(1+temp, rho-3);
end

if (rho==0)
    newr=1/r;
    G=newr*log(1+temp/newr);
    dG=1./(1+temp/newr);
    ddG=-1/newr./power(1+temp/newr, 2);
    d3G=2/newr^2./power(1+temp/newr,3);
end


