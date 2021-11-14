% Estep
% rho=0, (r1,r2) corresponds to gamma frailty (r1, r2):
% x^{r1-1}exp\{-x/r2\}
% rho>0 corresponds to 
function Exi=Estep(coef, lambda, Delta,Z,indY, rho, r)
Lambda=indY*lambda;

LambdabetaZ=Lambda.*exp(Z*coef);

[~, dG1, ddG1, ~]=Gtransform(LambdabetaZ, rho, r);

Exi=-Delta.*ddG1./dG1+dG1;

     