% covariance estimation using the Louis-formula
% rho=0, (r1,r2) corresponds to gamma frailty (r1, r2):
% x^{r1-1}exp\{-x/r2\}
% rho>0 corresponds to
function cov=Covest(coef, lambda, Delta,Z, indY, n,rho, r)

ncoef=length(coef); 

templambda1=lambda;
templambda1(Delta==0)=1;

Lambda=indY*lambda;
EZ1b=exp(Z*coef);

LambdabetaZ1=Lambda.*EZ1b;

[~, dG1, ddG1, d3G1]=Gtransform(LambdabetaZ1, rho, r);


s1=(Delta.*ddG1./dG1-dG1).*EZ1b;
t1=(Delta.*(d3G1./dG1-ddG1.^2./dG1.^2)-ddG1).*power(EZ1b,2);

Icoef=(1:1:ncoef)'; I1=ncoef+(1:1:n)'; 
ddL=zeros(ncoef+n, ncoef+n);
ddL(Icoef, Icoef)=Z'*(repmat(t1.*Lambda.^2+s1.*Lambda,1,ncoef).*Z);
ddL(Icoef, I1)=Z'*(repmat(t1.*Lambda,1,n).*indY)+Z'*(repmat(s1,1,n).*indY);
ddL(I1,I1)=indY'*(repmat(t1,1,n).*indY);
ddL(:, Icoef)=ddL(Icoef,:)';
ddL=ddL-diag([zeros(ncoef,1); Delta./templambda1.^2]);
information=-ddL;
index=[1+zeros(ncoef,1); Delta];
information=information(index==1, index==1);
cov=inv(information);



