% Mstep

function [newcoef, newlambda]=Mstep(oldcoef, Exi, Delta,Z, indY, n)

                                          
d=length(oldcoef);
EbetaZ=exp(Z*oldcoef);

score1=zeros(d,n);
dscore1=zeros(d,d);
for i=1:1:n
    temp=indY(:,i).*EbetaZ.*Exi;
    tempnom=Z'*temp;
    tempdenom=sum(temp);
    tempnom2=Z'*(repmat(temp,1,d).*Z);
    score1(:,i)=Delta(i)*(Z(i,:)'-tempnom/tempdenom);
    dscore1=dscore1-Delta(i)*(tempnom2/tempdenom-tempnom*tempnom'/tempdenom^2);
end

newcoef=oldcoef- dscore1\(score1*ones(n,1));
EnewbetaZ=exp(Z*newcoef);
newlambda=Delta./(indY'*(EnewbetaZ.*Exi));