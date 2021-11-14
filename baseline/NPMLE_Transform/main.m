function main(n, nsim)
% main program
setting = 2;
rho = 0; r = 1;
p=3;

npt=100; testpt=linspace(0.01, 1.5, npt).';
true_H0=2*testpt;
truepar=[1;1;1; true_H0];
parest=zeros(nsim, p+npt);
sdest=zeros(nsim,  p+npt);
runtime=zeros(nsim, 1);
imse=zeros(nsim,1);
cenrate=0; 
for simind=1:1:nsim
    % load data
    data_file = strcat('../../data/survode/simudata_N', num2str(n), '_seed', num2str(simind), ...
        '_setting', num2str(setting), '.mat');
    
    if isfile(data_file)
        load(data_file, 'x', 'time', 'delta');
    else
        generator(n, simind, setting, rho, r);
        load(data_file, 'x', 'time', 'delta');
    end
    
    indY=(repmat(time,1,n)>=repmat(time',n,1));
    n1=sum(delta);
    
    oldcoef=zeros(p, 1); 
    oldlambda=delta./(indY*ones(n,1));

    
    epsilon=0.001;maxiter=100;error=1;iter=0;
    tic
    while (error>epsilon && iter< maxiter)
        Exi=Estep(oldcoef, oldlambda, delta, x, indY, rho, r);
        [newcoef, newlambda]=Mstep(oldcoef,Exi, delta, x,indY, n);
        error=sum(abs(newcoef-oldcoef))+sum(abs(newlambda-oldlambda));
        iter=iter+1;
        oldcoef=newcoef;
        oldlambda=newlambda;
    end
    runtime(simind) = toc;
    
    cov=Covest(newcoef, newlambda, delta, x,indY, n, rho, r);
    predmatrix=[diag(1+zeros(p,1)), zeros(p, n1);
        zeros(npt, p), repmat(testpt, 1,n1)>=repmat(time(delta==1)', npt, 1)];

    tempest=predmatrix*[newcoef; newlambda(delta==1)];
    tempsd=sqrt(diag(predmatrix*cov*predmatrix'));
    parest(simind,:)=tempest';
    sdest(simind,:)=tempsd';
    cenrate=(cenrate*(simind-1)+n1/n)/simind;
    imse(simind) = sum(abs(true_H0 - tempest((p+1):end)).^2) * 1.5/99;
    
    if (mod(simind, 50)==0)
        Sind=(1:1:simind)';
        disp([simind, cenrate, iter])
        indLambda=p+(1:npt);
        output=[truepar'; mean(parest(Sind,:))-truepar'; 
            sqrt(var(parest(Sind,:)));
            mean(sdest(Sind,:));
            mean(abs(parest(Sind,:)-repmat(truepar', simind,1))<1.96*sdest(Sind,:))];
        output(5, indLambda)=...
            mean(abs(log(parest(Sind, indLambda))-repmat(log(truepar(indLambda)'), simind,1))...
                 <1.96*sdest(Sind,indLambda)./parest(Sind, indLambda));
        disp(output(:,1:p))
        disp(mean(runtime(Sind)))
        disp(mean(imse(Sind)))
        save(strcat('../../res/baseline/NPMLE_output_N', num2str(n),'_m', num2str(nsim), '.mat'),...
            'parest','sdest','output','runtime','imse');
    end
end

