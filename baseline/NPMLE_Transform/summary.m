Ns = [1 2 4 8]*1000;
nsim=1000;
p=3;
for n = Ns
    load(strcat('../../res/baseline/NPMLE_output_N', num2str(n),'_m', num2str(nsim), '.mat'),...
        'output','runtime','imse');
    disp(n)
    disp(output(:,1:p))
    disp(mean(runtime))
end