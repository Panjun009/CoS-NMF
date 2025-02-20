% This code runs the synthetic data experiment from Section 5.1 in the
% paper Junjun Pan and Michael K. Ng, Co-Separable Nonnegative
% Matrix Factorization, 2021. 

clear all; clc; close all; 
addpath(genpath('./utils/')); 

% noise levels
epsilon =  logspace(-7,-1,20); 
% In the paper: logspace(-7,-1,20) + average over 25 runs
% Synthetic data sets generated fully randomly.
 
delta=1e-6; it=100;  maxiter=1000; 

for i = 1 : length(epsilon) 
    fprintf('****** Run for epsilon = %2.2f *********\n', epsilon(i) ); 
    % Generate a Co-separable matrix M
   
        m = 100; 
        n = 100; 
        rr = 10; 
        rc = 3; 
        r=rr+rc;
        S=rand(rr,rc);W=rand(m-rr,rr);H= rand(rc,n-rc);
        M0=[S,S*H;W*S,W*S*H];
  
    
    % Add noise 
        Noise=randn(m,n);  
        Noise=epsilon(i)*(Noise/norm(Noise,'fro'))*norm(M0,'fro'); 
        M=max(M0+Noise,0);
    
    % Add random permutations
    Pc = eye(n);
    %permc = randperm(n);
    permc = 1:n; 
    Pc = Pc(:,permc);
    Pr = eye(m);
    %permr = randperm(m);
    permr = 1:m; 
    Pr = Pr(permr,:);
    M = Pr*M*Pc;
    
 
    [M,Cweight,Rweight] = scaleRC_new(M,1000);
    
    
    % Groundtruth
    ic = find(permc <= rc);
    ir = find(permr <= rr);
    
    % *** Run CoS-FGM ***
      [Jc_fgm,Jr_fgm]=alfgm(M,rc,rr,delta,maxiter,100);
      [Pr_fgm,S_fgm,Pc_fgm]=getweight_COSNMF(M,Jr_fgm,Jc_fgm,it,delta);


    % Compute accuracy
    accCoS(i) = (length(intersect(ic, Jc_fgm)) + length(intersect(ir, Jr_fgm))) / (rr+rc);
    fprintf('Selection accuracy by CoS-FGM is %2.2f\n', accCoS(i));
    % Compute Relative error:
    resCoS(i)=100-100*norm(M-Pr_fgm*S_fgm*Pc_fgm,'fro')/norm(M,'fro'); 
      
    
    
    % *** Run SPA+ to etxract rr rows and rc column *** 
    
     [Jc_spa,~] = FastSepNMF(M,rc,1); 
     [Jr_spa,~] =FastSepNMF(M',rr,1);
     [Pr_spa,S_spa,Pc_spa]=getweight_COSNMF(M,Jr_spa,Jc_spa,it,delta);

    
   % Compute accuracy
    accSPA(i) = (length(intersect(ic, Jc_spa)) + length(intersect(ir,Jr_spa))) / r;

    fprintf('Selection accuracy by SPA* is %2.2f\n', accSPA(i));
    % Compute Relative error:
    resSPA(i)=100-100*norm(M-Pr_spa*S_spa*Pc_spa,'fro')/norm(M,'fro');
end

% Display the accuracy
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 2);
set(0,'DefaultLineMarkerSize', 10);
figure; 
subplot(1,2,1); 
semilogx(epsilon, accCoS, 'rx--'); 
hold on; 
semilogx(epsilon, accSPA, 'ko--'); 
legend('CoS-FGM',  'SPA^+'); 
ylabel('Accuracy'); 
xlabel('Noise level (\epsilon)'); 
axis([1e-7 1e-1 0 1])
% Display the relative approximation quality 
subplot(1,2,2); 
loglog(epsilon, resCoS, 'rx--'); 
hold on; 
loglog(epsilon, resSPA, 'ko--');
legend('CoS-FGM', 'SPA*'); 
ylabel('Approximation error'); 
xlabel('Noise level (\epsilon)'); 
axis([1e-7 1e-1 0 100])
