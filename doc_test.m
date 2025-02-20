% This code runs the k1b document data experiment from Section 5.2 in the
% paper Junjun Pan and Michael Ng, Co-Separable Nonnegative
% Matrix Factorization, 2021.

clc
clear

maxiter=1000;  delta=1e-2;
it=10;

addpath('./projection-FGNSR');

load sub_k1b
load Datasub_k1b
data=sub;   
r=6; rr=6; rc=6;
 
MN=generalM_weight(M,IDXcol,IDXrow); 
  
  
%% original clustering ground truth

%%+++ for document +++
num_cluster=unique(classid);
m0=length(classid);n0=length(num_cluster);
mat_gt=zeros(m0,n0); 
for t=1:m0
mat_gt(t,classid(t))=1; % ground truth matrix
end
for i=1:n0
 xgt{i}=find(classid==i);%ground truth index
end

%%+++ for words +++
num_word_cluster=unique(word_id);
m1=length(word_id);n1=length(num_word_cluster);
mat_word_gt=zeros(m1,n1); 
for t=1:m1
mat_word_gt(t,word_id(t))=1; % ground truth matrix
end
for i=1:n1
 xgt_word{i}=find(word_id==i);%ground truth index
end


%% run CoS-FGM

   tic
   [C_fgm,R_fgm]=alfgm(MN,rc,rr,1e-3,maxiter, 100);
   toc
   
 % Compute Relative error of data:
 C_fgm=J_col(C_fgm);R_fgm=J_row(R_fgm);
  
tic
[Pr_fgm,S_fgm,Pc_fgm]=getweight_COSNMF(data,R_fgm,C_fgm,it,1e-3);
toc
 
disp('1-||M-MX-YM||_F / ||M||_F of data by FGNSR: %')
res_fgm=100-100*norm(data-Pr_fgm*S_fgm*Pc_fgm,'fro')/norm(data,'fro')

%Compute accuracy of clustering of data:

 
[mat_doc_fgm,mat_word_fgm]=clusteroutNMF(Pr_fgm,Pc_fgm);
[~,~,acc_doc_fgm] = Compare_clustering(mat_doc_fgm,mat_gt,2,0)
[~,~,acc_word_fgm] = Compare_clustering(mat_word_fgm,mat_word_gt,2,0)


%% run SPA+

%     tic
% % % %[C_spa,R_spa]=alspa(M,rc,rr,0.1,10)
%    [C_spa,R_spa] = bestrand_spa(MN,rc,rr);
%    toc

%  %+++ compute result of data +++
  
% Compute Relative error of data:
% tic
% C_spa=J_col(C_spa);R_spa=J_row(R_spa);
% [Pr_spa,S_spa,Pc_spa]=getweight_COSNMF(data,R_spa,C_spa,it,1e-3);
% toc 
%   
% disp('1-||M-XMY||_F / ||M||_F of data by SPA %')
% res_spa=100-100*norm(data-Pr_spa*S_spa*Pc_spa,'fro')/norm(data,'fro')
% 
% %Compute accuracy of clustering of data:
% 
%  [mat_doc_spa,mat_word_spa]=clusteroutNMF(Pr_spa,Pc_spa);
%  [~,~,acc_doc_spa] = Compare_clustering(mat_doc_spa,mat_gt,2,0)
%   [~,~,acc_word_spa] = Compare_clustering(mat_word_spa,mat_word_gt,2,0)
  

function [Wmat_out,Hmat_out]=clusteroutNMF(WL,WR)
  W= WL; H= WR'; %W: m*r; H: n*r
% 
% %normalize W: the maximum of each row of W  and H is set to 1 and other entries
% % are divided by the maximum. 
% 
% 
%  W=diag(1./max(W'))*W;
%  H=diag(1./max(H'))*H;
% 

overlap=0;
Wmat_out = Binarizing(W,overlap);
Hmat_out = Binarizing(H,overlap);
 end
 
 function H = Binarizing(H,overlap)
    
    % Get the number of rows of H
    [n,~] = size(H);
    
    if overlap==0
        
        % Iterate over the rows of H
        for i=1:n
            
            % Get the position of the maximum(s) of the ith row
            maxpos = find(H(i,:)==max(H(i,:)));
            
            % Take randomly one maximum and put it to 1
            newpos = randi(length(maxpos));
            
            % Initialize the row to zero
            H(i,:) = 0;
            
            % Put the randomly chosed position to 1
            H(i,maxpos(newpos)) = 1;
        end
    else
        % Not implemented yet
    end

end
 