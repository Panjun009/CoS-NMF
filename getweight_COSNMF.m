% note: for artificial data set: nnls --1000times, iter=100;
% note: for real dataset : nnls---100 times; iter=10;
function [Pr,S,Pc,i]=getweight_COSNMF(M,R,C,iter, delta)

S=M(R,C);
 
[m,n]=size(M); [rr,rc]=size(S);

Pr_pre=zeros(m,rr); 
Pc_pre=zeros(rc,n);

   Mr=M(R,:); Mc=M(:,C); 
% % % %  % initial 1st
%   Pr_aft=nnlsHALSupdt_new(S*Mc',S',[],1000)'; 
%   Pc_aft=nnlsHALSupdt_new(S'*Mr,S,[],1000);
  % % % % % %     
% % % %   %  initial 2nd
Pr_aft = nnlsHALSupdt_new(Mr*M',Mr',[],100)';
Pc_aft = nnlsHALSupdt_new(Mc'*M,Mc,[],100);

e1= norm(Pr_aft-Pr_pre,'fro');%/norm(Pr_pre,'fro')
e2= norm(Pc_aft-Pc_pre,'fro');%/norm(Pc_pre,'fro')
e=e1+e2;

i=1;
while i< iter && e >= delta
    
    Pr_pre = Pr_aft;     Pc_pre = Pc_aft;
    
    W = Pr_aft*S;
    Pc_aft =  nnlsHALSupdt_new(W'*M,W,[],100); % solve from W*Pc_aft~=M;
    
    H = S*Pc_aft;
    Pr_aft = nnlsHALSupdt_new(H*M',H',[],100)'; % solve from Pr_aft*H~=M;
    
    e1= norm(Pr_aft-Pr_pre,'fro');%/norm(Pr_pre,'fro')
    e2= norm(Pc_aft-Pc_pre,'fro');%/norm(Pc_pre,'fro')
    e=e1+e2;
    
    i=i+1;

end
Pr = Pr_aft;
Pc = Pc_aft;
% 

% function [Pr,S,Pc]=getweight_COSNMF(M,R,C)
%    S=M(R,C);
%    
%    Mr=M(R,:); Mc=M(:,C); 
%   Pr = nnlsHALSupdt_new(S*Mc',S',[],100)'; 
%   Pc = nnlsHALSupdt_new(S'*Mr,S,[],100);
  
  
%    % 2nd
% Pr_aft = nnlsHALSupdt_new(Mr*M',Mr',[],100)';
% Pc_aft = nnlsHALSupdt_new(Mc'*M,Mc,[],100);

