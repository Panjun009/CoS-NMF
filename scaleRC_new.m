% Scale columns and rows of an m-by-n matrix M so that the sum of the
% entries in each row is n, and in each column is m
% 
% The algorithm alternatively scales the columns and rows of M. 
% It converges under some appropriate conditions; see for example
%
% Knight, P. A. (2008). The Sinkhorn–Knopp algorithm: convergence and 
% applications. SIAM Journal on Matrix Analysis and Applications, 30(1), 
% 261-275.
% and 
% Olshen, R. A., & Rajaratnam, B. (2010). Successive normalization of 
% rectangular arrays. Annals of statistics, 38(3), 1638.
% 
% *** Input ***
% M = 
% maxiter = maximum numner of iterations 
%
% *** Output ***
% Ms   = scaled matrix
% ec   = evolution of the error for the columns: sum_j | ||M(:,j)||_1-n |
% er   = evolution of the error for the rows: sum_i | ||M(i,:)||_1-m |
% conv = 1 if the algorithm has converged, otherwise the matrix might not
% be scalable, or the number of iterations was insufficient. 

function [Ms,CWeight,RWeight,ec,er,conv] = scaleRC_new(M,maxiter) 

if nargin <= 1
    maxiter = 100;
end
[m,n] = size(M); 
Ms = M; 
i = 1; 
conv = 1; 
if min(Ms(:)) < 0 || min( sum(Ms) ) == 0 || min( sum(Ms') ) == 0
    error('The input matrix should be nonnegative and have non-zero rows and columns.');
end
CWeight=eye(n);RWeight=eye(m);
while (i == 1 || er(i-1) > 1e-4&& ec(i-1) > 1e-4) && i <= maxiter 
    % Scale columns 
    CWeight=CWeight*diag(m./sum(Ms));
    Ms = Ms./repmat( sum(Ms), m ,1) *m; 
  
    er(i) = sum( abs( sum(Ms') - n ) ); 
    
    % Scale rows 
    RWeight=RWeight*diag(n./sum(Ms'));
    Ms = Ms./repmat( sum(Ms')', 1 , n) *n; 
  
    ec(i) = sum( abs( sum(Ms) - m ) ); 
    i = i+1; 
end  
if i == maxiter+1
    conv = 0; 
    warning('The matrix might not be scalable, or you might need to increase the number of iterations.');
end