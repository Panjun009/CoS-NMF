 %% ****** Alternating fast gradient method *******

   function [Jc_fgm,Jr_fgm,i]=alfgm(M,rc,rr,delta,maxiter,iter)
  [m,n] = size(M);
  
   MXP=M; MYP=M;
   
  [~, Jr_fgm] = fgnsr_new(MXP',rr,delta, maxiter);
  [~, Jc_fgm] = fgnsr_new(MYP,rc,delta, maxiter);
   MX = MXP(Jr_fgm,:);   MY = MYP(:,Jc_fgm); 
 
   WrP = zeros(rr,n);
   WcP = zeros(m,rc);
   
   ep1=norm(MX-WrP,'fro');%/norm(MX);
   ep2=norm(MY-WcP,'fro');%/norm(MY);
   i=1;
while i <= iter && ep1+ep2 >= delta %&& ep2 >= delta
      
      MXP = MX; 
      MYP = MY;
      
 [~, Jr_fgm] = fgnsr_new(MY',rr,delta, maxiter) ;
  MX = M(Jr_fgm,:); 
  
 [~, Jc_fgm] = fgnsr_new(MX,rc,delta, maxiter) ;
  MY = M(:,Jc_fgm); 
  
  ep1=norm(MX-MXP,'fro');%/norm(MX);
  ep2=norm(MY-MYP,'fro');%/norm(MY);
  i=i+1;
end
end
