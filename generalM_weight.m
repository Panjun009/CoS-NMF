function cMr_wgt=generalM_weight(M,IDXcol, IDXrow)
%% this function is to general a new data matrix by multiplying the weight to the columns and rows.
 CW=tabulate(IDXcol); 
 RW=tabulate(IDXrow);
 %The second column of CW/RW is the number of pixels connected to each cluster. Ref: help tabulate
 cMr_wgt=diag(sqrt(RW(:,2)))*M*diag(sqrt(CW(:,2)));
