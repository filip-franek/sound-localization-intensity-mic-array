function [wout]=AdjustWeights(weights,level)
% This function defines user's weighting
% To apply this function within T-ABCD, the weighting type must be 
% set to 'Constrained', 'N.Constrained' or % 'F.Constrained' 
% (depending on the clustering method)

if nargin<2
   level=0.15;
end   

low=max(floor(level*length(weights)),1);
high=ceil((1-level)*length(weights));
len=length(weights);
wout=zeros(1,len);

[v,ind]=sort(weights,'ascend');
wout(ind(1:low))=0;
wout(ind(high:end))=1;
wout(ind(low:high))=weights(ind(low:high));