function [l]=lagrangeb(knots,x,x_j)
l=ones(size(x));
ihat = find((x_j-knots)<10^(-10));
knotshat = [knots(1:ihat-1),knots(ihat+1:end)];
if ~isempty(knotshat)
    l=l.*prod(bsxfun(@minus,x,knotshat),2)./prod(x_j-knotshat); 
end