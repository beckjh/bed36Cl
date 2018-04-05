% ===========================================================
%  knots and weights of the Clenshaw-Curtis quadrature formula
% ===========================================================
%
%   [x,w] = knots_CC(nn,x_a,x_b)
%   nn: number of knots
%   x_a,x_b: interval containing the nodes
%   x: knots   [row vector]
%   w: weights [row vector]

function [x] = knots_CC(nn,x_a,x_b)

if nn==1
    x=(x_a+x_b)/2;
elseif mod(nn,2)==0
    error('error in knots_CC: Clenshaw-Curtis formula \n use only odd number of points')
else
    n=nn-1;

    x=(cos([0:n]*pi/n));
    x = (x_b-x_a)/2*x + (x_a+x_b)/2;
    x = sort(x);
end
