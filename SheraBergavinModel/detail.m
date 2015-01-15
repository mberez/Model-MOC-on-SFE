function y=detail(x,n)
% y=detail(x,?n=1?)
% Returns all but the last n elements of x
% If x is a matrix, details each column.
if (nargin==1)
  n=1;
end
[M,N] = size(x);
if (M>1 & N>1)
  n = min ([abs(n),M]);
  y=x(1:length(x)-n,:);
else
  n = min ([abs(n),length(x)]);
  y=x(1:length(x)-n);
end

