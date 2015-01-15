function rx = reverse (x)
% Reverse the vector x. If x is a matrix, performs
% a column-wise reversal equivalent to flipud(x);

[m,n]=size(x);
if (m>1 & n>1)
  rx = flipud(x);
else
  rx = x(length(x):-1:1);
end
