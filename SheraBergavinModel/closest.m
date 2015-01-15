function [x_n,n] = closest (x,value)
% function [x_n,n] = closest (x,value)
% Return element x_n of x whose value is closest
% to the specified value. Value can be a 
% vector of elements, in which case x_n and
% n are also vectors.

  % wrapper for nearestpoint (which handles large arrays well)...
  n = nearestpoint(value,x);
  x_n = x(n);

  [rx,cx] = size(x);
  if (rx > 1)				% give n same shape as x_n
    n = n(:);
  end
  
  return
