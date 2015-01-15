function [x,dx] = linscale (x0,x1,npts,rowvec);
% function [x,dx] = linscale (x0,x1,?npts=128?,?rowvec=true?)  
% Create linearly spaced vector with npts elements ranging from x0 to
% x1.  Similar to linspace, but with a different default value of npts
% and an orientation flag.  The value of rowvec can also be either
% 'row' or 'col', depending on taste.  Implemented for compatibility
% with logscale.
%
% Christopher A. Shera (shera at mit.edu)
% Eaton-Peabody Laboratory  
%

if (nargin < 3), npts = 128; end;
if (nargin < 4), rowvec = true; end;

if (ischar(rowvec))
  if (rowvec(1)=='r'), rowvec = true;
  elseif (rowvec(1)=='c'), rowvec = false; end
end

x = linspace (x0,x1,npts);
x(1) = x0; x(end) = x1;		% make sure end points match exactly
if (nargout == 2)
  dx = (x1-x0)/(npts-1);
end

if (~rowvec)
  x = x(:);
end
