function I = numint(y,dx,dim)
% function I = numint(y, ?dx or x?,?dim?)
% Numerical integral. Points x are assumed equally spaced
% with spacing dx. A higher order version of trapz.
% Note that the order of the arguments is reversed
% with respect to cumtrapz!  See cumint for cumulative version.
  
  if (nargin<2 | isempty(dx))
    dx = 1;
  elseif (length(dx)>1)
    dx = dx(2)-dx(1);
  end
  if (nargin<3)
    dim = 1;
  elseif (dim>2)
    error ('dim>2 not implemented');
  end
  
  sz = size(y);
  npts = sz(dim);
  if (isvector(y))
    y = y(:);
    npts = length(y);
    dim = 1;
  end

  % assume for the moment that npts>=6
  if (npts<6)
    error ('npts<6 not yet implemented');
  end
  
  c6 = [3/8;7/6;23/24];			% Numerical Recipes, 4.1.14
  cn = ones([npts,1]);
  cn(1:3) = c6;
  cn((end-2):end) = reverse(c6);

  if (~isvector(y));
    if (dim==1)
      cn = repmat(cn,1,sz(2));
    else
      cn = cn';
      cn = repmat(cn,sz(1),1);
    end
  end
  
  I = dot(cn,y,dim)*dx;
  
  return




