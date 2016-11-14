function [ymax, xmax] = maxlocate(v,th)
% [X,Y] = quadmaxloc(V,T)  Quadratic-interpolated index and values for loc maxs
%      V is a uniformly-sampled vector.  Find all local maxes above T 
%      (absolute), then do a quadratic fit to interpolate the location and 
%      height of the maxima.   Return these as correspoding elements of X
%      and Y.
% 1998may02 dpwe@icsi.berkeley.edu $Header: $ again?

if (nargin < 2) 
  th = 0;
end

if (size(v,1) > 1) 
  error('v must be a row vector');
end

nr = size(v,2);

% catch for error case
if max(v) == -Inf
  % whole vector is -Inf, i.e. frame was digital zero
  xmax = [];
  ymax = [];
else
  % filter for local maxima; ensure edges don't win
  gtl = (v > [v(1), v(1:(nr-1))]);
  % allow greater-than-or-equal to catch plateaux
  gtu = (v >= [v(2:nr), 1+v(nr)]);
  vmax = v .* (v >= th) .* gtl .* gtu;
  maxixs = find(vmax);
  
  % Interpolate the max pos's
  xmax = maxixs;
  ymax = v(maxixs);

end