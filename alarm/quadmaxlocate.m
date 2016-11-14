function [ymax, xmax] = quadmaxlocate(v,th)
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
  xmax = zeros(size(maxixs));
  ymax = zeros(size(maxixs));
  for i = 1:size(maxixs,2)
    % Solve quadratic fit to 3 pts (as y = ax(x-b) with 0,0 as col(rmax-1))
    rmax = maxixs(i);    
    y1 = v(rmax)-v(rmax-1);
    y2 = v(rmax+1)-v(rmax-1);
    a = (y2 - 2*y1)/2;
    b = 1-y1/a;
    xmax(i) = rmax-1+b/2;
    ymax(i) = v(rmax-1)-a*b*b/4;
  end
end