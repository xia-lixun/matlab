function r = compZCR(x)
% compute the zero crossing rate

[L,M] = size(x);
y = [zeros(1,M);x(1:end-1,:)];

r = sum((x.*y)<0)/L;