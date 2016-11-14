function fea = acorr_fea(y, nframe, nshift)

nsps = floor((length(y) - nframe) / nshift) * nshift + nframe;
y = y(1:nsps);
y = buffer(y, nframe, nframe-nshift, 'nodelay');

fea = zeros(1, size(y,2)-1);
for p = 1:size(y,2)
    u = y(:,p) - mean(y(:,p));
    %[rxx, ~] = xcorr(u,u,'coeff');
    rxx = cconv(u,u,length(u));
    [~, xmax] = quadmaxlocate(rxx', 0);
    fea(p) = std(diff(xmax))/mean(diff(xmax));
end

end