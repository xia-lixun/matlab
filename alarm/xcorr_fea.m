function fea = xcorr_fea(y, nframe, nshift)

nsps = floor((length(y) - nframe) / nshift) * nshift + nframe;
y = y(1:nsps);
y = buffer(y, nframe, nframe-nshift, 'nodelay');

fea = zeros(1, size(y,2)-1);
for p = 1:size(y,2)-1
    u = y(:,p) - mean(y(:,p));
    v = y(:,p+1) - mean(y(:,p+1));  
    [rxx, lag] = xcorr(u,v,'coeff');
    [ymax, xmax] = quadmaxlocate(rxx(nframe+nshift:end)', min(rxx(nframe+nshift:end)));
    fea(p) = std(diff(xmax))/mean(diff(xmax));
end

% median filter
%fea = fea(1:floor(length(fea)/m)*m);
%fea = reshape(fea, m, length(fea)/m);
%fea = median(fea);

end