% n: frame length, 128 by default
% m: frame shift, 64 by default
function ysm = iso7731(x, n, m, fs)


% pre-emphsis to flatten the spectrum
% PolyB = [1 -0.99];
% PolyA = 1;
% x = filter(PolyB, PolyA, x); 


% cut integer number of sliding window from the time series
x = x(1:floor((length(x) - n) / m) * m + n);
x = buffer(x, n, n-m, 'nodelay');





% prepare filtering bank
bn = 23;
ba = zeros(bn,2);

%ba(:,1) = 187.5:375:8000-375;
%ba(:,2) = ba(:,1) + 375;

ba(1,:) = [430 1000]; 
ba(2,:) = [500 1400]; 
ba(3,:) = [500 1700]; 
ba(4,:) = [600 1300]; 
ba(5,:) = [700 1000];

ba(6,:) = [700 1600]; 
ba(7,:) = [700 1700]; 
ba(8,:) = [900 1200]; 
ba(9,:) = [1000 1400];
ba(10,:)= [1100 1500];

ba(11,:) = [1300 2700]; 
ba(12,:) = [1500 1900]; 
ba(13,:) = [1600 2000]; 
ba(14,:) = [1600 3200];
ba(15,:) = [2050 2450];
ba(16,:) = [2100 2700];

ba(17,:) = [1875 2400]; 
ba(18,:) = [2400 2800]; 
ba(19,:) = [2800 3200];
ba(20,:) = [3200 3600];
ba(21,:) = [3500 3700]; 
ba(22,:) = [3600 4000];
ba(23,:) = [4000 4400];



ba = ba./fs.*n;
ba(:,1) = floor(ba(:,1));
ba(:,2) = ceil(ba(:,2));


% extract the features
y = zeros(bn, size(x,2));
for k = 1:size(x,2)
    
    % zero centering
    u = fft(x(:,k));
    mu = sum(x(:,k));
    u(1) = u(1) - mu;
    
    % band selection
    v = zeros(size(u));
    for bi = 1:bn
        % v must be complex conjugate
        v(ba(bi,1):ba(bi,2)) = u(ba(bi,1):ba(bi,2));
        v(n+2-ba(bi,2):n+2-ba(bi,1)) = u(n+2-ba(bi,2):n+2-ba(bi,1));
        
        % circular convolution
        rvv = real(fft(v.*conj(v)));
        [~, vmax] = quadmaxlocate(rvv', 0);
        if numel(vmax) <= 1
            y(bi, k) = 1;
        else
            vd = diff(vmax);
            y(bi, k) = std(vd)/mean(vd);  
        end
    end
end
y = prod(y).^(1/bn);

tc = 0.016;
alpha = 1-exp(-m/fs/tc);
ysm = zeros(size(y));
for i = 2:length(y)
    ysm(i) = y(i-1) * alpha + (1 - alpha) * ysm(i-1);
end
ysm = ysm';

%end of function
end