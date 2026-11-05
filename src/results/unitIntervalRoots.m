function values = unitIntervalRoots(p, imaginaryTolerance)
%UNITINTERVALROOTS Real roots in [0,1], coefficients in ascending powers.
scale = max(abs(p));
values = [];
if scale == 0, return; end
p = p / scale;
while numel(p) > 1 && abs(p(end)) < 1e-14, p(end) = []; end
if numel(p) < 2, return; end
r = roots(fliplr(p));
values = real(r(abs(imag(r)) < imaginaryTolerance & real(r) >= 0 & real(r) <= 1)).';
end
