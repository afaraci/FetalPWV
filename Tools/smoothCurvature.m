function s = smoothCurvature( xvd, s, t, bDebug )

if nargin < 4
    bDebug = 0;
end

cur = abs(diff(s,2));
m = max(cur);
while m > t
    cur = abs(diff(s,2));
    [m,idx] = max(cur);
    idx = idx+1;
    s(idx) = mean(s(idx-1:idx+1));
end

if bDebug
    plot(xvd,s,'LineWidth',2)
end