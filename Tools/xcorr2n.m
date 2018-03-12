function C = xcorr2n(A,B)

A = double(A);
B = double(B);

dR = size(B,1) - size(A,1) + 1;
dC = size(B,2) - size(A,2) + 1;
v = A(:);
C = zeros(dR,dC);
for iR = 1 : dR
    ro = iR : size(A,1)+iR-1;
    for iC = 1 : dC
        co = iC : size(A,2)+iC-1;
        Bi = B(ro,co);
        w = Bi(:);
        C(iR,iC) = dot(v,w)/(norm(v)*norm(w));
    end
end