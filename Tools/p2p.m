function [r,c] = p2p(vFrame,A,r,c,Bbor)
B = vFrame( r-Bbor:r+Bbor, c-Bbor:c+Bbor );
C = xcorr2n(A(:,:,iP),B);
[Mcor(iFrame,iP),I] = max(C(:));
%                     if Mcor(iFrame,iP) < 0.95 && (Bbor-Abor) < 9
%                         Bbor = Bbor + 2;
%                         continue
%                     end
[I_row, I_col] = ind2sub(size(C),I);
r = r + I_row - (Bbor-Abor + 1);
c = c + I_col - (Bbor-Abor + 1);