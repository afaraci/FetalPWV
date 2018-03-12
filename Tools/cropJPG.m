function cropJPG( fileName, border, bDebug )

if nargin < 2, border = 0; end
if nargin < 3, bDebug = 0; end

I = imread(fileName);
r = mean(mean(I,3),2);
c = mean(mean(I,3),1);
r1 = find( ~(r==255), 1, 'first' );
re = find( ~(r==255), 1, 'last' );
c1 = find( ~(c==255), 1, 'first' );
ce = find( ~(c==255), 1, 'last' );

I = I(r1:re,c1:ce,:);
J = I;

perm1 = [3 2 1];
perm2 = [3 1 2];
perms = [perm1; perm2];
for iPerm = 1 : 2
    J = permute(J,perms(iPerm,:));
    for count = 1 : 2
        J(:,:,end+1:end+border) = ones(size(J,1),size(J,2),border)*255;
        J = flip(J,3);
    end
end
J = permute(J,[1 3 2]);
I = J;
imwrite(I,fileName);
if bDebug
    dos(fileName);
    pause
end