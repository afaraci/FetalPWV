function video = createFrame(pwv)

%pwv = -3;
fps = 148;

numXpix = 1024; % 
numYpix =  512; %

P = generatePWVphantom( pwv, fps, 2000);
amp = abs( randn(size(P,1),1)*4);
P(:,2,:) = [];
Q(:,1,:) = ( P(:,1,:) + 0.02 )/0.04*( numXpix - 1 ) + 1;
Q(:,2,:) = ( P(:,2,:) - 0.05 )/0.02*( numYpix - 1 ) + 1;
for iFrame = 1 : size(Q,3)
    P = Q(:,:,iFrame);
    for i = 1 : size(P,1); H(i,:) = [ P(i,2)-amp(i)   P(i,2,1)+amp(i)   ]; end
    for i = 1 : size(P,1); W(i,:) = [ P(i,1)-amp(i)*2 P(i,1,1)+amp(i)*2 ]; end
    H = floor( H );
    W = floor( W );
    H( H>numYpix ) = numYpix;
    H( H<1) = 1;
    W( W>numXpix ) = numXpix;
    W( W<1) = 1;
    I = zeros( numYpix, numXpix );
    for i = 1 : size(P,1)
        vH = H(i,1) : H(i,2);
        vW = W(i,1) : W(i,2);
        I( vH, vW ) = 1;
    end
    %figure,imshow(I)
    
    PSF = fspecial( 'gaussian', [12 60], 20 );
    % I = imfilter( I, PSF );
    frames(:,:,iFrame) = I;
end

%%
frames = frames - min( frames(:) );
frames = frames / max( frames(:) );
% frames = frames*255;
% frames = round(frames);
%frames = uint8(frames);

%%
video = frames;
% cmppxx = ( xSpace(end) - xSpace(1) ) / numXpix / 10; % [cm]
% cmppxy = ( ySpace(end) - ySpace(1) ) /  numYpix / 10; % [cm]
% save( fullfile( dataDir, fileName ),'video','fps','cmppxx')
implay(video);