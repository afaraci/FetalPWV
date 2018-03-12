function generatePWV(pwv)

% pwv in m/s

if nargin < 1
    pwv = 3;
end

close all
load simsig
h = 512;
w = 1024;
nFrames = 296;
fps = 148; % Frames per second
% lenAor = 30; % Length of Aorta in mm
lenAor = (w-1)*cmppxx*10; % Length of aorta [mm]
t = 1/fps : 1/fps : 2; % time
d = 1 : w;
d = d*cmppxx*10; d = d - d(1); % [mm]
for iD = 1 : w
    y(:,iD) = feval( fittedmodel, t + (lenAor/pwv/1000)/lenAor*d(iD) );
end
d = round( d/(cmppxx*10) ) + 1;
y = round( y/(cmppxx*10) );
yT = round(h/2 - y/2);
yB = round(h/2 + y/2);
video = ones(h,w,nFrames);
for iFrame = 1 : nFrames
    for iCul = 1 : w
        video(yT(iFrame,iCul):yB(iFrame,iCul),iCul,iFrame) = 0;
    end
end

save('simsig','video','-append');
% 
% writerObj = VideoWriter('generatePWV.avi','Uncompressed AVI');
% writerObj.FrameRate = fps;
% open(writerObj);
% 
% figure('units','pixels','position', [ 1 1 w h ], 'resize', 'off', 'toolbar', 'none' )
% 
% axis ij
% axis equal
% axis tight
% set(gca,'nextplot','replacechildren');
% set(gcf,'Renderer','zbuffer');
% for iFrame = 1 : nFrames
%     imagesc(video(:,:,iFrame)); colormap(gray)
%     axis equal
%     drawnow
%     frame = getframe;
%     writeVideo(writerObj,frame);
% end
% close(writerObj);