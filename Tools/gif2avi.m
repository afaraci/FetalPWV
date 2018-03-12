function gif2avi(fileName)

close all

info = dicominfo(fileName);
w = double( info.Width  );
h = double( info.Height );
cmppxx = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;
cmppxy = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY;
cmx = cmppxx*[0:(w-1)];
cmy = cmppxy*[0:(h-1)];

X = dicomread(fileName);
X = squeeze(X);

writerObj = VideoWriter([ fileName '.avi'],'Uncompressed AVI');

open(writerObj);

figure('units','pixels','position', [ 1 1 w h ], 'resize', 'off', 'toolbar', 'none' )

axis ij
axis equal
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
for i = 1:size(X,3)
    imagesc( cmx, cmy, X(:,:,i) ), colormap(gray)
    frame = getframe;
    writeVideo(writerObj,frame);
end

close(writerObj);

writerObj.FrameRate = info.FrameTime;