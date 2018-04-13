function mat2avi( filename )

video = [];
cmppxx = [];
fps = [];
load( filename, 'video', 'cmppxx', 'fps' )

nFrames = size( video, 3 );
for i = 0 : 2
    videotmp( :, :, nFrames*i+1 : nFrames*(i+1) ) = video;
end
video = videotmp;
clear videotmp

[ h, w, f ] = size( video );

cmppxy = cmppxx;

cmx = cmppxx*( 1 : w );
cmy = cmppxy*( 1 : h );

writerObj = VideoWriter( [ filename '.avi' ], 'Uncompressed AVI' );
writerObj.FrameRate = fps;

open( writerObj );

close all

figure(  'units',    'pixels', ... 
      'position', [ 1 1 w h ], ...
        'resize',       'off', ...
        'toolbar',     'none' );

axis ij
axis equal
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
for i = 1 : f
    imagesc( cmx, cmy, video(:,:,i) ), colormap(gray)
    frame = getframe;
    for j = 1 : 3
        new_frame.cdata(:,:,j) = imresize( frame.cdata(:,:,j), [512 1024] );
    end
    frame.cdata = new_frame.cdata;
    writeVideo(writerObj,frame);
end
close( writerObj );
