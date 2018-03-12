function [ rf_lines, tstart_lines, image_width ] = ... 
    generatePWVphantomImage( pwv, fps, nLines )

if nargin < 1
    nLines = 30;
end

%field_init(-1);
image_width = 40; % [mm]
[ Pt, phantom_amplitudes ] = generatePWVphantom( pwv, fps );
nFrames = size( Pt, 3 );

%%
    rf_lines =  cell( nLines, nFrames );
tstart_lines = zeros( nLines, nFrames );
parfor iFrame = 1 : nFrames
    %field_init(-1);
    disp(iFrame)
    phantom_positions = Pt(:,:,iFrame);
    [ rf_lines(:,iFrame), tstart_lines(:,iFrame) ] = ...
        sim_img( phantom_positions, phantom_amplitudes, ...
        image_width, nLines );
end