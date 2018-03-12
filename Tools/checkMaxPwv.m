function maxPwv = checkMaxPwv

mainDir = 'C:\Users\Alessandro\OneDrive - King''s College London\Matlab\FetalPWV';
dataDir = fullfile( mainDir, 'Clinical data', 'ToUse' );
d = dir( dataDir );
listD = 1:length(d);
listD = listD([d.isdir]);
listD(1:2) = [];
i = 1;
for iDir = listD
    d2 = dir(fullfile( d(iDir).folder, d(iDir).name ));
%     listD2 = 1:length(d2);
%     listD2 = listD2(~[d2.isdir]);
%     listD2(1:2) = [];
    for iDir2 = 3 : length( d2 )
        fileName = fullfile( d2( iDir2 ).folder, d2( iDir2 ).name );
        if strcmpi( d2( iDir2 ).name(end-3:end),'.dcm' )
            info = dicominfo( fileName );
            if isfield( info, 'SequenceOfUltrasoundRegions')
                cmppxx = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;
                fps = info.CineRate;
                maxL = norm(double( [ info.Height info.Width ] ));
                maxPwv(i).name = fileName;
                maxPwv(i).di = maxL*cmppxx*10*fps/1000;
                maxPwv(i).wi = info.Width*cmppxx*10*fps/1000;
                i = i+1;
            end
        end
    end
end