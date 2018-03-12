function filtSig = filterWave( sig, StBfF, ord )

if nargin < 3
    ord = 0;
end

%%
aveSig = mean( sig );
sig = sig - aveSig;

Fs = 1000;              % Sampling frequency
L = length( sig );       % Length of signal

Y = fft( sig );
P2 = abs(Y/L);
P1 = P2( 1 : round(L/2)+1 );
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
[ ~, I ] = max(P1(2:end)); % Main frequency
I = I + 1;

CoF = f(I);
% dFilt = designfilt('lowpassfir', 'FilterOrder',ord, 'CutoffFrequency', CoF, 'SampleRate', 1000);
dFilt = designfilt('lowpassfir', 'FilterOrder',ord, 'PassbandFrequency', CoF, 'StopbandFrequency',CoF*StBfF, 'DesignMethod','ls','SampleRate', 1000);

% dFilt = designfilt('bandpassfir', 'FilterOrder', ord, ....
%     'CutoffFrequency1', 0.000000000000000001, 'CutoffFrequency2', CoF, 'SampleRate', Fs);

filtSig = filtfilt( dFilt, sig );
% miF = filter( dFilt, mi );
%  gd = grpdelay( dFilt );
% miF = miF( gd(1)+1 : end );
filtSig = filtSig + aveSig;