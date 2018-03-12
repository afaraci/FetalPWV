% Program to take audio data and convert it into a sonogram
%
%  Version 3, 23/10-96, by JAJ, after version by Henrik Klebaek

%  Find the number of samples

lines=max(size(audio));

%  Set parameters for calculating the sonogram

dN=floor(lines/1000+0.5);       %  Number of new samples for spectrum
Ndft=256;                       %  Samples in one spectrum
Nspec=floor((lines-Ndft)/dN);   %  Number of spectra
f0=2e6;                         %  Transducer center frequency
fprf=5e3;                       %  Pulse repetition frequency

%  Initialize data structure

power_dft=zeros(Ndft,Nspec);
power_min=[];
power_max=[];

disp('Making spectrum...')
for i=1:Nspec
  if (rem(i,50)==0)
    i
    end  
  sig=( audio((i-1)*dN+(1:Ndft)) )';
  sig=sig.*hanning(Ndft);
  dft=fftshift(abs(fft(sig)).^2);
  power_dft(:,i)=dft;
end

%  Scale power density

power_min=min(min(power_dft));
power_max=max(max(power_dft));
scale=128/(power_max-power_min);
power_dft=(power_dft-power_min)*scale;

power_dft=log(power_dft+0.005);
power_min=min(min(power_dft));
power_max=max(max(power_dft));
scale=128/(power_max-power_min);
power_dft=(power_dft-power_min)*scale;

%  Cut everything below a threshold

%power_dft = (power_dft>5) .* power_dft;

%  Make frequency axis

f=fprf/Ndft*((Ndft/2):-1:-(Ndft/2-1));
tidsakse=(1:Nspec)*dN/fprf;

%  Make the figure

figure(1);
image(tidsakse,f/1000,power_dft(Ndft:-1:1,:));		% Sonogram-plot
map=gray(128);
colormap(map(128:-1:1,:))
brighten(0.2)
%colormap(map)
%brighten(0.2)
title('Sonogram');					
xlabel('Time [s]');
ylabel('Frequency [kHz]');
axis([0 max(tidsakse) -2.5 2.5])

%print -deps sonogram.eps

