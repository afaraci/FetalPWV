%  Load the calculated data and plot it


%  Physical data

f0=5e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
fprf=5e3;                %  Pulse emissions frequency  [Hz]
D=20;                    %  Sampling frequency decimation rate
no_lines=5000;           %  Number of lines for one direction

%  Load the data for one image line and a number of
%  pulse emissions

data=0;
data=zeros(250,no_lines);
for i=1:no_lines
  if (rem(i,20)==0)
    i
    end
  start_sample=6000;

  cmd=['load /data3/users_nbu/jaj/sim_flow_data/rf_data/rf_ln',num2str(i),'.mat'];
  eval(cmd);
  
    %  Decimate the data and store it in data
  
  rf_sig = rf_data (start_sample-tstart*fs:length(rf_data));
     
  rf_sig=hilbert(rf_sig(1:D:max(size(rf_sig))));
  data(1:250,i)=rf_sig(1:250);
  end

plot(real(data(125,:)))
drawnow
audio=data(125,:);
save audio.mat audio
sonogram