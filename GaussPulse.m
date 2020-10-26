%----------------------------------------------------
% Department of Engineering Science and Mechanics
% Penn State University
% by Yang Liu in Nov. 2011
%----------------------------------------------------
%From this program, we can get the excitation signal for FE model
%----------------------------------------------------
clear;
clc;
f=871700;        %frequency of Sin wave
freqRange=0.5;
cutoff=-80;
tc = gauspuls('cutoff',f,freqRange,[],cutoff); 
L=4092;
% dt=1e-8;
dt = 2*tc/(L-1);    %sampling rate
t = -tc : dt : tc; 

i=0:dt:2*tc;
[yi yq]= gauspuls(t,f,freqRange); 
plot(i,yi)


% [m n]=size(yi)
% freq = 1/(2*dt)*linspace(0,1,L/2);
% z=abs(fft(yi))/max(abs(fft(yi)));
% 
% 
% figure;
% plot(freq/1e3,z(1:L/2));

file = ['Gauss2_' num2str(f) 'Hz.txt']
fid = fopen(file, 'w'); % In order to create a text file of the simulated signal
for iLine = 1:4:length(yq)
    fprintf(fid, '%e,%e,%e,%e,%e,%e,%e,%e\r\n',i(iLine),yi(iLine),i(iLine+1),yi(iLine+1),i(iLine+2),yi(iLine+2),i(iLine+3),yi(iLine+3));
end
fclose(fid);
