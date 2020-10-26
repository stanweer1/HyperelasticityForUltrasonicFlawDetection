% 2.82*pi*3076.40792391/(2*0.1*2*pi) = 21688.6758636
% 1.988*pi*3076.40792391/(2*0.1*2*pi) = 15289.7473818
% 2.63*pi*3076.40792391/(2*0.1*2*pi) = 20227.3820997
% sqrt(26.5*10^9/2800) = 3076.40792391

clc;
clear;

dt = 1e-7;
f = 20227;
T = 5/f;
t = 0:dt:T;
w = 0.5*(1 - cos(2*pi*f.*(t-1)/5.6)).*sin(2*pi*f.*(t-1));
plot(t,w)

% save('./5kHz.txt','t','w','-ASCII')
% title('5cycle Hanning Windowed Toneburst - 5kHz')
% for i = 1:1:round((T/dt))+1
%     wave = [t(1,i) w(1,i)];
%     save('./30kHz.txt','wave','-ASCII','-append')
% end

file = ['Hanning' num2str(f) 'Hz.txt'];
fid = fopen(file, 'w'); % In order to create a text file of the simulated signal
for iLine = 1:4:round((T/dt))+1
    fprintf(fid, '%e,%e,%e,%e,%e,%e,%e,%e\r\n',t(iLine),w(iLine),t(iLine+1),w(iLine+1),t(iLine+2),w(iLine+2),t(iLine+3),w(iLine+3));
end
fclose all;