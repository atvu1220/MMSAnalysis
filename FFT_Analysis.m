%FFT Analysis of Bdata
%April 21, 2020
clear all
close all
data_type = 'brst';
date_start = strcat('2017-11-04 05:15:45.0');
date_end = strcat('2017-11-04 05:19:00.0');

formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(date_start,formatIn);
tend = datenum(date_end,formatIn);

[fgm_timedata, fgm_bdata] = load_fgm(date_start,date_end,1,data_type);
start_index = find(fgm_timedata > tstart, 1)-1;
end_index = find(fgm_timedata > tend, 1);

B = fgm_bdata(start_index:end_index,4);

if strcmp(data_type,'srvy')
    Fs = 16;
elseif strcmp(data_type,'brst')
    Fs = 128;
end

T = 1/Fs;
L = length(B);
t = (0:L-1)*T;

plot(t,B)

title('Magnetic Field Strength')
xlabel('t (s)')
ylabel('B(t)')

Y = fft(B);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
set(gcf,'color','w');
plot(f,P1,'linewidth',2); set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2) 
title('Fast Fourier Transform of B(t)','FontSize', 18, 'FontWeight', 'normal')
xlabel('f (Hz)','FontSize', 14)
ylabel('|P(f)|','FontSize', 14)
xlim([0 0.2])
peak_frequency = f(find(P1==max(P1)))