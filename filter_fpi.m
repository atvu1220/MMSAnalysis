function [signalFiltered] =  filter_fpi(signal,specie,hampelWindow,hampelSTD)
     if strcmp(specie,'e')
         frequency = 1/0.03; %hz
         LowPassNoiseCutoff = 0.15; %hz
     elseif strcmp(specie,'i')
         LowPassNoiseCutoff = 1; %hz
         frequency = 1/0.15; %hz
     end
     y = fft(signal);
     ymag = (y.*conj(y)).^(1/2);
%      ymag = real(y);
     yphase = atan2(imag(y),real(y));
     [yymag,j,xmedian,xsigma] = hampel(ymag,hampelWindow,hampelSTD);
     yyymag = lowpass(yymag,LowPassNoiseCutoff,frequency);
%      yyymag = yyymag.^(1/2);
     [x,y] = pol2cart(yphase,yyymag);
     
     yFiltered = complex(x,y);
     signalFiltered = ifft(yFiltered,'symmetric');
     
%      [xx,yy] = pol2cart(yphase,ymag)
%      ystar = complex(xx,yy)
%      close all
%      figure
%      plot(signal)
%      hold on
%      plot(signalFiltered)
 end
 