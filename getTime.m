function [time1,time2] = getTime(time1,time2)
    %Converts times: datestring to datenum and datenum to Datestring
    format='yyyy-mm-dd HH:MM:SS.FFF';
    
    if nargin == 1
        if isa(time1,'double')
            time1 = datestr(time1,format);
        else
            time1 = datenum(time1,format);
        end
    elseif nargin==2
        if isa(time1,'double')
            time1 = datestr(time1,format);
            time2 = datestr(time2,format);
        else
            time1 = datenum(time1,format);
            time2 = datestr(time2,format);
        end
    end
    
    
end


