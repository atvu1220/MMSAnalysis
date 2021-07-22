function [timeStart,timeEnd] = getTimeRange(timeCenter,range,timeType)
    %Inputs center time and time range (in s or ms) and finds the start and end times, in timeType format
    format='yyyy-mm-dd HH:MM:SS.FFF';
    
    %Check what is the time range input
    if range > 30
        %Probably input is milliseconds
    else
        range = range * 1000; %Input is seconds, convert to milliseconds
    end
    
    
    
    if ~isa(timeCenter,'double') %check to see if it is datenum, if not, convert it.
        timeCenter = datenum(timeCenter,format);
    end
    
    timeStart = timeCenter - milliseconds(range/2);
    timeEnd = timeCenter + milliseconds(range/2);
    
    
    if ~strcmp(timeType,'datenum') %Convert to datatype if not datenum.
        timeStart = datestr(timeStart,format);
        timeEnd = datestr(timeEnd,format);
    end
    
end


