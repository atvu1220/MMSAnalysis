function [duration] = findDuration(event_start,event_end)
    %Takes in two datetime strings and outputs the number of seconds as 
    %double type between the two times.
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    tstart = datenum(event_start,formatIn);
    tend = datenum(event_end,formatIn);
    
    duration = (tend-tstart)*86400;
end

