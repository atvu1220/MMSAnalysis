function [n] = tdnormal(date_start,date_end,fgm_timedata,fgm_bdata,flag)
    %Calculates the tangential discontinuity normal by the cross product of
    %the upstream and downstream magnetic fields divided by their
    %magnitude.
    
    if strcmp(flag,'all')
        %this is for the entire burst mode data file, 128 is per second of brst
        %data; 5 seconds average at beginning and end.
        n_cs = cross(mean(fgm_bdata(1:128*1,1:3)),mean(fgm_bdata(end-128*1:end,1:3)));
        n = n_cs/norm(n_cs);
        %n(1) = abs(n(1))*(-1);
        n(1) = abs(n(1));
    elseif strcmp(flag,'event')
        formatIn='yyyy-mm-dd HH:MM:SS.FFF';
        tstart = datenum(date_start,formatIn);
        tend = datenum(date_end,formatIn);
        
        %Find the start and end limits of the event in the data
        start_index = find(fgm_timedata > tstart, 1);
        end_index = find(fgm_timedata > tend, 1);
        
        %crop data
        fgm_bdata = fgm_bdata(start_index:end_index,:);
        %n_cs = bd cross bu / magnitude (bd, bu)
        n_cs = cross(mean(fgm_bdata(end-128*5:end,1:3)),mean(fgm_bdata(1:128*5,1:3))); %16 data points is srvy, 128 is brst
       
        n = n_cs/norm(n_cs); 
        
        %-1 for X component corresponds to pointing away from the sun as the
        %discontinuity convects across the bow shock in the antisunward direction, passing the
        %spacecraft. 
        %we force this in accordance to thomsen et al. 1993 as it is ambiguous for
        %which direction the normal should point, thus -x points towards
        %Earth.
        %n(1) = abs(n(1))*(-1);
        n(1) = abs(n(1)); %schwartz et al 2018 +x gse
    end
    
end

