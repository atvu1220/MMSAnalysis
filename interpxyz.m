function [subordinate_time_interp,subordinate_dataxyz_interp] = interpxyz(subordinate_time,subordinate_dataxyz,main_time)
    %interpolates all components, xyz
    %input r_time(:), rdata(:,xyz), maintime(:)
    columns = size(subordinate_dataxyz,2);
    
    subordinate_time_interp = zeros(length(main_time),1);
    subordinate_dataxyz_interp = zeros(length(main_time),columns);
    
    parfor i=1:columns
        %[x,index] = unique(fgm_timedata);
        %fgm_bdata_interp(:,i) = interp1(x,fgm_bdata(index,i),fpi_timedata);
        [x,index]= unique(subordinate_time);
        subordinate_dataxyz_interp(:,i) = interp1(x,subordinate_dataxyz(index,i),main_time);
        %subordinate_dataxyz_interp(:,i) = interp1(subordinate_time,subordinate_dataxyz(:,i),main_time);
    end
    
end

