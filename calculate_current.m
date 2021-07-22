function [current] = calculate_current(fpi_e_timedata,fpi_i_timedata,fpi_e_ndata,fpi_e_vdata,fpi_i_vdata)
    
    %Downsample e
    [~,fpi_e_ndata] = interpxyz(fpi_e_timedata,fpi_e_ndata,fpi_i_timedata);
    [~,fpi_e_vdata] = interpxyz(fpi_e_timedata,fpi_e_vdata,fpi_i_timedata);
    
    q = 1.60e-19;
    %Units of mA/km^2
    current = fpi_e_ndata.*q.*(fpi_i_vdata-fpi_e_vdata).*10^18;
end

