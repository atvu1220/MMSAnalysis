function [magnetic_pressure] = calculate_magnetic_pressure(fgm_bdata)
    %inputs a single time step of data or multiple time steps of data and
    %calculates magnetic pressure
    magnetic_pressure = 10^-9*fgm_bdata(:,4).^2 / (2*1.25663706e-6); %in nPa
end

