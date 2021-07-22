function [dynamic_pressure] = calculate_dynamic_pressure(fpi_ndata,fpi_vdata,specie)
    %inputs a single time step of data or multiple time steps of data and
    %calculates dynamic pressure
    
    if strcmp(specie,'i')
        
        mass_ion = 1.673*10^-27;%kg
        dynamic_pressure = 1/2*mass_ion.*fpi_ndata.*vecnorm(fpi_vdata,2,2).^2;
        dynamic_pressure = dynamic_pressure * 10^12*10^9; %to kg/ms^2,pascals, and then to nPA
    else
        
        mass_electron = 9.109*10^-31;%kg
        dynamic_pressure = 1/2*mass_electron.*fpi_ndata.*vecnorm(fpi_vdata,2,2).^2;
        dynamic_pressure = dynamic_pressure * 10^12*10^9; %to kg/ms^2,pascals, and then to nPA
        
    end
end

