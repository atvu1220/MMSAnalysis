function [thermal_pressure] = calculate_thermal_pressure(fpi_pressure)
    %inputs a single time step of data or multiple time steps of data and
    %calculates thermal pressure
    
    
    
    %Diagonalize thermal pressure tensor
    fpi_diagpressure = zeros(length(fpi_pressure),3);
    for i=1:length(fpi_pressure)
        fpi_diagpressure(i,:) = (eig(fpi_pressure(:,:,i)))';
    end
    
    %Trace of square matrix or sum of eigenvalues times 1/3
    thermal_pressure = (1/3).*sum(fpi_diagpressure,2);
    
end

