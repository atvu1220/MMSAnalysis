

function [normal,r_sc] = calculate_bowshocknormal(date_start,mec_rdata,fpi_timedata,fpi_ndata,fpi_vdata)
    %Inputs the position of the spacecraft in GSE and solar wind speed and uses the slavin
    %holzer mean model to shift the bow shock to locally at the current
    %spacecraft observation position of the event. Then it calculates the
    %gradient for the bow shock normal. Outputs the bow shock normal.
    %refer to schwartz 1998 or Analysis Method for Multispacecraft data
    %10.4.6
    %uses the dynamic pressure for scaling the bow shock model to the
    %observation point
    
    %Perhaps will complete later if needed for accuracy
    %Using the Slavin Holzer Mean Method, established parameters from
    %stastical studies of bow shock
    % Alternatively scale these
    % parameters so that the model passes through a given position vector r crossing. This
    % process can be reduced to the substitutions L ! L, ro ! ro, and r ! r crossing
    % in equations 10.22 and 10.23 and solving the resulting quadratic equation for . In
    % the case of hyperbolic ( > 1) models, the larger of the two roots corresponds to the
    % correct branch of the hyperbola.
    
    
    
    
    %%%%%%%%%%%%%%%%%%%Dynamic Pressure before event Calculation%%%%%%%%%%%%%%%
    %Calculate the dynamic pressure to scale model
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
   
% %     %Find the start and end limits of the event in the data
% %     start_index = find(fpi_timedata > tstart, 1);
% %     
% %     %crop data to just the first second of data  and then calculate the
% %     %mean
% %     fpi_ndata = mean(fpi_ndata(start_index:start_index+128,:));
% %     fpi_vdata = mean(fpi_vdata(start_index:start_index+128,:));
    
    start_index=1; %use the first data point in data, brst
    fpi_ndata = fpi_ndata(start_index,:);
    fpi_vdata = fpi_vdata(start_index,:);
    
    %Calculate dynamic pressure rho*v^2 with mean n and v
    dynamic_pressure = calculate_dynamic_pressure(fpi_ndata,fpi_vdata,'i');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %Calculate needed parameters for bow shock model fit
    r_sc= mec_rdata(1,:)/6371.2; %spacecraft position at observation
    v_sw = norm(fpi_vdata); %solar wind speed at observation, mean
   
    
    %Begin bow shock model fit
    p_dam = 2.1;%nPA
    sigma = (dynamic_pressure/p_dam)^(-1/6); %Scaling Factor
    
    epsilon= 1.16;
    L=sigma*23.3; %Re
    x_0 = sigma*3.0; %Re
    y_0 = 0; %Re
    z_0 = 0; %Re
    
    alpha = atan(30/v_sw); % solar wind in km/s
    
    %Variable transformation from r_sc to r^{abd}
    %Transformation Matrix
    rotation_matrix = ...
        [cos(alpha) -sin(alpha)  0;
        sin(alpha)  cos(alpha)  0;
        0           0       1];
    
    %New position vector in the model coordinate system, xyz
    r_abd = rotation_matrix * [r_sc(1);r_sc(2);r_sc(3)] - [x_0; y_0; z_0];
    
    %Surface equation
    r_abd_mag = norm(r_abd);
    S = (r_abd_mag + epsilon*r_abd(1))^2 - L^2;
    
    %calculate gradient and rotated back into the unaberrated coordinate
    %frame
    grad_S = ( 2*L / r_abd_mag) * ...
        [(r_abd(1)*(1-epsilon^2) + epsilon*L)*cos(alpha) + r_abd(2)*sin(alpha);...
        -(r_abd(1)*(1-epsilon^2) + epsilon*L)*sin(alpha) + r_abd(2)*cos(alpha);...
        r_abd(3)
        ];
    
    %Magnitude of vector
    grad_S_mag = norm(grad_S);
    
    %Calculate the normal of the curved surface, S(r^{abd})
    normal = (grad_S/grad_S_mag)';
    
end

