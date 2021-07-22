%Calculates the Merka Bow Shock model from Merka et al. 2005 using Mach Numbers
% and finds the itnersection of B vector on BS to get the shock angle
function [shockNormal,shockAngle] = calculate_BBSShockAngle(...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms1_fgm_timedata_raw,...
        outsideTime,B,...
        n,v,M_A)
    %% Positions of Event
    
    %Interpolate the Data
    [~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    %Find center position of MMS within leading boundary
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    outsideTime_index = find(mms1_fgm_timedata_raw >= datenum(outsideTime,formatIn), 1);
    
    r0 = mms1_mec_rdata_interp(outsideTime_index,:)./6371.2;
    B = B./norm(B);
    
    r0 = rotate_GSE2GPE(v,r0);
    B_rotated = rotate_GSE2GPE(v,B);
    syms t
    Bvector_x = r0(1) + B_rotated(1)*t;
    Bvector_y = r0(2) + B_rotated(2)*t;
    Bvector_z = r0(3) + B_rotated(3)*t;
    
    %% Dynamic Pressure for Standoff Distance
    ScalingFactor = (n * norm(v)^2 / (7 * 457.5^2))^(1/6);
    
    %% Merka Bow Shock
    %Use fitted model from Merka for the Coefficients of B
    b1 = [0.0063, -0.0098];
    b3 = [0.8351,0.0102];
    b4 = [-0.0298,0.0040];
    b7 = [16.39, 0.2087,108.3];
    b8 = [-0.9241,0.0721];
    b10 = [-444,-2.935,-1930];
    
    
    a1  =  b1(1)  +  b1(2)*M_A;
    a3  =  b3(1)  +  b3(2)*M_A;
    a4  =  b4(1)  +  b4(2)*M_A;
    a7  =  b7(1)  +  b7(2)*M_A  + b7(3)/(M_A-1)^2;
    a8  =  b8(1)  +  b8(2)*M_A;
    a10 =  b10(1) +  b10(2)*M_A + b10(3)/(M_A-1)^2;
  
    syms x y z
    merka_bowshock = a1.*(ScalingFactor.*x).^2 + (ScalingFactor.*y).^2 + a3.*(ScalingFactor.*z).^2 + 2.*a4.*(ScalingFactor.*x).*(ScalingFactor.*y) + 2.*a7.*(ScalingFactor.*x) + 2.*a8.*(ScalingFactor.*y) + a10;

   
    
    %% Calculate Intersection of B and the Bow Shock
    F = subs(merka_bowshock,[x y z],[Bvector_x Bvector_y Bvector_z]);
    T = solve(F==0,t);
    
    BS_pos_x = double(subs(Bvector_x,t,T));
    BS_pos_y = double(subs(Bvector_y,t,T));
    BS_pos_z = double(subs(Bvector_z,t,T));
    BS_pos = [BS_pos_x BS_pos_y BS_pos_z];
    
    %Check the closer to HFA observation position.
    [~,I] = min(vecnorm((BS_pos - r0'),2,2));
    BS_pos = BS_pos(I,:);
     
    
    %% Calculate Bow Shock Normal at Connection Point
    
    % Calculation from Chu 2017 Thesis, p.62
    shockNormal = [2*a1*BS_pos(1) + 2*a4*BS_pos(2) + 2*a7,...
        2*BS_pos(2) + 2*a4*BS_pos(1) + 2*a8,...
        2*a3*BS_pos(3)];
    shockNormal = shockNormal ./ norm(shockNormal);
    
    shockNormal = rotate_GPE2GSE(v,shockNormal);
    BS_pos = rotate_GPE2GSE(v,BS_pos);
    
    BS_pos = BS_pos';
    shockNormal = shockNormal';
    shockAngle = acosd(dot(B,shockNormal)/norm(B));
end