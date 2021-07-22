%Calculates the Merka Bow Shock model from Merka et al. 2005 using Mach Numbers
function [leading_DistancetoConnectionPoint,leading_shock_normal,...
        trailing_DistancetoConnectionPoint,trailing_shock_normal,...
        leading_DistanceFromIntersectiontoBS,trailing_DistanceFromIntersectiontoBS,...
        event_width,event_height] = calculate_MerkaBS_EventEdges(...
        mms1_fgm_timedata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
        event_start,event_end,...
        left_InnerEdge,right_InnerEdge,...
        leading_start,leading_end,...
        trailing_start,trailing_end,...
        timing_n1,timing_n2,...
        n,v,M_A,timing_v1)
    %% Positions of Event
    
    %Interpolate the Data
    [~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw);
    [~,mms2_mec_rdata_interp] = interpxyz(mms2_mec_timedata_raw,mms2_mec_rdata_raw,mms1_fgm_timedata_raw);
    [~,mms3_mec_rdata_interp] = interpxyz(mms3_mec_timedata_raw,mms3_mec_rdata_raw,mms1_fgm_timedata_raw);
    [~,mms4_mec_rdata_interp] = interpxyz(mms4_mec_timedata_raw,mms4_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    %Find center position of MMS within leading boundary
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    leading_start_index = find(mms1_fgm_timedata_raw >= datenum(leading_start,formatIn), 1);
    leading_end_index = find(mms1_fgm_timedata_raw >= datenum(leading_end,formatIn), 1);
    
    CenterofLeading_TetrahedronPosition = mean([...
        mms1_mec_rdata_interp(leading_start_index:leading_end_index,:);...
        mms2_mec_rdata_interp(leading_start_index:leading_end_index,:);...
        mms3_mec_rdata_interp(leading_start_index:leading_end_index,:);...
        mms4_mec_rdata_interp(leading_start_index:leading_end_index,:)...
        ]);
    
    %Find center position of MMS within core of event
    left_InnerEdge_index = find(mms1_fgm_timedata_raw >= datenum(left_InnerEdge,formatIn), 1);
    right_InnerEdge_index = find(mms1_fgm_timedata_raw >= datenum(right_InnerEdge,formatIn), 1);
    
    CenterofCore_TetrahedronPosition = mean([...
        mms1_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms2_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms3_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms4_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:)...
        ]);
    
    %Find center position of MMS within leading boundary
    trailing_start_index = find(mms1_fgm_timedata_raw >= datenum(trailing_start,formatIn), 1);
    trailing_end_index = find(mms1_fgm_timedata_raw >= datenum(trailing_end,formatIn), 1);
    
    CenterofTrailing_TetrahedronPosition = mean([...
        mms1_mec_rdata_interp(trailing_start_index:trailing_end_index,:);...
        mms2_mec_rdata_interp(trailing_start_index:trailing_end_index,:);...
        mms3_mec_rdata_interp(trailing_start_index:trailing_end_index,:);...
        mms4_mec_rdata_interp(trailing_start_index:trailing_end_index,:)...
        ]);
    
    
    
    
    trailing_position=CenterofTrailing_TetrahedronPosition./6371.2;
    leading_position=CenterofLeading_TetrahedronPosition./6371.2;
%     core_position=CenterofCore_TetrahedronPosition./6371.2;
    
    
    %Include new position of leading edge at time of trailing edge position observation
    coreDuration = findDuration(left_InnerEdge,right_InnerEdge);
    leading_velocity = timing_n1.*timing_v1./6371.2;
%     leading_position = leading_position + coreDuration.*leading_velocity;
    core_position = mean([leading_position; trailing_position]);
    
    %% Dynamic Pressure for Standoff Distance
    ScalingFactor = (n * norm(v)^2 / (7 * 457.5^2))^(1/6);
    
    %% Merka Bow Shock
    %Use fitted model from Merka for the Coefficients of B
    % %     b1 = [0.0063, -0.0098];
    % %     b3 = [0.8351,0.0102];
    % %     b4 = [-0.0298,0.0040];
    % %     b7 = [16.39, 0.2087,108.3];
% %     b8 = [-0.9241,0.0721];
% %     b10 = [-444,-2.935,-1930];
% %     
% %     
% %     a1  =  b1(1)  +  b1(2)*M_A;
% %     a3  =  b3(1)  +  b3(2)*M_A;
% %     a4  =  b4(1)  +  b4(2)*M_A;
% %     a7  =  b7(1)  +  b7(2)*M_A  + b7(3)/(M_A-1)^2;
% %     a8  =  b8(1)  +  b8(2)*M_A;
% %     a10 =  b10(1) +  b10(2)*M_A + b10(3)/(M_A-1)^2;
    
    %Manually a-coefficients to those found in Merka
    a1 = [-0.12292,-0.02996,-0.21644,0.04615];
    a3 = [0.81092,0.93743,0.91729,1.00468];
    a4 = [-0.20902,0.00227,0.01680,0.01286];
    a7 = [28.27,21.196,19.570,20.280];
    a8 = [-0.75985,-0.57751,0.02222,-0.03024];
    a10 = [-649.29,-536.81,-471.21,-512.18];
    
    if M_A <= 5
        a_index = 1;
    elseif M_A <= 9
        a_index = 2;
    elseif M_A <= 14
        a_index = 3;
    else
        a_index = 4;
    end
    a1 = a1(a_index);
    a3 = a3(a_index);
    a4 = a4(a_index);
    a7 = a7(a_index);
    a8 = a8(a_index);
    a10 = a10(a_index);
    
    syms x y z
    %merka_bowshock = a1.*x.^2 + y.^2 + a3.*z.^2 + 2.*a4.*x.*y + 2.*a7.*x + 2.*a8.*y + a10;
    
    merka_bowshock = a1.*(ScalingFactor.*x).^2 + (ScalingFactor.*y).^2 + a3.*(ScalingFactor.*z).^2 + 2.*a4.*(ScalingFactor.*x).*(ScalingFactor.*y) + 2.*a7.*(ScalingFactor.*x) + 2.*a8.*(ScalingFactor.*y) + a10;

    %% Spacecraft Trajectory Plane
    %Rotate Vectors from GSE to GPE
    trailing_position=rotate_GSE2GPE(v,trailing_position);
    core_position = rotate_GSE2GPE(v,core_position);
    leading_position = rotate_GSE2GPE(v,leading_position);
    
    trailing_core_vector = trailing_position - core_position;
    leading_core_vector = leading_position - core_position;
    
    cross_trailing_leading = cross(trailing_core_vector,leading_core_vector);
    
    spacecraft_trajectory_plane = ...
        cross_trailing_leading(1).*(x-core_position(1))...
        + cross_trailing_leading(2).*(y-core_position(2))...
        + cross_trailing_leading(3).*(z-core_position(3));
    
    %% Event Edges Plane
    timing_n1 = rotate_GSE2GPE(v,timing_n1);
    timing_n2 = rotate_GSE2GPE(v,timing_n2);
    
    leading_plane = ...
        timing_n1(1).*(x-leading_position(1))...
        + timing_n1(2).*(y-leading_position(2))...
        + timing_n1(3).*(z-leading_position(3));
    
    trailing_plane = ...
        timing_n2(1).*(x-trailing_position(1))...
        + timing_n2(2).*(y-trailing_position(2))...
        + timing_n2(3).*(z-trailing_position(3));
    
    %% Calculate connection point on BS
    
%     [leading_DistancetoConnectionPoint,leading_bs_pos] = calculate_ConnectionPoint (merka_bowshock,spacecraft_trajectory_plane,leading_plane,leading_position,a1,a3,a4,a7,a8);
%     [trailing_DistancetoConnectionPoint,trailing_bs_pos] = calculate_ConnectionPoint (merka_bowshock,spacecraft_trajectory_plane,trailing_plane,trailing_position,a1,a3,a4,a7,a8);
   
    [leading_DistancetoConnectionPoint,leading_bs_pos] = calculate_ConnectionPoint (merka_bowshock,leading_plane,timing_n1,leading_position,a1,a3,a4,a7,a8,ScalingFactor);
    [trailing_DistancetoConnectionPoint,trailing_bs_pos] = calculate_ConnectionPoint (merka_bowshock,trailing_plane,timing_n2,trailing_position,a1,a3,a4,a7,a8,ScalingFactor);
    
    %% Calculate Bow Shock Normal at Connection Point
    [leading_shock_normal] = calculate_shockNormal(a1,a3,a4,a7,a8,leading_bs_pos,ScalingFactor);
    [trailing_shock_normal] = calculate_shockNormal(a1,a3,a4,a7,a8,trailing_bs_pos,ScalingFactor);
    
    %% Calculate plane defined by the two (closest) connection points
    leading_bs_core_vector = leading_bs_pos' - core_position;
    trailing_bs_core_vector = trailing_bs_pos' - core_position;
    
    
    cross_trailing_leading_bs = cross(leading_bs_core_vector,trailing_bs_core_vector);
    
    bs_connections_plane = ...
        cross_trailing_leading_bs(1).*(x-core_position(1))...
        + cross_trailing_leading_bs(2).*(y-core_position(2))...
        + cross_trailing_leading_bs(3).*(z-core_position(3));
    
    %% Calculate Intersection Point between Event Edges Plane
    %      [intersection_point,leading_DistanceFromIntersectiontoBS,trailing_DistanceFromIntersectiontoBS] = ...
    %          calculate_EdgesPlaneIntersection (spacecraft_trajectory_plane,leading_plane,leading_bs_pos,trailing_plane,trailing_bs_pos);
    
    [intersection_point,leading_DistanceFromIntersectiontoBS,trailing_DistanceFromIntersectiontoBS] = ...
        calculate_EdgesPlaneIntersection (bs_connections_plane,leading_plane,leading_bs_pos,trailing_plane,trailing_bs_pos);
    
    %% Calculate Height of Event from BS
    midpoint_bs_pos = (leading_bs_pos + trailing_bs_pos)./2;
    event_width = norm(leading_bs_pos - trailing_bs_pos);
    %event_height = norm(intersection_point-midpoint_bs_pos);
    
    % Substract portion of height that is within the bow shock
    if subs(merka_bowshock,[x,y,z],midpoint_bs_pos) < 0 %If midpoint is inside, proceed
        
        
        
        if subs(merka_bowshock,[x,y,z],intersection_point) > 0 %If intersection is outside, proceed %norm(intersection_point) > norm(midpoint_bs_pos)
            dt = 0.001;
            midpointToIntersect = (intersection_point-midpoint_bs_pos)./norm(intersection_point-midpoint_bs_pos);
            currentPosition = midpoint_bs_pos;
            
            while subs(merka_bowshock,[x,y,z],currentPosition) < 0
                currentPosition = currentPosition + midpointToIntersect*dt;
            end
            
            event_height = norm(intersection_point-currentPosition); %Intersection outside,midpoint inside, take distance from BS to Intersection point
            
        else
            event_height = norm(intersection_point-midpoint_bs_pos); %Intersection point is inside
        end
        
        
        
    else
        event_height = norm(intersection_point-midpoint_bs_pos); %Midpoint is outside
    end
    
    
    
%     %% Rotate Vectors back to GSE Coordinates
%     leading_shock_normal = rotate_GPE2GSE(v,leading_shock_normal);
%     leading_bs_pos = rotate_GPE2GSE(v,leading_bs_pos);
%     leading_bs_pos = leading_bs_pos';
%     leading_shock_normal = leading_shock_normal';
%     
%     trailing_shock_normal = rotate_GPE2GSE(v,trailing_shock_normal);
%     trailing_bs_pos = rotate_GPE2GSE(v,trailing_bs_pos); 
%     trailing_bs_pos = trailing_bs_pos';
%     trailing_shock_normal = trailing_shock_normal';
%     
%     core_position = rotate_GPE2GSE(v,core_position);
%     intersection_point = rotate_GPE2GSE(v,intersection_point);
    %% Plot Configuration
    figure
    Zmerka_bowshock = solve(merka_bowshock==0,z);
    fsurf(matlabFunction(Zmerka_bowshock(1)))
    hold on
    fsurf(matlabFunction(Zmerka_bowshock(2)))
    
        zbs_connections_plane = solve(bs_connections_plane==0,z);
        fsurf(matlabFunction(zbs_connections_plane))
    
    
    plot3(core_position(1),core_position(2),core_position(3),'.r','markersize',100)
    plot3(leading_bs_pos(1),leading_bs_pos(2),leading_bs_pos(3),'.m','markersize',150)
    
    zleading_plane = solve(leading_plane==0,z);
    fsurf(matlabFunction(zleading_plane))
    plot3(leading_bs_pos(1),leading_bs_pos(2),leading_bs_pos(3),'.m','markersize',150)
    
    ztrailing_plane = solve(trailing_plane==0,z);
    fsurf(matlabFunction(ztrailing_plane))
    plot3(trailing_bs_pos(1),trailing_bs_pos(2),trailing_bs_pos(3),'.c','markersize',150)
    
    plot3(intersection_point(1),intersection_point(2),intersection_point(3),'.y','markersize',150)
    
    hold off
    axis([-10 50 -50 50 -50 50])
    view([150 47])
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Event Edges - BS Configuration')
    plot_name =  strcat('EventEdges_Configuration');
    %print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
    print(gcf,'-dpng','-r300', '-loose', plot_name);
    savefig(plot_name)
    
    
    %% Rotate Vectors back to GSE Coordinates
    leading_shock_normal = rotate_GPE2GSE(v,leading_shock_normal);
    leading_bs_pos = rotate_GPE2GSE(v,leading_bs_pos);
    leading_bs_pos = leading_bs_pos';
    leading_shock_normal = leading_shock_normal';
    
    trailing_shock_normal = rotate_GPE2GSE(v,trailing_shock_normal);
    trailing_bs_pos = rotate_GPE2GSE(v,trailing_bs_pos);
    trailing_bs_pos = trailing_bs_pos';
    trailing_shock_normal = trailing_shock_normal';
    
    core_position = rotate_GPE2GSE(v,core_position);
    intersection_point = rotate_GPE2GSE(v,intersection_point);
end
%% Calculate Event Edges Intersection
function [intersection_point,Leading_DistanceFromIntersectiontoBS,Trailing_DistanceFromIntersectiontoBS] = calculate_EdgesPlaneIntersection (bs_connections_plane,leading_plane,leading_bs_pos,trailing_plane,trailing_bs_pos)
    syms x y z 
    X = solve(leading_plane==trailing_plane,x);
    A = subs(leading_plane,x,X);
    Y = solve(A,y);
    
    B = subs(bs_connections_plane,x,X);
    C = subs(B,y,Y);
    ZZ = solve(C,z);
    
    YY = subs(Y,z,ZZ);
    XX = subs(X,z,ZZ);
    
    XX = subs(XX,y,YY);
    
    x_intersect = double(XX);
    y_intersect = double(YY);
    z_intersect = double(ZZ);
    
    intersection_point = [x_intersect,y_intersect,z_intersect];
    
    %Check the closer to HFA observation position.
    Leading_DistanceFromIntersectiontoBS = norm(intersection_point - leading_bs_pos);
    Trailing_DistanceFromIntersectiontoBS = norm(intersection_point - trailing_bs_pos);
end
%% Calculate Event Edge Connection with BS
% function [DistancetoConnectionPoint,bs_pos] = calculate_ConnectionPoint (merka_bowshock,spacecraft_trajectory_plane,event_edge_plane,edge_position,a1,a3,a4,a7,a8)
function [DistancetoConnectionPoint,bs_pos] = calculate_ConnectionPoint (merka_bowshock, event_edge_plane, edge_normal,core_position,a1,a3,a4,a7,a8,ScalingFactor)

        
    syms x y z lambda
    dLdx = 2*x - 2*core_position(1) + 2*lambda*a1*(ScalingFactor*x) + 2*lambda*a4*(ScalingFactor*y) + 2*lambda*a7 - edge_normal(1);
    dLdy = 2*y - 2*core_position(2) + 2*lambda*(ScalingFactor*y) + 2*lambda*a4*(ScalingFactor*x) + 2*lambda*a8 - edge_normal(2);
    dLdz = 2*z - 2*core_position(3) + 2*lambda*a3*(ScalingFactor*z) - edge_normal(3);
    
    
    
    Z = solve(dLdz == 0,z);
    Y = solve(dLdy == 0,y);
    X = subs(dLdx,y,Y);
    X = solve(X==0,x);
    
    A = subs(merka_bowshock-event_edge_plane,z,Z);
    A = subs(A,y,Y);
    A = subs(A,x,X);
    Lambda = solve(A==0,lambda);
    Y = subs(Y,x,X);
    
    ZZ = subs(Z,lambda,Lambda);
    YY = subs(Y,lambda,Lambda);
    XX = subs(X,lambda,Lambda);
    
    %Get rid of imaginary values
    bs_pos = [double(XX),double(YY),double(ZZ)];
    
    
    bs_pos2 = [];
    for i=1:size(bs_pos,1)
        if imag(bs_pos(i,:)) == 0
            bs_pos2 = [bs_pos2; bs_pos(i,:)];
        end
    end
    
    bs_pos = bs_pos2;
%     bs_pos = bs_pos( imag(bs_pos)==0);
    
    
%     bs_pos = reshape(bs_pos,[],3);
    
    %Get the closest point on the bow shock
    [DistancetoConnectionPoint,I] = min(vecnorm((bs_pos - core_position'),2,2));
    bs_pos = bs_pos(I,:);
    
    
    
% % % % %     %%%
% % % % %     syms x y z
% % % % %     X = solve(event_edge_plane==spacecraft_trajectory_plane,x);
% % % % %     A = subs(event_edge_plane,x,X);
% % % % %     Y = solve(A,y);
% % % % %     
% % % % %     B = subs(merka_bowshock,x,X);
% % % % %     C = subs(B,y,Y);
% % % % %     ZZ = solve(C,z);
% % % % %     
% % % % %     YY = subs(Y,z,ZZ);
% % % % %     XX = subs(X,z,ZZ);
% % % % %     
% % % % %     XX(1) = subs(XX(1),y,YY(1));
% % % % %     XX(2) = subs(XX(2),y,YY(2));
% % % % %     
% % % % %     x_bs = double(XX);
% % % % %     y_bs = double(YY);
% % % % %     z_bs = double(ZZ);
% % % % %     
% % % % %     bs_pos = [x_bs,y_bs,z_bs];
% % % % %     
% % % % %     %Check the closer to HFA observation position.
% % % % %     [DistancetoConnectionPoint,I] = min(vecnorm((bs_pos - edge_position'),2,2));
% % % % %     bs_pos = bs_pos(I,:);
% % % % %     
% % % % %     
% % % % %     if imag(bs_pos) ~= 0
% % % % %         
% % % % %         ZZ = solve(diff(Y,z)==diff(C,z));
% % % % %         YY = subs(Y,z,ZZ);
% % % % %         XX = subs(X,[y z],[YY ZZ]);
% % % % %         
% % % % %         edge_closest_position = double([XX;YY;ZZ]);
% % % % %         
% % % % %         %Solve for the minimum distance from the SHFA to the BS with lagrange multipliers
% % % % %         syms lambda
% % % % %         dLdx = 2*x - 2*edge_closest_position(1) + 2*lambda*a1*x + 2*lambda*a4*y + 2*lambda*a7;
% % % % %         dLdy = 2*y - 2*edge_closest_position(2) + 2*lambda*y + 2*lambda*a4*x + 2*lambda*a8;
% % % % %         dLdz = 2*z - 2*edge_closest_position(3) + 2*lambda*a3*z;
% % % % %         
% % % % %         Z = solve(dLdz == 0,z);
% % % % %         Y = solve(dLdy == 0,y);
% % % % %         X = subs(dLdx,y,Y);
% % % % %         X = solve(X==0,x);
% % % % %         
% % % % %         
% % % % %         A = subs(merka_bowshock,z,Z);
% % % % %         A = subs(A,y,Y);
% % % % %         A = subs(A,x,X);
% % % % %         Lambda = solve(A==0,lambda);
% % % % %         Y = subs(Y,x,X);
% % % % %         
% % % % %         ZZ = subs(Z,lambda,Lambda);
% % % % %         YY = subs(Y,lambda,Lambda);
% % % % %         XX = subs(X,lambda,Lambda);
% % % % %         
% % % % %         %Get rid of imaginary values
% % % % %         bs_pos = [double(XX),double(YY),double(ZZ)];
% % % % %         bs_pos = bs_pos( imag(bs_pos)==0);
% % % % %         
% % % % %         
% % % % %         bs_pos = reshape(bs_pos,[],3);
% % % % %         
% % % % %         %Get the closest point on the bow shock
% % % % %         [DistancetoConnectionPoint,I] = min(vecnorm((bs_pos - edge_position'),2,2));
% % % % %         bs_pos = bs_pos(I,:);
% % % % %     end
    
end
%% Calculate Bow Shock Normal at Connection Point
function [shock_normal] = calculate_shockNormal(a1,a3,a4,a7,a8,bs_pos,ScalingFactor)
    % Calculation from Chu 2017 Thesis, p.62
    shock_normal = [2*a1*(ScalingFactor*bs_pos(1)) + 2*a4*(ScalingFactor*bs_pos(2)) + 2*(ScalingFactor*a7),...
        2*bs_pos(2) + 2*a4*(ScalingFactor*bs_pos(1)) + 2*(ScalingFactor*a8),...
        2*a3*(ScalingFactor*bs_pos(3))];
    shock_normal = shock_normal ./ norm(shock_normal);
end
%% Calculate Closest Distance to BS
function [ClosestDistance,closest_point] = calculate_closestDistance(merka_bowshock, core_position,a1,a3,a4,a7,a8,ScalingFactor)

    %Solve for the minimum distance from the HFA to the BS with lagrange multipliers
    syms x y z lambda
    dLdx = 2*x - 2*core_position(1) + 2*lambda*a1*(ScalingFactor*x) + 2*lambda*a4*(ScalingFactor*y) + 2*lambda*a7;
    dLdy = 2*y - 2*core_position(2) + 2*lambda*(ScalingFactor*y) + 2*lambda*a4*(ScalingFactor*x) + 2*lambda*a8;
    dLdz = 2*z - 2*core_position(3) + 2*lambda*a3*(ScalingFactor*z);
    
    Z = solve(dLdz == 0,z);
    Y = solve(dLdy == 0,y);
    X = subs(dLdx,y,Y);
    X = solve(X==0,x);
    
    
    A = subs(merka_bowshock,z,Z);
    A = subs(A,y,Y);
    A = subs(A,x,X);
    Lambda = solve(A==0,lambda);
    Y = subs(Y,x,X);
    
    ZZ = subs(Z,lambda,Lambda);
    YY = subs(Y,lambda,Lambda);
    XX = subs(X,lambda,Lambda);
    
    %Get rid of imaginary values
    closest_point = [double(XX),double(YY),double(ZZ)];
    closest_point = closest_point( imag(closest_point)==0);
    
    
    closest_point = reshape(closest_point,[],3);
    
    %Get the closest point on the bow shock
    [ClosestDistance,I] = min(vecnorm((closest_point - core_position'),2,2));
    closest_point = closest_point(I,:);
    
    
end