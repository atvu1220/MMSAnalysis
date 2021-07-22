%Calculates the Merka Bow Shock model from Merka et al. 2005 using Mach Numbers
function [shock_normal,HFAtoBS_Distance,bs_pos,closest_shock_normal,ClosestDistance,closest_point] = calculate_MerkaBS_CS(...
        mms1_fgm_timedata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
        left_InnerEdge,right_InnerEdge,...
        n_cs,n,v,M_A)
    %% Positions of Event
    
    %Interpolate the Data
    [~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw);
    [~,mms2_mec_rdata_interp] = interpxyz(mms2_mec_timedata_raw,mms2_mec_rdata_raw,mms1_fgm_timedata_raw);
    [~,mms3_mec_rdata_interp] = interpxyz(mms3_mec_timedata_raw,mms3_mec_rdata_raw,mms1_fgm_timedata_raw);
    [~,mms4_mec_rdata_interp] = interpxyz(mms4_mec_timedata_raw,mms4_mec_rdata_raw,mms1_fgm_timedata_raw);
    
    %Find center position of MMS within leading boundary
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
%     
%     leading_start_index = find(mms1_fgm_timedata_raw >= datenum(leading_start,formatIn), 1);
%     leading_end_index = find(mms1_fgm_timedata_raw >= datenum(leading_end,formatIn), 1);
%     
%     CenterofLeading_TetrahedronPosition = mean([...
%         mms1_mec_rdata_interp(leading_start_index:leading_end_index,:);...
%         mms2_mec_rdata_interp(leading_start_index:leading_end_index,:);...
%         mms3_mec_rdata_interp(leading_start_index:leading_end_index,:);...
%         mms4_mec_rdata_interp(leading_start_index:leading_end_index,:)...
%         ]);
    
    %Find center position of MMS within core of event
    left_InnerEdge_index = find(mms1_fgm_timedata_raw >= datenum(left_InnerEdge,formatIn), 1);
    right_InnerEdge_index = find(mms1_fgm_timedata_raw >= datenum(right_InnerEdge,formatIn), 1);
    
    CenterofCore_TetrahedronPosition = mean([...
        mms1_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms2_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms3_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:);...
        mms4_mec_rdata_interp(left_InnerEdge_index:right_InnerEdge_index,:)...
        ]);
%     
    %Find center position of MMS within trailing boundary
%     trailing_start_index = find(mms1_fgm_timedata_raw >= datenum(trailing_start,formatIn), 1);
%     trailing_end_index = find(mms1_fgm_timedata_raw >= datenum(trailing_end,formatIn), 1);
%     
%     CenterofTrailing_TetrahedronPosition = mean([...
%         mms1_mec_rdata_interp(trailing_start_index:trailing_end_index,:);...
%         mms2_mec_rdata_interp(trailing_start_index:trailing_end_index,:);...
%         mms3_mec_rdata_interp(trailing_start_index:trailing_end_index,:);...
%         mms4_mec_rdata_interp(trailing_start_index:trailing_end_index,:)...
%         ]);
%     
    
%     trailing_position=CenterofTrailing_TetrahedronPosition./6371.2;
%     leading_position=CenterofLeading_TetrahedronPosition./6371.2;
    core_position=CenterofCore_TetrahedronPosition./6371.2;
    %% Dynamic Pressure for Standoff Distance
    ScalingFactor = (n * norm(v)^2 / (7 * 457.5^2))^(1/6);
    
    %% Merka Bow Shock [2005]
    %Use fitted model from Merka for the Coefficients of B
    %     b1 = [0.0063, -0.0098];
    %     b3 = [0.8351,0.0102];
    %     b4 = [-0.0298,0.0040];
    %     b7 = [16.39, 0.2087,108.3];
    %     b8 = [-0.9241,0.0721];
    %     b10 = [-444,-2.935,-1930];
    %
    %
    %     a1  =  b1(1)  +  b1(2)*M_A;
    %     a3  =  b3(1)  +  b3(2)*M_A;
    %     a4  =  b4(1)  +  b4(2)*M_A;
    %     a7  =  b7(1)  +  b7(2)*M_A  + b7(3)/(M_A-1)^2;
    %     a8  =  b8(1)  +  b8(2)*M_A;
    %     a10 =  b10(1) +  b10(2)*M_A + b10(3)/(M_A-1)^2;
    
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
%     trailing_position=rotate_GSE2GPE(v,trailing_position);
    core_position = rotate_GSE2GPE(v,core_position);
%     leading_position = rotate_GSE2GPE(v,leading_position);
    
%     %Create two outward pointing vectors from event core to event edges.
%     trailing_core_vector = trailing_position - core_position;
%     leading_core_vector = leading_position - core_position;
%     
%     %The cross product of these two vectors gives the normal of the MMS Trajectory Plane
%     cross_trailing_leading = cross(trailing_core_vector,leading_core_vector);
%     
%     %The MMS Trajectory Plane Equation using the normal and its position at the event core.
%     spacecraft_trajectory_plane = ...
%         cross_trailing_leading(1).*(x-core_position(1))...
%         + cross_trailing_leading(2).*(y-core_position(2))...
%         + cross_trailing_leading(3).*(z-core_position(3));
%     
    %% Current Sheet Plane
    %Rotate the given current sheet normal to GPE coordinates
    n_cs = rotate_GSE2GPE(v,n_cs);
    
    %The Current Sheet PLane Equation using its normal and position at the event core.
    current_sheet_plane = ...
        n_cs(1).*(x-core_position(1))...
        + n_cs(2).*(y-core_position(2))...
        + n_cs(3).*(z-core_position(3));
    

    %% Calculate Bow Shock Surface and Current Sheet Intersection Curve
    syms x y z lambda
%     assume(x,'real')
%     assume(y,'real')
%     assume(z,'real')
%     
%     BCS_Z_xy = solve(merka_bowshock==current_sheet_plane,z);
%     Lagrange = (x-core_position(1))^2 + (y-core_position(2))^2 + (z-core_position(3))^2 + lambda.*(merka_bowshock);
%     Lagrange_xy = subs(Lagrange,z,BCS_Z_xy);
%     
%     dLdx = diff(Lagrange_xy,x);
%     dLdy = diff(Lagrange_xy,y);
%     
%     X = solve(dLdx==0,x);
%     Y = subs(dLdy,x,X);
%     Y = solve(Y==0,y);
%     X = subs(X,y,Y);
%     
%     Z = subs(BCS_Z_xy,x,X);
%     Z = subs(BCS_Z_xy,y,Y);
%     
%     
%             
    
    
    
    
    
    
    
    dLdx = 2*x - 2*core_position(1) + 2*lambda*a1*(ScalingFactor*x) + 2*lambda*a4*(ScalingFactor*y) + 2*lambda*a7 - n_cs(1);
    dLdy = 2*y - 2*core_position(2) + 2*lambda*(ScalingFactor*y) + 2*lambda*a4*(ScalingFactor*x) + 2*lambda*a8 - n_cs(2);
    dLdz = 2*z - 2*core_position(3) + 2*lambda*a3*(ScalingFactor*z) - n_cs(3);
    
    Z = solve(dLdz == 0,z);
    Y = solve(dLdy == 0,y);
    X = subs(dLdx,y,Y);
    X = solve(X==0,x);
    
    
    A = subs(merka_bowshock-current_sheet_plane,z,Z);
    A = subs(A,y,Y);
    A = subs(A,x,X);
    Lambda = solve(A==0,lambda);
    Y = subs(Y,x,X);
    
    ZZ = subs(Z,lambda,Lambda);
    YY = subs(Y,lambda,Lambda);
    XX = subs(X,lambda,Lambda);
    
    %Get rid of imaginary values
    bs_pos = [double(XX),double(YY),double(ZZ)];
%     bs_pos = bs_pos( imag(bs_pos)==0);

        bs_pos2 = [];
    for i=1:size(bs_pos,1)
        if imag(bs_pos(i,:)) == 0
            bs_pos2 = [bs_pos2; bs_pos(i,:)];
        end
    end
    
    bs_pos = bs_pos2;
    
%     bs_pos = reshape(bs_pos,[],3);
    
    %Get the closest point on the bow shock
    [HFAtoBS_Distance,I] = min(vecnorm((bs_pos - core_position'),2,2));
    bs_pos = bs_pos(I,:);
    
    
    
    
% % %      %% Calculate HFA Current Sheet Connection with BS in the plane of the MMS Trajectory   
% % %     X = solve(current_sheet_plane==spacecraft_trajectory_plane,x); %Solve f(x,y,z) for x. x=> x(y,z)
% % %     A = subs(current_sheet_plane,x,X); %Substitute x(y,z) back into the current sheet plane equation, f(x,y,z) = f(y,z)
% % %     Y = solve(A,y); %Solve for y(z)
% % %     
% % %     B = subs(merka_bowshock,x,X); %Substitute x(y,z) back into the bow shock surface equation
% % %     C = subs(B,y,Y); %Substitute y(z) back into the bow shock surface equation
% % %     ZZ = solve(C,z); %Solve for the last variable, z.
% % %     
% % %     YY = subs(Y,z,ZZ); %Substitute z back into y(z) for y
% % %     XX = subs(X,z,ZZ); %Substitute z back into x(y,z) for x(y)
% % %     
% % %     XX(1) = subs(XX(1),y,YY(1)); %Substitute y1 into x(y) for 1 solution
% % %     XX(2) = subs(XX(2),y,YY(2)); %Substitute y2 into x(y) for the 2nd solution
% % %     
% % %     x_bs = double(XX);
% % %     y_bs = double(YY);
% % %     z_bs = double(ZZ);
% % %     
% % %     bs_pos = [x_bs,y_bs,z_bs]; %Position on the bow shock that is the intersection of the current sheet and MMS trajectory planes
% % %     
% % %     %When there's two intersections, choose the intersection that is closer to the HFA observation position.
% % %     [HFAtoBS_Distance,I] = min(vecnorm((bs_pos - core_position'),2,2));
% % %     bs_pos = bs_pos(I,:);

%% Calculate Bow Shock Normal at Connection Point

% Calculation from Chu 2017 Thesis, p.62
%Gradient of the Merka Bow Shock at the Bow Shock Position gives its Shock Normal
    shock_normal = [2*a1*(ScalingFactor*bs_pos(1)) + 2*a4*(ScalingFactor*bs_pos(2)) + 2*(ScalingFactor*a7),...
        2*bs_pos(2) + 2*a4*(ScalingFactor*bs_pos(1)) + 2*(ScalingFactor*a8),...
        2*a3*(ScalingFactor*bs_pos(3))];
    %Normalize to unity
    shock_normal = shock_normal ./ norm(shock_normal);
    

    
    
    %% Calculate closest point from event observation to BS
    [ClosestDistance,closest_point] = calculate_closestDistance(merka_bowshock, core_position,a1,a3,a4,a7,a8,ScalingFactor);
    
    closest_shock_normal = [2*a1*closest_point(1) + 2*a4*closest_point(2) + 2*a7,...
        2*closest_point(2) + 2*a4*closest_point(1) + 2*a8,...
        2*a3*closest_point(3)];
    
    closest_shock_normal = closest_shock_normal ./ norm(closest_shock_normal);
    
    closest_point = rotate_GPE2GSE(v,closest_point);
    closest_point = closest_point';
    
    closest_shock_normal = rotate_GPE2GSE(v,closest_shock_normal);
    closest_shock_normal = closest_shock_normal';
    
    %Check if event is inside the bow shock model, due to the misapproximation of the bow shock model.
    insideOrOutsideBS = subs(merka_bowshock,[x y z],core_position');
    if insideOrOutsideBS < 0
        ClosestDistance = -ClosestDistance;
        HFAtoBS_Distance = -HFAtoBS_Distance;
    end
    
    %Rotate Back to GSE for plotting
    core_position = rotate_GPE2GSE(v,core_position);
    core_position = core_position';
    
% %     %% If the BS and CS don't conect, find the closest point of the CS with the BS
% %     if imag(bs_pos) ~= 0
% %         Intersection = -1;
% %         
% %         ZZ = solve(diff(Y,z)==diff(C,z));
% %         YY = subs(Y,z,ZZ);
% %         XX = subs(X,[y z],[YY ZZ]);
% %         
% %         CS_closest_position = double([XX;YY;ZZ]);
% %         
% %         %Solve for the minimum distance from the HFA to the BS with lagrange multipliers
% %         syms lambda
% %         dLdx = 2*x - 2*CS_closest_position(1) + 2*lambda*a1*x + 2*lambda*a4*y + 2*lambda*a7;
% %         dLdy = 2*y - 2*CS_closest_position(2) + 2*lambda*y + 2*lambda*a4*x + 2*lambda*a8;
% %         dLdz = 2*z - 2*CS_closest_position(3) + 2*lambda*a3*z;
% %         
% %         Z = solve(dLdz == 0,z);
% %         Y = solve(dLdy == 0,y);
% %         X = subs(dLdx,y,Y);
% %         X = solve(X==0,x);
% %         
% %         
% %         A = subs(merka_bowshock,z,Z);
% %         A = subs(A,y,Y);
% %         A = subs(A,x,X);
% %         Lambda = solve(A==0,lambda);
% %         Y = subs(Y,x,X);
% %         
% %         ZZ = subs(Z,lambda,Lambda);
% %         YY = subs(Y,lambda,Lambda);
% %         XX = subs(X,lambda,Lambda);
% %         
% %         %Get rid of imaginary values
% %         bs_pos = [double(XX),double(YY),double(ZZ)];
% %         bs_pos = bs_pos( imag(bs_pos)==0);
% %         
% %         bs_pos = reshape(bs_pos,[],3);
% %         
% %         %Get the closest point on the bow shock
% %         [HFAtoBS_Distance,I] = min(vecnorm((bs_pos - core_position'),2,2));
% %         bs_pos = bs_pos(I,:);
% %         
% %         % Calculation from Chu 2017 Thesis, p.62 for the bow shock normal at its position
% %         shock_normal = [2*a1*bs_pos(1) + 2*a4*bs_pos(2) + 2*a7,...
% %             2*bs_pos(2) + 2*a4*bs_pos(1) + 2*a8,...
% %             2*a3*bs_pos(3)];
% %         
% %         shock_normal = shock_normal ./ norm(shock_normal);
% %         
% %         shock_normal = rotate_GPE2GSE(v,shock_normal);
% %         
% %         bs_pos = rotate_GPE2GSE(v,bs_pos);
% %         
% %         bs_pos = bs_pos';
% %         shock_normal = shock_normal';
% %         
% %     else
% %         %If the CS does connect with the Bow Shock
% %         Intersection = 1;
% %     end
    
    
%% Plot Configuration
%Plot Bow shock Surface
Zmerka_bowshock = solve(merka_bowshock==0,z); %Solve for z
figure
fsurf(matlabFunction(Zmerka_bowshock(1)),'facecolor','none') %Plot the +z surface
hold on
fsurf(matlabFunction(Zmerka_bowshock(2)),'facecolor','none') %Plot the -z surface

%Plot the MMS Trajectory Plane
%     zspacecraft_trajectory_plane = solve(spacecraft_trajectory_plane==0,z);
%     fsurf(matlabFunction(zspacecraft_trajectory_plane))

%Plot the Current Sheet Plane
zcurrent_sheet_plane = solve(current_sheet_plane==0,z);
fsurf(matlabFunction(zcurrent_sheet_plane))


plot3(core_position(1),core_position(2),core_position(3),'.r','markersize',100)
plot3(bs_pos(1),bs_pos(2),bs_pos(3),'.b','markersize',150)

%Plot Properties
hold off
axis([-05 30 -30 30 -30 30])
view([45 15])
% view([0 90])
xlabel('x');
ylabel('y');
zlabel('z');
title('BS-CS Configuration')
plot_name =  strcat('BS-CS_Configuration');
%print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
print(gcf,'-dpng','-r300', '-loose', plot_name);
savefig(plot_name)

    %Rotate Back to GSE
    shock_normal = rotate_GPE2GSE(v,shock_normal);
    bs_pos = rotate_GPE2GSE(v,bs_pos);
    
    bs_pos = bs_pos';
    shock_normal = shock_normal';
    
end


function [ClosestDistance,closest_point] = calculate_closestDistance(merka_bowshock, core_position,a1,a3,a4,a7,a8,ScalingFactor)
    %% Calculate Closest Distance to BS
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