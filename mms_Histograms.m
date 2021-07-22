close all
clear;
tic
% %Shear Angle
% HFA_shear = get_data('HFA','Shear Angle');
% SHFA_shear = get_data('SHFA','Shear Angle');
% FB_shear = get_data('FB', 'Shear Angle');
% plot_histogram({'Magnetic Shear Angle','\theta'},10,HFA_shear,'HFA',SHFA_shear,'SHFA',FB_shear,'FB')

cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Histograms'
%% HFA, SHFA, FB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %X-Y
% plot_scatter({'X','Y'},...
%     get_data('HFA','Rx'),get_data('HFA','Ry'),'HFA',...
%     get_data('SHFA','Rx'),get_data('SHFA','Ry'),'SHFA',...
%     get_data('FB', 'Rx'), get_data('FB', 'Ry'),'FB')
%
% %X-Z
% plot_scatter({'X','Z'},...
%     get_data('HFA','Rx'),get_data('HFA','Rz'),'HFA',...
%     get_data('SHFA','Rx'),get_data('SHFA','Rz'),'SHFA',...
%     get_data('FB', 'Rx'),get_data('FB', 'Rz'),'FB')
%
% %Y-Z
% plot_scatter({'Y','Z'},...
%     get_data('HFA','Ry'),get_data('HFA','Rz'),'HFA',...
%     get_data('SHFA','Ry'),get_data('SHFA','Rz'),'SHFA',...
%     get_data('FB', 'Ry'),get_data('FB', 'Rz'),'FB')
%
%
% %CS-BS Angle
% plot_histogram({'CS-BS','[\Theta]'},10,...
%     get_data('HFA','CS-BS Angle'),'HFA',...
%     get_data('SHFA','CS-BS Angle'),'SHFA',...
%     get_data('FB', 'CS-BS Angle'),'FB')
%
% %Shear Angle
% plot_histogram({'Magnetic Shear Angle','[\theta]'},10,...
%     get_data('HFA','Shear Angle'),'HFA',...
%     get_data('SHFA','Shear Angle'),'SHFA',...
%     get_data('FB', 'Shear Angle'),'FB')
%
% %Shock Angle Before
% plot_histogram({'Shock Angle Before','[\theta]'},15,...
%     get_data('HFA','Shock Angle Before'),'HFA',...
%     get_data('SHFA','Shock Angle Before'),'SHFA',...
%     get_data('FB', 'Shock Angle Before'),'FB')
%
% %Shock Angle After
% plot_histogram({'Shock Angle After','[\theta]'},15,...
%     get_data('HFA','Shock Angle After'),'HFA',...
%     get_data('SHFA','Shock Angle After'),'SHFA',...
%     get_data('FB', 'Shock Angle After'),'FB')
%
%
%
% % %E-CS Before Angle
% % plot_histogram({'E-CS Angle Before','\theta'},15,...
% %     get_data('HFA','E-CS Angle Before'),'HFA',...
% %     get_data('SHFA','E-CS Angle Before'),'SHFA',...
% %     get_data('FB', 'E-CS Angle Before'),'FB')
% %
% % %E-CS After Angle
% % plot_histogram({'E-CS Angle After','\theta'},15,...
% %     get_data('HFA','E-CS Angle After'),'HFA',...
% %     get_data('SHFA','E-CS Angle After'),'SHFA',...
% %     get_data('FB', 'E-CS Angle After'),'FB')
%
%
% %% Xiao et al. 2015, expanding HFAs
% % n1_HFA = [get_data('HFA','Leading Timing Normal_x'), get_data('HFA','Leading Timing Normal_y'), get_data('HFA','Leading Timing Normal_z')];
% % n2_HFA = [get_data('HFA','Trailing Timing Normal_x'), get_data('HFA','Trailing Timing Normal_y'), get_data('HFA','Trailing Timing Normal_z')];
% % v1_HFA = get_data('HFA','Leading Timing Speed');
% % v2_HFA = get_data('HFA','Trailing Timing Speed');
% % preV_HFA = [get_data('HFA','preV_x'), get_data('HFA','preV_y'), get_data('HFA','preV_z')];
% % postV_HFA = [get_data('HFA','postV_x'), get_data('HFA','postV_y'), get_data('HFA','postV_z')];
% %
% % V1_HFA = v1_HFA - dot(preV_HFA,n1_HFA,2);
% % V2_HFA = v2_HFA - dot(postV_HFA,n2_HFA,2);
% %
% % Boundary_Movement_HFA = dot(V2_HFA.*n2_HFA,n1_HFA,2) - V1_HFA
% % %Boundary_Exp = vecnorm(Boundary_Movement,2,2);
% % %Boundary_Exp_HFA = v2_HFA-v1_HFA;
% % for i=1:length(n1_HFA)
% %     theta_n1_HFA_sunearth(i) = angle(n1_HFA(i,:),[1,0,0]);
% %     theta_n2_HFA_sunearth(i) = angle(n2_HFA(i,:),[1,0,0]);
% % end
% %
% %
% % n1_SHFA = [get_data('SHFA','Leading Timing Normal_x'), get_data('SHFA','Leading Timing Normal_y'), get_data('SHFA','Leading Timing Normal_z')];
% % n2_SHFA = [get_data('SHFA','Trailing Timing Normal_x'), get_data('SHFA','Trailing Timing Normal_y'), get_data('SHFA','Trailing Timing Normal_z')];
% % v1_SHFA = get_data('SHFA','Leading Timing Speed');
% % v2_SHFA = get_data('SHFA','Trailing Timing Speed');
% % preV_SHFA = [get_data('SHFA','preV_x'), get_data('SHFA','preV_y'), get_data('SHFA','preV_z')];
% % postV_SHFA = [get_data('SHFA','postV_x'), get_data('SHFA','postV_y'), get_data('SHFA','postV_z')];
% %
% % V1_SHFA = v1_SHFA - dot(preV_SHFA,n1_SHFA,2);
% % V2_SHFA = v2_SHFA - dot(postV_SHFA,n2_SHFA,2);
% %
% % Boundary_Movement_SHFA = dot(V2_SHFA.*n2_SHFA,n1_SHFA,2) - V1_SHFA
% % %Boundary_Movement_SHFA = v2_SHFA.*n2_SHFA-v1_SHFA.*n1_SHFA; %Needs to substract, V1 = U1(timing) - V_SW dot n1;
% % %Boundary_Exp = vecnorm(Boundary_Movement,2,2);
% % %Boundary_Exp_SHFA = v2_SHFA-v1_SHFA;
% % for i=1:length(n1_SHFA)
% %     theta_n1_SHFA_sunearth(i) = angle(n1_SHFA(i,:),[1,0,0]);
% %     theta_n2_SHFA_sunearth(i) = angle(n2_SHFA(i,:),[1,0,0]);
% % end
% %
% %
% %
% % plot_scatter({'V_{leading}';'V_{trailing}'},...
% %     v1_HFA,v2_HFA,'HFA',...
% %     v1_SHFA,v2_SHFA,'SHFA'); hold on
% % refline(1,0)
% % plot_scatter({'V_{leading}';'Delta V'},...
% %     v1_HFA,Boundary_Movement_HFA,'HFA',...
% %     v1_SHFA,Boundary_Movement_SHFA,'SHFA')
% % plot_scatter({'\theta_{leading}';'\theta_{trailing}'},...
% %     theta_n1_HFA_sunearth,theta_n2_HFA_sunearth,'HFA',...
% %     theta_n1_SHFA_sunearth,theta_n2_SHFA_sunearth,'SHFA'); hold on
% % refline(1,0)
% %
% %
%
%
% %% All Event Types
% %CS Speed
% plot_histogram({'CS Speed','[km/s]'},50,...
%     get_data('HFA','CS Speed'),'HFA',...
%     get_data('SHFA','CS Speed'),'SHFA',...
%     get_data('FB', 'CS Speed'),'FB')
%
% %Transversal Speed
% plot_histogram({'Transversal Speed','[km/s]'},50,...
%     get_data('HFA','Transversal Speed'),'HFA',...
%     get_data('SHFA','Transversal Speed'),'SHFA',...
%     get_data('FB', 'Transversal Speed'),'FB')
%
% %Gyration Ratio Before
% plot_histogram({'Gyration Ratio Before',''},0.1,...
%     get_data('HFA','Gyration Ratio Before'),'HFA',...
%     get_data('SHFA','Gyration Ratio Before'),'SHFA',...
%     get_data('FB', 'Gyration Ratio Before'),'FB')
%
% %Gyration Ratio After
% plot_histogram({'Gyration Ratio After',''},0.1,...
%     get_data('HFA','Gyration Ratio After'),'HFA',...
%     get_data('SHFA','Gyration Ratio After'),'SHFA',...
%     get_data('FB', 'Gyration Ratio After'),'FB')
%
% %Size
% plot_histogram({'Size','[Re]'},1,...
%     get_data('HFA','Size'),'HFA',...
%     get_data('SHFA','Size'),'SHFA',...
%     get_data('FB', 'Size'),'FB')
%
%
% %Expansion Speed
% plot_histogram({'Expansion Speed','[km/s]'},25,...
%     get_data('HFA','Expansion Speed'),'HFA',...
%     get_data('SHFA','Expansion Speed'),'SHFA',...
%     get_data('FB', 'Expansion Speed'),'FB')
%
% %Duration
% plot_histogram({'Duration','[s]'},5,...
%     get_data('HFA','Duration'),'HFA',...
%     get_data('SHFA','Duration'),'SHFA',...
%     get_data('FB', 'Duration'),'FB')
%
% %Age
% plot_histogram({'Age','[s]'},10,...
%     get_data('HFA','Age'),'HFA',...
%     get_data('SHFA','Age'),'SHFA',...
%     get_data('FB', 'Age'),'FB')
%
%
% %Distance Traveled
% plot_histogram({'Distance Traveled','[Re]'},1,...
%     get_data('HFA','Distance Traveled'),'HFA',...
%     get_data('SHFA','Distance Traveled'),'SHFA',...
%     get_data('FB', 'Distance Traveled'),'FB')
%
%
% %% HFA, SHFA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X-Y
plot_scatter({'X','Y'},...
    get_data('HFA','Rx'),get_data('HFA','Ry'),'HFA',...
    get_data('SHFA','Rx'),get_data('SHFA','Ry'),'SHFA')

%X-Z
plot_scatter({'X','Z'},...
    get_data('HFA','Rx'),get_data('HFA','Rz'),'HFA',...
    get_data('SHFA','Rx'),get_data('SHFA','Rz'),'SHFA')

%Y-Z
plot_scatter({'Y','Z'},...
    get_data('HFA','Ry'),get_data('HFA','Rz'),'HFA',...
    get_data('SHFA','Ry'),get_data('SHFA','Rz'),'SHFA')


%CS-BS Angle
plot_histogram({'CS-BS','[\Theta]'},10,...
    get_data('HFA','CS-BS Angle'),'HFA',...
    get_data('SHFA','CS-BS Angle'),'SHFA')

%Shear Angle
plot_histogram({'Magnetic Shear Angle','[\theta]'},10,...
    get_data('HFA','Shear Angle'),'HFA',...
    get_data('SHFA','Shear Angle'),'SHFA')

%Shock Angle Before
plot_histogram({'Shock Angle Before','[\theta]'},15,...
    get_data('HFA','Shock Angle Before'),'HFA',...
    get_data('SHFA','Shock Angle Before'),'SHFA')

%Shock Angle After
plot_histogram({'Shock Angle After','[\theta]'},15,...
    get_data('HFA','Shock Angle After'),'HFA',...
    get_data('SHFA','Shock Angle After'),'SHFA')



% %E-CS Before Angle
% plot_histogram({'E-CS Angle Before','\theta'},15,...
%     get_data('HFA','E-CS Angle Before'),'HFA',...
%     get_data('SHFA','E-CS Angle Before'),'SHFA',...
%     get_data('FB', 'E-CS Angle Before'),'FB')
%
% %E-CS After Angle
% plot_histogram({'E-CS Angle After','\theta'},15,...
%     get_data('HFA','E-CS Angle After'),'HFA',...
%     get_data('SHFA','E-CS Angle After'),'SHFA',...
%     get_data('FB', 'E-CS Angle After'),'FB')



%CS Speed
plot_histogram({'CS Speed','[km/s]'},50,...
    get_data('HFA','CS Speed'),'HFA',...
    get_data('SHFA','CS Speed'),'SHFA')

%Transversal Speed
plot_histogram({'Transversal Speed','[km/s]'},50,...
    get_data('HFA','Transversal Speed'),'HFA',...
    get_data('SHFA','Transversal Speed'),'SHFA')

%Gyration Ratio Before
plot_histogram({'Gyration Ratio Before',''},0.1,...
    get_data('HFA','Gyration Ratio Before'),'HFA',...
    get_data('SHFA','Gyration Ratio Before'),'SHFA')

%Gyration Ratio After
plot_histogram({'Gyration Ratio After',''},0.1,...
    get_data('HFA','Gyration Ratio After'),'HFA',...
    get_data('SHFA','Gyration Ratio After'),'SHFA')

%Size
plot_histogram({'Size','[Re]'},1,...
    get_data('HFA','Size'),'HFA',...
    get_data('SHFA','Size'),'SHFA')


%Expansion Speed
plot_histogram({'Expansion Speed','[km/s]'},25,...
    get_data('HFA','Expansion Speed'),'HFA',...
    get_data('SHFA','Expansion Speed'),'SHFA')

%Duration
plot_histogram({'Duration','[s]'},5,...
    get_data('HFA','Duration'),'HFA',...
    get_data('SHFA','Duration'),'SHFA')

%Age
plot_histogram({'Age','[s]'},10,...
    get_data('HFA','Age'),'HFA',...
    get_data('SHFA','Age'),'SHFA')


%% Substructures (SS) & No Substructures (NS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X-Y
plot_scatter({'X','Y'},...
    get_data('SS','Rx'),get_data('SS','Ry'),'Substructures',...
    get_data('NS','Rx'),get_data('NS','Ry'),'No Substructures')

%X-Z
plot_scatter({'X','Z'},...
    get_data('SS','Rx'),get_data('SS','Rz'),'Substructures',...
    get_data('NS','Rx'),get_data('NS','Rz'),'No Substructures')

%Y-Z
plot_scatter({'Y','Z'},...
    get_data('SS','Ry'),get_data('SS','Rz'),'Substructures',...
    get_data('NS','Ry'),get_data('NS','Rz'),'No Substructures')


%CS-BS Angle
plot_histogram({'CS-BS','[\Theta]'},10,...
    get_data('SS','CS-BS Angle'),'Substructures',...
    get_data('NS','CS-BS Angle'),'No Substructures')

%Shear Angle
plot_histogram({'Magnetic Shear Angle','[\theta]'},10,...
    get_data('SS','Shear Angle'),'Substructures',...
    get_data('NS','Shear Angle'),'No Substructures')

%Shock Angle Before
plot_histogram({'Shock Angle Before','[\theta]'},15,...
    get_data('SS','Shock Angle Before'),'Substructures',...
    get_data('NS','Shock Angle Before'),'No Substructures')

%Shock Angle After
plot_histogram({'Shock Angle After','[\theta]'},15,...
    get_data('SS','Shock Angle After'),'Substructures',...
    get_data('NS','Shock Angle After'),'No Substructures')



% %E-CS Before Angle
% plot_histogram({'E-CS Angle Before','\theta'},15,...
%     get_data('SS','E-CS Angle Before'),'Substructures',...
%     get_data('NS','E-CS Angle Before'),'No Substructures')
%
% %E-CS After Angle
% plot_histogram({'E-CS Angle After','\theta'},15,...
%     get_data('SS','E-CS Angle After'),'Substructures',...
%     get_data('NS','E-CS Angle After'),'No Substructures')



%CS Speed
plot_histogram({'CS Speed','[km/s]'},50,...
    get_data('SS','CS Speed'),'Substructures',...
    get_data('NS','CS Speed'),'No Substructures')

%Transversal Speed
plot_histogram({'Transversal Speed','[km/s]'},50,...
    get_data('SS','Transversal Speed'),'Substructures',...
    get_data('NS','Transversal Speed'),'No Substructures')

%Gyration Ratio Before
plot_histogram({'Gyration Ratio Before',''},0.1,...
    get_data('SS','Gyration Ratio Before'),'Substructures',...
    get_data('NS','Gyration Ratio Before'),'No Substructures')

%Gyration Ratio After
plot_histogram({'Gyration Ratio After',''},0.1,...
    get_data('SS','Gyration Ratio After'),'Substructures',...
    get_data('NS','Gyration Ratio After'),'No Substructures')

%Size
plot_histogram({'Size','[Re]'},1,...
    get_data('SS','Size'),'Substructures',...
    get_data('NS','Size'),'No Substructures')


%Expansion Speed
plot_histogram({'Expansion Speed','[km/s]'},25,...
    get_data('SS','Expansion Speed'),'Substructures',...
    get_data('NS','Expansion Speed'),'No Substructures')

%Duration
plot_histogram({'Duration','[s]'},5,...
    get_data('SS','Duration'),'Substructures',...
    get_data('NS','Duration'),'No Substructures')

%Age
plot_histogram({'Age','[s]'},10,...
    get_data('SS','Age'),'Substructures',...
    get_data('NS','Age'),'No Substructures')


%Distance Traveled
plot_histogram({'Distance Traveled','[Re]'},1,...
    get_data('SS','Distance Traveled'),'Substructures',...
    get_data('NS','Distance Traveled'),'No Substructures')

%Distance Traveled
plot_histogram({'Distance Traveled','[Re]'},1,...
    get_data('SS','Distance Traveled'),'HFA',...
    get_data('NS','Distance Traveled'),'SHFA')







%Leading and Trailing Velocities (SW Frame)
plot_scatter({'V_{leading}','V_{trailing}'},...
    get_data('SS','V1'),get_data('SS','V2'),'Substructures',...
    get_data('NS','V1'),get_data('NS','V2'),'No Substructures')

%Leading and Trailing Velocities (Not SW Frame)
plot_scatter({'Vt_{leading}','Vt_{trailing}'},...
    get_data('SS','Vt1'),get_data('SS','Vt2'),'Substructures',...
    get_data('NS','Vt1'),get_data('NS','Vt2'),'No Substructures')

%Angles of Leading and Trailing Boundaries with Sun-Earth
plot_scatter({'Theta_{leading}','Theta_{trailing}'},...
    get_data('SS','N1'),get_data('SS','N2'),'Substructures',...
    get_data('NS','N1'),get_data('NS','N2'),'No Substructures')

%Angles of Leading and Trailing Boundaries with Sun-Earth, difference
plot_scatter({'Magnetic Shear','Theta_{diff}'},...
    get_data('SS','Shear Angle'),get_data('SS','N2')-get_data('SS','N1'),'Substructures',...
    get_data('NS','Shear Angle'),get_data('NS','N2')-get_data('NS','N1'),'No Substructures')

%Magnetic Shear vs delta_V (SW Frame)
plot_scatter({'Magnetic Shear','V_{diff}'},...
    get_data('SS','Shear Angle'),get_data('SS','BoundaryExpansion'),'Substructures',...
    get_data('NS','Shear Angle'),get_data('NS','BoundaryExpansion'),'No Substructures')

%Magnetic Shear vs delta_V (Not SW Frame)
plot_scatter({'Magnetic Shear','Vt_{diff}'},...
    get_data('SS','Shear Angle'),get_data('SS','BoundaryExpansion2'),'Substructures',...
    get_data('NS','Shear Angle'),get_data('NS','BoundaryExpansion2'),'No Substructures')



%CS angle (with sun Earth) vs delta_V
currentSheetNormal_SS = [get_data('SS','CS_x'),get_data('SS','CS_y'),get_data('SS','CS_z')];
for i=1:length(currentSheetNormal_SS)
    CSSE_Angle_SS(i) = angle(currentSheetNormal_SS(i,:),[1,0,0]);
end

currentSheetNormal_NS = [get_data('NS','CS_x'),get_data('NS','CS_y'),get_data('NS','CS_z')];
for i=1:length(currentSheetNormal_NS)
    CSSE_Angle_NS(i) = angle(currentSheetNormal_NS(i,:),[1,0,0]);
end

%CS angle vs delta_V (SW Frame)
plot_scatter({'CS Angle','V_{diff}'},...
    CSSE_Angle_SS,get_data('SS','BoundaryExpansion2'),'Substructures',...
    CSSE_Angle_NS,get_data('NS','BoundaryExpansion2'),'No Substructures')


%CS-BS angle with delta_v (sw frame)
plot_scatter({'CS-BS Angle','V_{diff}'},...
    get_data('SS','CS-BS Angle'),get_data('SS','BoundaryExpansion2'),'Substructures',...
    get_data('NS','CS-BS Angle'),get_data('NS','BoundaryExpansion2'),'No Substructures')


%Event Size with delta_v (sw frame)
plot_scatter({'Event Size','V_{diff}'},...
    get_data('SS','Size'),get_data('SS','BoundaryExpansion2'),'Substructures',...
    get_data('NS','Size'),get_data('NS','BoundaryExpansion2'),'No Substructures')


%Boundary Edges Angle with Delta_v (swFrame)
plot_scatter({'Edges Angle','V_{diff}'},...
    180-get_data('SS','Event Edges Angle'),get_data('SS','BoundaryExpansion2'),'Substructures',...
    180-get_data('NS','Event Edges Angle'),get_data('NS','BoundaryExpansion2'),'No Substructures')

currentSheetNormal_NS = [get_data('NS','CS_x'),get_data('NS','CS_y'),get_data('NS','CS_z')];
%BS Normal with Delta_v (swFrame)
plot_scatter({'Edges Angle','V_{diff}'},...
    180-get_data('SS','Event Edges Angle'),get_data('SS','BoundaryExpansion2'),'Substructures',...
    180-get_data('NS','Event Edges Angle'),get_data('NS','BoundaryExpansion2'),'No Substructures')



%Vdiff histograms
plot_histogram({'V_{diff}','[km/s]'},10,...
    get_data('SS','BoundaryExpansion'),'Substructure',get_data('NS','BoundaryExpansion'),'No Substructure')
plot_histogram({'Vt_{diff}','[km/s]'},10,...
    get_data('SS','BoundaryExpansion2'),'Substructure',get_data('NS','BoundaryExpansion2'),'No Substructure')



%Core Densities Histogram
plot_histogram({'Core Density STD','[\sigma]'},0.2,...
    [get_data('SS','Core Density STD');get_data('NS','Core Density STD')],'HFA')

%Core Density STD with delta_v (swframe)
plot_scatter({'Core Density STD','V_{diff}'},...
    [get_data('SS','Core Density STD');get_data('NS','Core Density STD')],[get_data('SS','BoundaryExpansion');get_data('NS','BoundaryExpansion')],'HFA')

%Core Density STD with shear
plot_scatter({'Core Density STD','Shear Angle'},...
    [get_data('SS','Core Density STD');get_data('NS','Core Density STD')],[get_data('SS','Shear Angle');get_data('NS','Shear Angle')],'HFA')

%Core Density STD with CS-BS Angle
plot_scatter({'Core Density STD','CS-BS Angle'},...
    [get_data('SS','Core Density STD');get_data('NS','Core Density STD')],[get_data('SS','CS-BS Angle');get_data('NS','CS-BS Angle')],'HFA')








%Event Edges shock angle
plot_scatter({'Leading Shock Angle','Trailing Shock Angle'},...
    get_data('SS','Leading Shock Angle'),get_data('SS','Trailing Shock Angle'),'Substructures',...
    get_data('NS','Leading Shock Angle'),get_data('NS','Trailing Shock Angle'),'No Substructures')


%Event Shock Edges Angle with B
plot_histogram({'Leading Shock Angle','[\theta]'},15,...
    [get_data('SS','Leading Shock Angle');get_data('NS','Leading Shock Angle')],'HFA')

plot_histogram({'Leading Shock Angle','[\theta]'},15,...
    get_data('SS','Leading Shock Angle'),'Substructure',get_data('NS','Leading Shock Angle'),'No Substructure')

plot_histogram({'Trailing Shock Angle','[\theta]'},15,...
    [get_data('SS','Trailing Shock Angle');get_data('NS','Trailing Shock Angle')],'HFA')

plot_histogram({'Trailing Shock Angle','[\theta]'},15,...
    get_data('SS','Trailing Shock Angle'),'Substructure',get_data('NS','Trailing Shock Angle'),'No Substructure')





XY_mag_SS = (get_data('SS','Rx').^2+get_data('SS','Ry').^2).^(1/2);
XY_mag_NS = (get_data('NS','Rx').^2+get_data('NS','Ry').^2).^(1/2);

%XY distance vs delta_v
plot_scatter({'XY Distance','V_{diff}'},...
    XY_mag_SS,get_data('SS','BoundaryExpansion2'),'Substructures',...
    XY_mag_NS,get_data('NS','BoundaryExpansion2'),'No Substructures')

%XY distance vs delta_v
plot_scatter({'Core Density STD','XY Distance'},...
    [get_data('SS','Core Density STD');get_data('NS','Core Density STD')],[XY_mag_SS;XY_mag_NS],'HFA')


%E-CS Angles
plot_scatter({'E-CS Angle Before','E-CS Angle After'},...
    get_data('SS','E-CS Angle Before'),get_data('SS','E-CS Angle After'),'Substructures',...
    get_data('NS','E-CS Angle Before'),get_data('NS','E-CS Angle After'),'No Substructures')

plot_histogram({'E-CS Angle Before','[\theta]'},15,...
    get_data('SS','E-CS Angle Before'),'Substructure',get_data('NS','E-CS Angle Before'),'No Substructure')
plot_histogram({'E-CS Angle After','[\theta]'},15,...
    get_data('SS','E-CS Angle After'),'Substructure',get_data('NS','E-CS Angle After'),'No Substructure')


% vs core density STD
plot_scatter({'Core Density STD','Core Temp Corr'},...
    [get_data('SS','Core Density STD');get_data('NS','Core Density STD')],[get_data('SS','Core Temp Corr');get_data('NS','Core Temp Corr')],'HFA')

% vs core density CV
plot_scatter({'Core Density CV','Core Temp Corr'},...
    [get_data('SS','Core Density CV');get_data('NS','Core Density CV')],[get_data('SS','Core Temp Corr');get_data('NS','Core Temp Corr')],'HFA')


%f vs core density STD
plot_histogram({'Core Density STD'},0.5,...
    get_data('SS','Core Density STD'),'Substructure',get_data('NS','Core Density STD'),'No Substructure')

%f vs core density CV
plot_histogram({'Core Density CV'},0.05,...
    get_data('SS','Core Density CV'),'Substructure',get_data('NS','Core Density CV'),'No Substructure')

%f vs core density CV
plot_histogram({'Core Density CV'},0.05,...
    get_data('HFA','Core Density CV'),'HFA',...
    get_data('SHFA','Core Density CV'),'SHFA')

% %XY distance vs age
% plot_scatter({'XY Distance','Age'},...
%     XY_mag_SS,get_data('SS','Age'),'Substructures',...
%     XY_mag_NS,get_data('NS','Age'),'No Substructures')
%
% %XY distance vs size
% plot_scatter({'XY Distance','Size'},...
%     XY_mag_SS,get_data('SS','Size'),'Substructures',...
%     XY_mag_NS,get_data('NS','Size'),'No Substructures')
%
% %XY distance vs Distance Traveled
% plot_scatter({'XY Distance','Distance Traveled'},...
%     XY_mag_SS,get_data('SS','Distance Traveled'),'Substructures',...
%     XY_mag_NS,get_data('NS','Distance Traveled'),'No Substructures')



%Core Density STD with delta_V (swframe)




% % %Event Shock Angle
% plot_scatter({'Theta_{leading}','Theta_{trailing}'},...
%     get_data('SS','N1'),get_data('SS','N2'),'Substructures',...
%     get_data('NS','N1'),get_data('NS','N2'),'No Substructures')
%
% B_pre_SS = [get_data('SS','Bpre_x'),get_data('SS','Bpre_y'),get_data('SS','Bpre_z')];
% B_post_SS = [get_data('SS','Bpost_x'),get_data('SS','Bpost_y'),get_data('SS','Bpost_z')];
%

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%Plot the histogram, up to 3 data points
function [] = plot_histogram(parameter,bin_size,data1,label1,data2,label2,data3,label3)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 350 350])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    if min(data1) < 0
        min_bin = min(data1)-bin_size;
    else
        min_bin = 0;
    end
    
    if nargin == 4
        histogram(data1,min_bin:bin_size:max(data1)+bin_size,'FaceColor','r','FaceAlpha',0.25,'LineWidth',2,'Displayname',label1)
        fileName = strcat(cell2mat(parameter(1)),'_',label1);
    elseif nargin == 6
        histogram(data1,min_bin:bin_size:max(data1)+bin_size,'FaceColor','r','FaceAlpha',0.25,'LineWidth',2,'Displayname',label1)
        hold on
        histogram(data2,min_bin:bin_size:max(data2)+bin_size,'FaceColor','y','FaceAlpha',0.25,'LineWidth',2,'Displayname',label2)
        fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2);
    elseif nargin == 8
        histogram(data1,min_bin:bin_size:max(data1)+bin_size,'FaceColor','r','FaceAlpha',0.25,'LineWidth',2,'Displayname',label1)
        hold on
        histogram(data2,min_bin:bin_size:max(data2)+bin_size,'FaceColor','y','FaceAlpha',0.25,'LineWidth',2,'Displayname',label2)
        histogram(data3,min_bin:bin_size:max(data3)+bin_size,'FaceColor','c','FaceAlpha',0.25,'LineWidth',2,'Displayname',label3)
        fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2,'_',label3);
    end
    
    colormap('winter');
    xlabel(parameter,'FontSize',14)
    ylabel({'# of Events'},'FontSize',14)
    title(strcat('f vs.', {' '}, parameter(1)),'FontSize',16,'FontWeight', 'normal')
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    
    legend
    
    print(gcf,'-dpng','-r300', '-loose', strcat(fileName));
    
    
end

%Plot the scatter, up to 3 data points
function [] = plot_scatter(parameter,data1a,data1b,label1,data2a,data2b,label2,data3a,data3b,label3)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 600 400])
    co = [0 0 1;
        0 1 0;
        1 0 0;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    if nargin == 4
        scatter(data1a,data1b,'filled','Displayname',label1)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1);
    elseif nargin == 7
        scatter(data1a,data1b,'filled','Displayname',label1)
        hold on
        scatter(data2a,data2b,'filled','Displayname',label2)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1,'_',label2);
    elseif nargin == 10
        scatter(data1a,data1b,'filled','Displayname',label1)
        hold on
        scatter(data2a,data2b,'filled','Displayname',label2)
        scatter(data3a,data3b,'filled','Displayname',label3)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1,'_',label2,'_',label3);
    end
    
    xlabel(parameter(1),'FontSize',14)
    ylabel(parameter(2),'FontSize',14)
    %title(strcat('Observation Position', { ' ' }, parameter(1),'-',parameter(2)),'FontSize',16,'FontWeight', 'normal')
    title(strcat(parameter(1), { ' ' },'vs.',{ ' ' },parameter(2)),'FontSize',16,'FontWeight', 'normal')
    set(gca,'XMinorTick','on','TickDir','out','YMinorTick','on','linewidth',2)
    legend('Location','SouthEast')
    grid on
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    box on
    
    if contains(parameter(1),'V') && contains(parameter(2),'V')
        refline(1,0);
    end
    print(gcf,'-dpng','-r300', '-loose', strcat(fileName));
    
end

%Get the data from the Spreadsheet
function [data] = get_data(Event_Type,Event_Parameter)
    Page_number = 1;
    Column_number = 'A1';
    %     Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabaseFebruary26.xlsx';
    Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabase.xlsx';
    
    switch Event_Type
        case 'SS'
            Page_number = 2;
        case 'NS'
            Page_number = 3;
        case 'HFA'
            Page_number = 4;
        case 'SHFA'
            Page_number = 5;
        case 'FB'
            Page_number = 6;
    end
    
    
    switch Event_Parameter
        case 'Rx'
            Column_number = 'F';
        case 'Ry'
            Column_number = 'G';
        case 'Rz'
            Column_number = 'H';
            
        case 'CS_x'
            Column_number = 'I';
        case 'CS_y'
            Column_number = 'J';
        case 'CS_z'
            Column_number = 'K';
            
            
        case 'BS_x'
            Column_number = 'L';
        case 'BS_y'
            Column_number = 'M';
        case 'BS_z'
            Column_number = 'N';
            
            
        case 'CS-BS Angle'
            Column_number = 'O';
        case 'Shear Angle'
            Column_number = 'P';
            
            
        case 'Leading Timing Normal_x'
            Column_number = 'V';
        case 'Leading Timing Normal_y'
            Column_number = 'W';
        case 'Leading Timing Normal_z'
            Column_number = 'X';
        case 'Leading Timing Speed'
            Column_number = 'Y';
            
            
            
            
        case 'Trailing Timing Normal_x'
            Column_number = 'AF';
        case 'Trailing Timing Normal_y'
            Column_number = 'AG';
        case 'Trailing Timing Normal_z'
            Column_number = 'AH';
        case 'Trailing Timing Speed'
            Column_number = 'AI';
            
            
        case 'Bpre_x'
            Column_number = 'AK';
        case 'Bpre_y'
            Column_number = 'AL';
        case 'Bpre_z'
            Column_number = 'AM';
        case 'Bpost_x'
            Column_number = 'AN';
        case 'Bpost_y'
            Column_number = 'AO';
        case 'Bpost_z'
            Column_number = 'AP';
            
            
            
        case 'preV_x'
            Column_number = 'AR';
        case 'preV_y'
            Column_number = 'AS';
        case 'preV_z'
            Column_number = 'AT';
            
        case 'postV_x'
            Column_number = 'AU';
        case 'postV_y'
            Column_number = 'AV';
        case 'postV_z'
            Column_number = 'AW';
            
            
            
            
        case 'Shock Angle Before'
            Column_number = 'AX';
        case 'Shock Angle After'
            Column_number = 'AY';
        case 'E-CS Angle Before'
            Column_number = 'AZ';
        case 'E-CS Angle After'
            Column_number = 'BA';
        case 'CS Speed'
            Column_number = 'BB';
        case 'Transversal Speed'
            Column_number = 'BF';
        case 'Gyration Ratio Before'
            Column_number = 'BH';
        case 'Gyration Ratio After'
            Column_number = 'BJ';
        case 'Size'
            Column_number = 'BT';
        case 'Expansion Speed'
            Column_number = 'BU';
        case 'Duration'
            Column_number = 'BV';
        case 'Core Duration'
            Column_number = 'BW';
        case 'Age'
            Column_number = 'BX';
        case 'Distance Traveled'
            Column_number = 'BY';
        case 'V1'
            Column_number = 'BZ';
        case 'V2'
            Column_number = 'CA';
        case 'BoundaryExpansion'
            Column_number = 'CB';
            
            
            %         case 'BoundaryExpansion2'
            %             Column_number = 'CC';
            %         case 'N1'
            %             Column_number = 'CD';
            %         case 'N2'
            %             Column_number = 'CE';
            
            
        case 'Vt1'
            Column_number = 'CC';
        case 'Vt2'
            Column_number = 'CD';
        case 'BoundaryExpansion2'
            Column_number = 'CE';
        case 'N1'
            Column_number = 'CF';
        case 'N2'
            Column_number = 'CG';
            
        case 'Leading Shock Angle'
            Column_number = 'CH';
        case 'Trailing Shock Angle'
            Column_number = 'CI';
        case 'Event Edges Angle'
            Column_number = 'CJ';
        case 'Core Density STD'
            Column_number = 'CK';
            
        case 'Core Temp Corr'
            Column_number = 'CL';
            
        case 'Core Density CV'
            Column_number = 'CM';
            
            
            
            
    end
    
    Column_range = strcat(Column_number,'2:',Column_number,'200');
    
    data = xlsread(Database_Directory,Page_number,Column_range);
    data = data(~isnan(data));
    
end
