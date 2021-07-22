close all
clear;
tic
% %Shear Angle
% HFA_shear = get_data('HFA','Shear Angle');
% SHFA_shear = get_data('SHFA','Shear Angle');
% FB_shear = get_data('FB', 'Shear Angle');
% plot_histogram({'Magnetic Shear Angle','\theta'},10,HFA_shear,'HFA',SHFA_shear,'SHFA',FB_shear,'FB')

cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Histograms'

%Strings for types Of events
type1 = 'SS';
name1 = 'Substructures (2 sigma)';
type2 = 'NS';
name2 = 'No Substructures';

% type1 = 'HFA';
% name1 = 'HFA';
% type2 = 'SHFA';
% name2 = 'SHFA';
HFAnumber = size(get_data(type1,'Rx'))
SHFAnumber = size(get_data(type2,'Rx'))

% %% Substructures (SS) & No Substructures (NS)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %X-Y
% plot_scatter({'X','Y'},...
%     get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
%     get_data(type2,'Rx'),get_data(type2,'Ry'),name2)
% 
% %X-Z
% plot_scatter({'X','Z'},...
%     get_data(type1,'Rx'),get_data(type1,'Rz'),name1,...
%     get_data(type2,'Rx'),get_data(type2,'Rz'),name2)
% 
% %Y-Z
% plot_scatter({'Y','Z'},...
%     get_data(type1,'Ry'),get_data(type1,'Rz'),name1,...
%     get_data(type2,'Ry'),get_data(type2,'Rz'),name2)
% 
% 
% % %CS-BS Angle
% plot_histogram({'CS-BS','[\Theta]'},10,...
%     get_data(type1,'CS-BS Angle'),name1,...
%     get_data(type2,'CS-BS Angle'),name2)
% 
% %CS-BS Angle
% % plot_histogram({'CS-BS','[\Theta]'},10,...
% %     get_data('HFA','CS-BS Angle'),'HFA')
% 
% %Shear Angle
% plot_histogram({'Magnetic Shear Angle','[\theta]'},15,...
%     get_data(type1,'Shear Angle'),name1,...
%     get_data(type2,'Shear Angle'),name2)
% 
%Shock Angle Before
plot_histogram({'Shock Angle Downstream','[\theta]'},15,...
    get_data(type1,'Shock Angle Before'),name1,...
    get_data(type2,'Shock Angle Before'),name2)

% plot_scatter({'Shock Angle Before','Duration'},...
%     get_data(type1,'Shock Angle Before'),get_data(type1,'Duration'),name1,...
%     get_data(type2,'Shock Angle Before'),get_data(type2,'Duration'),name2)

plot_scatter({'Shock Angle Downstream','Size'},...
    get_data(type1,'Shock Angle Before'),get_data(type1,'Size'),name1,...
    get_data(type2,'Shock Angle Before'),get_data(type2,'Size'),name2)

%Shock Angle After
plot_histogram({'Shock Angle Upstream','[\theta]'},15,...
    get_data(type1,'Shock Angle After'),name1,...
    get_data(type2,'Shock Angle After'),name2)

plot_scatter({'Shock Angle Upstream','Duration'},...
    get_data(type1,'Shock Angle After'),get_data(type1,'Duration'),name1,...
    get_data(type2,'Shock Angle After'),get_data(type2,'Duration'),name2)

plot_scatter({'Shock Angle Upstream','Size'},...
    get_data(type1,'Shock Angle After'),get_data(type1,'Size'),name1,...
    get_data(type2,'Shock Angle After'),get_data(type2,'Size'),name2)


plot_scatter({'Shock Angle Downstream','Shock Angle Upstream'},...
    get_data(type1,'Shock Angle Before'),get_data(type1,'Shock Angle After'),name1,...
    get_data(type2,'Shock Angle Before'),get_data(type2,'Shock Angle After'),name2)

% 
% 
% 

%E Fields
%E-CS Before Mag
ECSMag_type1 = get_data(type1,'E-CS Magnitude Before');
ECSAngle_type1 = get_data(type1,'E-CS Angle Before');
ECSAngle_type1 = ECSAngle_type1(abs(ECSMag_type1) > 200);
ECSMag_type1 = ECSMag_type1(abs(ECSMag_type1) <= 4000);

ECSMag_type2 = get_data(type2,'E-CS Magnitude Before');
ECSAngle_type2 = get_data(type2,'E-CS Angle Before');
ECSAngle_type2 = ECSAngle_type2(abs(ECSMag_type2) > 200);
ECSMag_type2 = ECSMag_type2(abs(ECSMag_type2) <= 4000);

plot_histogram({'E-CS Magnitude Downstream','[\theta]'},250,...
    ECSMag_type1,name1,...
    ECSMag_type2,name2)

%E-CS Before Angle
plot_histogram({'E-CS Angle Before','\theta'},22.5,...
    get_data(type1,'E-CS Angle Before'),name1,...
    get_data(type2,'E-CS Angle Before'),name2)


%E-CS After Mag
ECSMag_type1 = get_data(type1,'E-CS Magnitude After');
ECSAngle_type1 = get_data(type1,'E-CS Angle After');
ECSAngle_type1 = ECSAngle_type1(abs(ECSMag_type1) >= 150);
ECSMag_type1 = ECSMag_type1(abs(ECSMag_type1) <= 4000);

ECSMag_type2 = get_data(type2,'E-CS Magnitude After');
ECSAngle_type2 = get_data(type2,'E-CS Angle After');
ECSAngle_type2 = ECSAngle_type2(abs(ECSMag_type2) >= 150);
ECSMag_type2 = ECSMag_type2(abs(ECSMag_type2) <= 4000);

plot_histogram({'E-CS Magnitude Upstream','[\theta]'},250,...
    ECSMag_type1,name1,...
    ECSMag_type2,name2)

%E-CS After Angle
plot_histogram({'E-CS Angle After','\theta'},22.5,...
    ECSAngle_type1,name1,...
    ECSAngle_type2,name2)









% 
% % %E-CS Before Angle
% % plot_histogram({'E-CS Angle Before','\theta'},15,...
% %     get_data('HFA','E-CS Angle Before'),'HFA')
% % 
% % %E-CS After Angle
% % plot_histogram({'E-CS Angle After','\theta'},15,...
% %     get_data('HFA','E-CS Angle After'),'HFA')
% 
% 
% %E-CS Angles
% plot_scatter({'E-CS Angle Before','E-CS Angle After'},...
%     get_data(type1,'E-CS Angle Before'),get_data(type1,'E-CS Angle After'),name1,...
%     get_data(type2,'E-CS Angle Before'),get_data(type2,'E-CS Angle After'),name2)
% 
% % plot_scatter({'E-CS Angle Before','E-CS Angle After'},...
% %     get_data('HFA','E-CS Angle Before'),get_data('HFA','E-CS Angle After'),'HFA')
% 
% %Combined
% plot_histogram({'E-CS Angle Downstream','[\theta]'},30,...
%     get_data(type1,'E-CS Angle Before'),name1,get_data(type2,'E-CS Angle Before'),name2)
% 
% plot_histogram({'E-CS Angle Upstream','[\theta]'},30,...
%     get_data(type1,'E-CS Angle After'),name1,get_data(type2,'E-CS Angle After'),name2)
% 
% 
% 
% %CS Speed
% plot_histogram({'CS Speed','[km/s]'},50,...
%     get_data(type1,'CS Speed'),name1,...
%     get_data(type2,'CS Speed'),name2)
% plot_histogram({'CS Speed','[km/s]'},50,...
%     get_data('HFA','CS Speed'),'HFA')


% %Transversal Speed
% plot_histogram({'Transversal Speed','[km/s]'},50,...
%     get_data(type1,'Transversal Speed'),name1,...
%     get_data(type2,'Transversal Speed'),name2)



%             
%         case 'Cone Angle Before'
%             Column_number = 'BB';
%         case 'Cone Angle After'
%             Column_number = 'BC';
%         case 'Alfven V Before'
%             Column_number = 'BD';
%         case 'Mach Number Before'
%             Column_number = 'BE';
%         case 'Alfven V After'
%             Column_number = 'BF';
%         case 'Mach umber After'
%             Column_number = 'BG';


%Cone Angle
plot_histogram({'Cone Angle Downstream','[\theta]'},15,...
    get_data(type1,'Cone Angle Before'),name1,...
    get_data(type2,'Cone Angle Before'),name2)

plot_scatter({'Cone Angle Downstream','Size'},...
    get_data(type1,'Cone Angle Before'),get_data(type1,'Size'),name1,...
    get_data(type2,'Cone Angle Before'),get_data(type2,'Size'),name2)



% plot_scatter({'Cone Angle Downstream','Duration'},...
%     get_data(type1,'Cone Angle Before'),get_data(type1,'Duration'),name1,...
%     get_data(type2,'Cone Angle Before'),get_data(type2,'Duration'),name2)

plot_histogram({'Cone Angle Upstream','[\theta]'},15,...
    get_data(type1,'Cone Angle After'),name1,...
    get_data(type2,'Cone Angle After'),name2)

plot_scatter({'Cone Angle Upstream','Size'},...
    get_data(type1,'Cone Angle After'),get_data(type1,'Size'),name1,...
    get_data(type2,'Cone Angle After'),get_data(type2,'Size'),name2)

% plot_scatter({'Cone Angle Upstream','Duration'},...
%     get_data(type1,'Cone Angle After'),get_data(type1,'Duration'),name1,...
%     get_data(type2,'Cone Angle After'),get_data(type2,'Duration'),name2)

plot_scatter({'Cone Angle Upstream','Cone Angle Downstream'},...
    get_data(type1,'Cone Angle Before'),get_data(type1,'Cone Angle After'),name1,...
    get_data(type2,'Cone Angle Before'),get_data(type2,'Cone Angle After'),name2)


%Alfven V
plot_histogram({'Alfven V Downstream','[km/s]'},20,...
    get_data(type1,'Alfven V Before'),name1,...
    get_data(type2,'Alfven V Before'),name2)

plot_histogram({'Alfven V Upstream','[km/s]'},20,...
    get_data(type1,'Alfven V After'),name1,...
    get_data(type2,'Alfven V After'),name2)


%Mach Number
plot_histogram({'Mach Number Downstream',''},5,...
    get_data(type1,'Mach Number Before'),name1,...
    get_data(type2,'Mach Number Before'),name2)

plot_histogram({'Mach Number Upstream',''},5,...
    get_data(type1,'Mach Number After'),name1,...
    get_data(type2,'Mach Number After'),name2)



%Event Edges Shock Angle
plot_histogram({'Leading Edge Shock Angle','[\theta]'},15,...
    get_data(type1,'Leading Shock Angle'),name1,...
    get_data(type2,'Leading Shock Angle'),name2)
plot_histogram({'Trailing Shock Angle','[\theta]'},15,...
    get_data(type1,'Trailing Shock Angle'),name1,...
    get_data(type2,'Trailing Shock Angle'),name2)


plot_scatter({'Leading Edge Shock Angle','Trailing Shock Angle'},...
    get_data(type1,'Leading Shock Angle'),get_data(type1,'Trailing Shock Angle'),name1,...
    get_data(type2,'Leading Shock Angle'),get_data(type2,'Trailing Shock Angle'),name2)



plot_histogram({'Delta B over B',''},0.1,...
    get_data(type1,'dBoverB'),name1,...
    get_data(type2,'dBoverB'),name2)

plot_histogram({'Bn over B',''},0.1,...
    get_data(type1,'BnoverB'),name1,...
    get_data(type2,'BnoverB'),name2)


plot_scatter({'Bn over B','Delta B over B'},...
    get_data(type1,'BnoverB'),get_data(type1,'dBoverB'),name1,...
    get_data(type2,'BnoverB'),get_data(type2,'dBoverB'),name2)

% %Gyration Ratio Before
% type1gyroBefore = get_data(type1,'Gyration Ratio Before');
% type2gyroBefore = get_data(type2,'Gyration Ratio Before');
% type1gyroBefore = type1gyroBefore(type1gyroBefore < 2.0);
% type2gyroBefore = type2gyroBefore(type2gyroBefore < 2.0);
% 
% plot_histogram({'Gyration Ratio Before',''},0.1,...
%     type1gyroBefore,name1,...
%     type2gyroBefore,name2)
% 
% plot_histogram({'Gyration Ratio Before',''},0.1,...
%     [type1gyroBefore;type2gyroBefore],'All')
% 
% %Gyration Ratio After
% type1gyroAfter = get_data(type1,'Gyration Ratio After');
% type2gyroAfter = get_data(type2,'Gyration Ratio After');
% type1gyroAfter = type1gyroAfter(type1gyroAfter < 2.0);
% type2gyroAfter = type2gyroAfter(type2gyroAfter < 2.0);
% 
% plot_histogram({'Gyration Ratio After',''},0.1,...
%     type1gyroAfter,name1,...
%     type2gyroAfter,name2)
% 
% plot_histogram({'Gyration Ratio After',''},0.1,...
%     [type1gyroAfter;type2gyroAfter],'All')

% 
% %Size
% type1Size = get_data(type1,'Size');
% type2Size = get_data(type2,'Size');
% type1Size = type1Size(type1Size < 5.0);
% type2Size = type2Size(type2Size < 5.0);
% 
% plot_histogram({'Size','[Re]'},0.25,...
%     type1Size,name1,...
%     type2Size,name2)
% 
% 
% %Expansion Speed
% plot_histogram({'Expansion Speed','[km/s]'},100,...
%     get_data(type1,'Expansion Speed'),name1,...
%     get_data(type2,'Expansion Speed'),name2)
% 
% 
% %Duration
% type1Duration = get_data(type1,'Duration');
% type2Duration =  get_data(type2,'Duration');
% type1Duration = type1Duration(type1Duration < 240);
% type2Duration = type2Duration(type2Duration < 240);
% 
% plot_histogram({'Duration','[s]'},30,...
%    type1Duration,name1,...
%    type2Duration,name2)
% 
% 
% % %Age
% % type1Age = get_data(type1,'Age');
% % type2Age = get_data(type2,'Age');
% % type1Age = type1Age(type1Age < 300);
% % type2Age = type2Age(type2Age < 300);
% % 
% % plot_histogram({'Age','[s]'},30,...
% %     type1Age,name1,...
% %     type2Age,name2)
% 
% 
% %Distance Traveled
% type1Distance =  get_data(type1,'Distance Traveled');
% type2Distance =  get_data(type2,'Distance Traveled');
% type1Distance = type1Distance(type1Distance < 15);
% type2Distance = type2Distance(type2Distance < 15);
% 
% plot_histogram({'Distance Traveled','[Re]'},1,...
%    type1Distance,name1,...
%    type2Distance,name2)
% 
% %Leading and Trailing Velocities (SW Frame)
% plot_scatter({'V_{leading}','V_{trailing}'},...
%     get_data(type1,'V1'),get_data(type1,'V2'),name1,...
%     get_data(type2,'V1'),get_data(type2,'V2'),name2)
% 
% %Leading and Trailing Velocities (Not SW Frame)
% plot_scatter({'Vt_{leading}','Vt_{trailing}'},...
%     get_data(type1,'Vt1'),get_data(type1,'Vt2'),name1,...
%     get_data(type2,'Vt1'),get_data(type2,'Vt2'),name2)
% 
% %Angles of Leading and Trailing Boundaries with Sun-Earth
% plot_scatter({'Theta_{leading}','Theta_{trailing}'},...
%     get_data(type1,'N1'),get_data(type1,'N2'),name1,...
%     get_data(type2,'N1'),get_data(type2,'N2'),name2)
% 
% %Angles of Leading and Trailing Boundaries with Sun-Earth, difference
% plot_scatter({'Magnetic Shear','Theta_{diff}'},...
%     get_data(type1,'Shear Angle'),get_data(type1,'N2')-get_data(type1,'N1'),name1,...
%     get_data(type2,'Shear Angle'),get_data(type2,'N2')-get_data(type2,'N1'),name2)
% 
% %Magnetic Shear vs delta_V (SW Frame)
% plot_scatter({'Magnetic Shear','V_{diff}'},...
%     get_data(type1,'Shear Angle'),get_data(type1,'BoundaryExpansion'),name1,...
%     get_data(type2,'Shear Angle'),get_data(type2,'BoundaryExpansion'),name2)
% 
% %Magnetic Shear vs delta_V (Not SW Frame)
% plot_scatter({'Magnetic Shear','Vt_{diff}'},...
%     get_data(type1,'Shear Angle'),get_data(type1,'BoundaryExpansion2'),name1,...
%     get_data(type2,'Shear Angle'),get_data(type2,'BoundaryExpansion2'),name2)
% 
% %CS angle (with sun Earth) vs delta_V
% currentSheetNormal_type1 = [get_data(type1,'CS_x'),get_data(type1,'CS_y'),get_data(type1,'CS_z')];
% for i=1:length(currentSheetNormal_type1)
%     CSSE_Angle_type1(i) = angle(currentSheetNormal_type1(i,:),[1,0,0]);
% end
% 
% currentSheetNormal_type2 = [get_data(type2,'CS_x'),get_data(type2,'CS_y'),get_data(type2,'CS_z')];
% for i=1:length(currentSheetNormal_type2)
%     CSSE_Angle_type2(i) = angle(currentSheetNormal_type2(i,:),[1,0,0]);
% end
% 
% %CS angle vs delta_V (SW Frame)
% plot_scatter({'CS Angle','V_{diff}'},...
%     CSSE_Angle_type1,get_data(type1,'BoundaryExpansion2'),name1,...
%     CSSE_Angle_type2,get_data(type2,'BoundaryExpansion2'),name2)
% % plot_scatter({'CS Angle','V_{diff}'},...
% %     CSSE_Angle_type1,get_data(type1,'BoundaryExpansion2'),name1)
% 
% 
% 
% %CS-BS angle with delta_v (sw frame)
% plot_scatter({'CS-BS Angle','V_{diff}'},...
%     get_data(type1,'CS-BS Angle'),get_data(type1,'BoundaryExpansion2'),name1,...
%     get_data(type2,'CS-BS Angle'),get_data(type2,'BoundaryExpansion2'),name2)
% 
% % plot_scatter({'CS-BS Angle','V_{diff}'},...
% %     get_data('HFA','CS-BS Angle'),get_data('HFA','BoundaryExpansion2'),'HFA')
% 
% 
% %Event Size with delta_v (sw frame)
% plot_scatter({'Event Size','V_{diff}'},...
%     log10(get_data(type1,'Size')),get_data(type1,'BoundaryExpansion2'),name1,...
%     log10(get_data(type2,'Size')),get_data(type2,'BoundaryExpansion2'),name2)
% 
% plot_scatter({'Event Size','V_{diff}'},...
%     [get_data(type1,'Size');get_data(type2,'Size')],[get_data(type1,'BoundaryExpansion2');get_data(type2,'BoundaryExpansion2')],'All')
% 
% 
% %Boundary Edges Angle with Delta_v (swFrame)
% plot_scatter({'Edges Angle','V_{diff}'},...
%     180-get_data(type1,'Event Edges Angle'),get_data(type1,'BoundaryExpansion2'),name1,...
%     180-get_data(type2,'Event Edges Angle'),get_data(type2,'BoundaryExpansion2'),name2)
% 
% 
% %Theta_Bn Angle with Delta_v (swFrame)
% plot_scatter({'Theta_{Bn} Before','V_{diff}'},...
%     get_data(type1,'Shock Angle Before'),get_data(type1,'BoundaryExpansion2'),name1,...
%     get_data(type2,'Shock Angle Before'),get_data(type2,'BoundaryExpansion2'),name2)
% 
% plot_scatter({'Theta_{Bn} After','V_{diff}'},...
%     get_data(type1,'Shock Angle After'),get_data(type1,'BoundaryExpansion2'),name1,...
%     get_data(type2,'Shock Angle After'),get_data(type2,'BoundaryExpansion2'),name2)
% 
% 
% plot_scatter({'Edges Angle','V_{diff}'},...
%     [180-get_data(type1,'Event Edges Angle');180-get_data(type2,'Event Edges Angle')],[get_data(type1,'BoundaryExpansion2');get_data(type2,'BoundaryExpansion2')],'All')
% 
% %Vdiff histograms
% plot_histogram({'V_{diff}','[km/s]'},100,...
%     get_data(type1,'BoundaryExpansion'),name1,...
%     get_data(type2,'BoundaryExpansion'),name2)
% 
% plot_histogram({'Vt_{diff}','[km/s]'},100,...
%     get_data(type1,'BoundaryExpansion2'),name1,...
%     get_data(type2,'BoundaryExpansion2'),name2)
% 
% 
% 
% %Core Densities Histogram
% plot_histogram({'Core Density STD','[\sigma]'},0.5,...
%     [get_data(type1,'Core Density STD');get_data(type2,'Core Density STD')],'All')
% plot_histogram({'Core Density STD','[\sigma]'},0.5,...
%     get_data(type1,'Core Density STD'),name1,...
%     get_data(type2,'Core Density STD'),name2)
% 
% %Core Densities CV Histogram
% plot_histogram({'Core Density CV','[\sigma]'},0.1,...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],'All')
% plot_histogram({'Core Density CV','[\sigma]'},0.1,...
%     get_data(type1,'Core Density CV'),name1,...
%     get_data(type2,'Core Density CV'),name2)
% 
% %Core Density CV with delta_v (swframe)
% plot_scatter({'Core Density CV','V_{exp}'},...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],[get_data(type1,'BoundaryExpansion');get_data(type2,'BoundaryExpansion')],'All')
% plot_scatter({'Core Density CV','V_{exp}'},...
%     get_data(type1,'Core Density CV'),get_data(type1,'BoundaryExpansion'),name1,...
%     get_data(type2,'Core Density CV'),get_data(type2,'BoundaryExpansion'),name2)
% 
% 
% 
% %Core Density CV with delta_v (swframe)
% plot_scatter({'Core Density CV','V_{diff}'},...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],[get_data(type1,'BoundaryExpansion2');get_data(type2,'BoundaryExpansion2')],'All')
% 
% 
% %Core Density CV with delta_v (swframe)
% plot_scatter({'Core Density CV','V_{diff}'},...
%     get_data(type1,'Core Density CV'),get_data(type1,'BoundaryExpansion2'),name1,...
%     get_data(type2,'Core Density CV'),get_data(type2,'BoundaryExpansion2'),name2)
% 
% 
% Core Density CV with shear
% plot_scatter({'Core Density CV','Shear Angle'},...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],[get_data(type1,'Shear Angle');get_data(type2,'Shear Angle')],'All')
% plot_scatter({'Core Density CV','Shear Angle'},...
%     get_data(type1,'Core Density CV'),get_data(type1,'Shear Angle'),name1,...
%     get_data(type2,'Core Density CV'),get_data(type2,'Shear Angle'),name2)
% 
% %Core Density CV with CS-BS Angle
% plot_scatter({'Core Density CV','CS-BS Angle'},...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],[get_data(type1,'CS-BS Angle');get_data(type2,'CS-BS Angle')],'All')
% % plot_scatter({'Core Density CV','CS-BS Angle'},...
% %     get_data(type1,'Core Density CV'),get_data(type1,'CS-BS Angle'),'HFA')
% 
% %Core Density CV with CS-BS Angle
% plot_scatter({'Core Density CV','CS-BS Angle'},...
%    get_data(type1,'Core Density CV'),get_data(type1,'CS-BS Angle'),name1,...
%    get_data(type2,'Core Density CV'),get_data(type2,'CS-BS Angle'),name2)
% % plot_scatter({'Core Density CV','CS-BS Angle'},...
% %     get_data('HFA','Core Density CV'),get_data('HFA','CS-BS Angle'),'HFA')
% 
% 
% 
% 
% 
% 
% 
% 
% %Event Edges shock angle
% plot_scatter({'Leading Shock Angle','Trailing Shock Angle'},...
%     get_data(type1,'Leading Shock Angle'),get_data(type1,'Trailing Shock Angle'),name1,...
%     get_data(type2,'Leading Shock Angle'),get_data(type2,'Trailing Shock Angle'),name2)
% 
% 
% %Event Shock Edges Angle with B
% plot_histogram({'Leading Shock Angle','[\theta]'},15,...
%     [get_data(type1,'Leading Shock Angle');get_data(type2,'Leading Shock Angle')],'All')
% 
% plot_histogram({'Leading Shock Angle','[\theta]'},15,...
%     get_data(type1,'Leading Shock Angle'),name1,get_data(type2,'Leading Shock Angle'),name2)
% 
% plot_histogram({'Trailing Shock Angle','[\theta]'},15,...
%     [get_data(type1,'Trailing Shock Angle');get_data(type2,'Trailing Shock Angle')],'All')
% 
% plot_histogram({'Trailing Shock Angle','[\theta]'},15,...
%     get_data(type1,'Trailing Shock Angle'),name1,get_data(type2,'Trailing Shock Angle'),name2)
% 
% 
% %Solar Wind Speed
% Vpre_mag_type1 = (get_data(type1,'preV_x').^2+get_data(type1,'preV_y').^2+get_data(type1,'preV_z').^2).^(1/2);
% Vpost_mag_type1 = (get_data(type1,'postV_x').^2+get_data(type1,'postV_y').^2+get_data(type1,'postV_z').^2).^(1/2);
% 
% Vpre_mag_type2 = (get_data(type2,'preV_x').^2+get_data(type2,'preV_y').^2+get_data(type2,'preV_z').^2).^(1/2);
% Vpost_mag_type2 = (get_data(type2,'postV_x').^2+get_data(type2,'postV_y').^2+get_data(type2,'postV_z').^2).^(1/2);
% % 
% plot_histogram({'Solar Wind Speed Before','[km/s]'},100,...
%      Vpre_mag_type1,name1,Vpre_mag_type2,name2)
%  
% plot_histogram({'Solar Wind Speed After','[km/s]'},100,...
%      Vpost_mag_type1,name1,Vpost_mag_type2,name2)
% plot_scatter({'Core Density CV','Solar Wind Speed Downstream'},...
%     get_data(type1,'Core Density CV'),Vpre_mag_type1,name1,...
%     get_data(type2,'Core Density CV'),Vpre_mag_type2,name2)
% 
% plot_scatter({'Core Density CV','Solar Wind Speed Upstream'},...
%     get_data(type1,'Core Density CV'),Vpost_mag_type1,name1,...
%     get_data(type2,'Core Density CV'),Vpost_mag_type2,name2)
% 
% 
% 
% %R distance vs delta_v
% XYZ_mag_type1 = (get_data(type1,'Rx').^2+get_data(type1,'Ry').^2+get_data(type1,'Rz').^2).^(1/2);
% XYZ_mag_type2 = (get_data(type2,'Rx').^2+get_data(type2,'Ry').^2+get_data(type2,'Rz').^2).^(1/2);
% 
% plot_scatter({'XY Distance','V_{diff}'},...
%     XYZ_mag_type1,get_data(type1,'BoundaryExpansion2'),name1,...
%     XYZ_mag_type2,get_data(type2,'BoundaryExpansion2'),name2)
% 
% %XY distance vs Core Density CV
% plot_scatter({'Core Density CV','XY Distance'},...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],[XYZ_mag_type1;XYZ_mag_type2],'All')
% 
% %XY distance vs Core Density CV
% plot_scatter({'Core Density CV','XY Distance'},...
%     get_data(type1,'Core Density CV'),XYZ_mag_type1,name1,...
%     get_data(type2,'Core Density CV'),XYZ_mag_type2,name2)
% 
% 
% 
% 
% 
% % vs core density STD
% plot_scatter({'Core Density STD','Core Temp Corr'},...
%     [get_data(type1,'Core Density STD');get_data(type2,'Core Density STD')],[get_data(type1,'Core Temp Corr');get_data(type2,'Core Temp Corr')],'All')
% 
% % vs core density CV
% plot_scatter({'Core Density CV','Core Temp Corr'},...
%     [get_data(type1,'Core Density CV');get_data(type2,'Core Density CV')],[get_data(type1,'Core Temp Corr');get_data(type2,'Core Temp Corr')],'All')
% 
% 
% %f vs core density STD
% plot_histogram({'Core Density STD'},0.5,...
%     get_data(type1,'Core Density STD'),name1,get_data(type2,'Core Density STD'),name2)
% 
% %f vs core density CV
% plot_histogram({'Core Density CV'},0.1,...
%     get_data(type1,'Core Density CV'),name1,get_data(type2,'Core Density CV'),name2)
% % 
% % %f vs core density CV
% % plot_histogram({'Core Density CV'},0.1,...
% %     get_data('HFA','Core Density CV'),'HFA',...
% %     get_data('SHFA','Core Density CV'),'SHFA')
% 
% %XY distance vs age
% plot_scatter({'Position','Age'},...
%     XYZ_mag_type1(get_data(type1,'Age') < 300),type1Age,name1,...
%     XYZ_mag_type2(get_data(type2,'Age') < 300),type2Age,name2)
% %
% %XY distance vs size
% plot_scatter({'Position','Size'},...
%     XYZ_mag_type1(get_data(type1,'Size') < 5),type1Size,name1,...
%     XYZ_mag_type2(get_data(type2,'Size') < 5),type2Size,name2)
% 
% %XY distance vs Distance Traveled
% plot_scatter({'Position','Distance Traveled'},...
%     XYZ_mag_type1(get_data(type1,'Distance Traveled') < 15),log(type1Distance),name1,...
%     XYZ_mag_type2(get_data(type2,'Distance Traveled') < 15),log(type2Distance),name2)

%Core Density STD with delta_V (swframe)




% % %Event Shock Angle
% plot_scatter({'Theta_{leading}','Theta_{trailing}'},...
%     get_data(type1,'N1'),get_data(type1,'N2'),name1,...
%     get_data(type2,'N1'),get_data(type2,'N2'),name2)
%
% B_pre_SS = [get_data(type1,'Bpre_x'),get_data(type1,'Bpre_y'),get_data(type1,'Bpre_z')];
% B_post_SS = [get_data(type1,'Bpost_x'),get_data(type1,'Bpost_y'),get_data(type1,'Bpost_z')];
%

% toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        min_bin = floor(min(data1)/bin_size)*bin_size;
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
    
    
    
    
%     binRange = 0:15:180;
%     A = histcounts(ECSAngle_type1,[binRange Inf]);
%     B = histcounts(ECSAngle_type2,[binRange Inf]);
%     figure
%     bar(binRange,[A;B]')
%     bar(binRange+15/2,[A;B]')
%     xlim([0,180])
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
    legend('Location','Northwest')
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
            
        case 'preN'
            Column_number = 'AX';
        case 'postN'
            Column_number = 'AY';
        case 'Shock Angle Before'
            Column_number = 'AZ';
        case 'Shock Angle After'
            Column_number = 'BA';
            
            
        case 'Cone Angle Before'
            Column_number = 'BB';
        case 'Cone Angle After'
            Column_number = 'BC';
        case 'Alfven V Before'
            Column_number = 'BD';
        case 'Mach Number Before'
            Column_number = 'BE';
        case 'Alfven V After'
            Column_number = 'BF';
        case 'Mach Number After'
            Column_number = 'BG';

        case 'E-CS Angle Before'
            Column_number = 'BH';
        case 'E-CS Angle After'
            Column_number = 'BI';
        case 'E-CS Magnitude Before'
            Column_number = 'BJ';
        case 'E-CS Magnitude After'
            Column_number = 'BK';
            
        case 'CS Speed'
            Column_number = 'BL';
        case 'Transversal Speed'
            Column_number = 'BP';
        case 'Gyration Ratio Before'
            Column_number = 'BR';
        case 'Gyration Ratio After'
            Column_number = 'BT';
        case 'Size'
            Column_number = 'CD';
        case 'Expansion Speed'
            Column_number = 'CE';
        case 'Duration'
            Column_number = 'CF';
        case 'Core Duration'
            Column_number = 'CG';
        case 'Age'
            Column_number = 'CH';
        case 'Distance Traveled'
            Column_number = 'CI';
        case 'V1'
            Column_number = 'CJ';
        case 'V2'
            Column_number = 'CK';
        case 'BoundaryExpansion'
            Column_number = 'CL';
            
            
            %         case 'BoundaryExpansion2'
            %             Column_number = 'CC';
            %         case 'N1'
            %             Column_number = 'CD';
            %         case 'N2'
            %             Column_number = 'CE';
            
            
        case 'Vt1'
            Column_number = 'CM';
        case 'Vt2'
            Column_number = 'CN';
        case 'BoundaryExpansion2'
            Column_number = 'CO';
        case 'N1'
            Column_number = 'CP';
        case 'N2'
            Column_number = 'CQ';
            
        case 'Leading Shock Angle'
            Column_number = 'CR';
        case 'Trailing Shock Angle'
            Column_number = 'CS';
        case 'Event Edges Angle'
            Column_number = 'CT';
        case 'Core Density STD'
            Column_number = 'CU';
            
        case 'Core Temp Corr'
            Column_number = 'CV';
            
        case 'Core Density CV'
            Column_number = 'CW';
            
            
        case 'dBoverB'
            Column_number = 'CX';   
            
        case 'BnoverB'
            Column_number = 'CY';
            
    end
    
    Column_range = strcat(Column_number,'2:',Column_number,'200');
    
    data = xlsread(Database_Directory,Page_number,Column_range);
    data = data(~isnan(data));
    
end
