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
name1 = 'Substructures';
type2 = 'NS';
name2 = 'No Substructures';

% type1 = 'HFA SS';
% name1 = 'Substructured HFAs';
% type2 = 'HFA NS';
% name2 = 'Non-Substructured HFAs';

% type1 = 'HFA';
% name1 = 'HFA';
% type2 = 'SHFA';
% name2 = 'SHFA';
HFAnumber = size(get_data('HFA','Rx'))
SHFAnumber = size(get_data('SHFA','Rx'))
SSnumber = size(get_data('SS','Rx'))
NSnumber = size(get_data('NS','Ry'))

%Load WIND Parameters
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/WIND_2017_2019.mat')
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/WIND_2017_2019_TD4.mat')
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/MMS_2017_2019_TD4.mat')
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/OMNI_2017_2019.mat')

%% Plots
% %% Substructures (SS) & No Substructures (NS)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Omni Averaged Data

plot_histogram({'OMNI Magnetic Field Magnitude','[nT]'},0:1:12,...
    OMNI_Bmag,'OMNI')

plot_histogram({'Magnetic Field Magnitude','[nT]'},0:1.5:6,...
    get_data(type1,'Omni B'),name1,...
    get_data(type2,'Omni B'),name2)

plot_histogram({'Normalized Magnetic Field Magnitude','[nT]'},0:1.5:6,...
    get_data(type1,'Omni B'),name1,...
    get_data(type2,'Omni B'),name2,...
    OMNI_Bmag)


plot_histogram({'OMNI Flow Speed','[km/s]'},300:100:600,...
    OMNI_Vmag,'OMNI')

plot_histogram({'Flow Speed','[km/s]'},300:100:600,...
    get_data(type1,'Omni Speed'),name1,...
    get_data(type2,'Omni Speed'),name2)

plot_histogram({'Normalized Flow Speed','[km/s]'},300:100:600,...
    get_data(type1,'Omni Speed'),name1,...
    get_data(type2,'Omni Speed'),name2,...
    OMNI_Vmag)


plot_histogram({'OMNI Ion Density','[cm^{-3}]'},0:4:20,...
    OMNI_n,'OMNI')

plot_histogram({'Ion Density','[cm^{-3}]'},0:4:20,...
    get_data(type1,'Omni Density'),name1,...
    get_data(type2,'Omni Density'),name2)

plot_histogram({'Normalized Ion Density','[cm^{-3}]'},0:5:20,...
    get_data(type1,'Omni Density'),name1,...
    get_data(type2,'Omni Density'),name2,...
    OMNI_n)


plot_histogram({'OMNI Plasma Beta','[#]'},0:1:12,...
    OMNI_Beta,'OMNI')

plot_histogram({'Plasma Beta','[#]'},0:3:12,...
    get_data(type1,'Omni Beta'),name1,...
    get_data(type2,'Omni Beta'),name2)

plot_histogram({'Normalized Plasma Beta','[#]'},0:3:09,...
    get_data(type1,'Omni Beta'),name1,...
    get_data(type2,'Omni Beta'),name2,...
    OMNI_Beta)


plot_histogram({'OMNI Alfven Mach Number','[#]'},0:2:24,...
    OMNI_Mach,'OMNI')

plot_histogram({'Alfven Mach Number','[#]'},6:6:24,...
    get_data(type1,'Omni Mach Number'),name1,...
    get_data(type2,'Omni Mach Number'),name2)

plot_histogram({'Normalized Alfven Mach Number','[#]'},6:6:18,...
    get_data(type1,'Omni Mach Number'),name1,...
    get_data(type2,'Omni Mach Number'),name2,...
    OMNI_Mach)




plot_histogram({'Core Duration','[s]'},0:30:90,...
    get_data(type1,'Core Duration'),name1,...
    get_data(type2,'Core Duration'),name2)

plot_histogram({'Total Duration','[s]'},0:30:90,...
    get_data(type1,'Duration'),name1,...
    get_data(type2,'Duration'),name2)

plot_histogram({'BoundaryExpansion','[km/s]'},-0:150:600,...
    get_data(type1,'BoundaryExpansion'),name1,...
    get_data(type2,'BoundaryExpansion'),name2)

plot_histogram({'Boundaries Expansion Size','[Re]'},0:2.5:5,...
    abs(get_data(type1,'Boundaries Expansion Size')),name1,...
    abs(get_data(type2,'Boundaries Expansion Size')),name2)

%% MMS Position
%X-Y

plot_stackedHistogram({'Size','[Re]','Closest Distance'},...
    get_data('all','Size'),1,'Size',...
    get_data('all','ClosestDistance'),2,'ClosestDistance')


plot_scatter({'Closest Distance','Duration'},...
    get_data(type1,'ClosestDistance'),get_data(type1,'Duration'),name1,...
    get_data(type2,'ClosestDistance'),get_data(type2,'Duration'),name2)


plot_scatter({'Closest Distance','Duration','Size'},...
    get_data(type1,'ClosestDistance'),get_data(type1,'Duration'),name1,...
    get_data(type2,'ClosestDistance'),get_data(type2,'Duration'),name2,...
    120*get_data(type1,'Size'),120*get_data(type2,'Size'),'Size')

plot_scatter({'X','Y'},...
    get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
    get_data(type2,'Rx'),get_data(type2,'Ry'),name2)

plot_scatter({'X','Y','Size'},...
    get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
    get_data(type2,'Rx'),get_data(type2,'Ry'),name2,...
    50*get_data(type1,'Size'),50*get_data(type2,'Size'),'size')

plot_scatter({'X','Y','Event Width'},...
    get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
    get_data(type2,'Rx'),get_data(type2,'Ry'),name2,...
    40*get_data(type1,'Event Width'),40*get_data(type2,'Event Width'),'Event Width')

plot_scatter({'X','Y','Event Height'},...
    get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
    get_data(type2,'Rx'),get_data(type2,'Ry'),name2,...
    40*get_data(type1,'Event Height'),40*get_data(type2,'Event Height'),'Event Height')

plot_scatter({'X','Y','Event Area'},...
    get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
    get_data(type2,'Rx'),get_data(type2,'Ry'),name2,...
    20*get_data(type1,'Event Area'),20*get_data(type2,'Event Area'),'Event Area')

plot_scatter({'X','Y','Cone Angle'},...
    get_data(type1,'Rx'),get_data(type1,'Ry'),name1,...
    get_data(type2,'Rx'),get_data(type2,'Ry'),name2,...
    get_data(type1,'Cone Angle'),get_data(type2,'Cone Angle'),'Cone Angle',90)

% %X-Z
% plot_scatter({'X','Z'},...
%     get_data(type1,'Rx'),get_data(type1,'Rz'),name1,...
%     get_data(type2,'Rx'),get_data(type2,'Rz'),name2)
%
% %Y-Z
% plot_scatter({'Y','Z'},...
%     get_data(type1,'Ry'),get_data(type1,'Rz'),name1,...
%     get_data(type2,'Ry'),get_data(type2,'Rz'),name2)

%% Fractional Solar Wind to Bulk Flows
%Fractional Solar Wind to Bulk Flow
% plot_histogram({'Fractional Downstream SW Speed Deflected','[km/s]'},0:0.075:0.6,...
%     get_data(type1,'Fractional downV'),name1,...
%     get_data(type2,'Fractional downV'),name2)
%
% plot_histogram({'Fractional Upstream SW Speed Deflected','[km/s]'},0:0.075:0.6,...
%     get_data(type1,'Fractional upV'),name1,...
%     get_data(type2,'Fractional upV'),name2)

% plot_scatter({'Shock Angle Before','Fractional Downstream SW Speed Deflected'},...
%     angle90(get_data(type1,'Shock Angle Before')),get_data(type1,'Fractional downV'),name1,...
%     angle90(get_data(type2,'Shock Angle Before')),get_data(type2,'Fractional downV'),name2)
%
% plot_scatter({'Shock Angle Upstream','Fractional Upstream SW Speed Deflected'},...
%     angle90(get_data(type1,'Shock Angle After')),get_data(type1,'Fractional upV'),name1,...
%     angle90(get_data(type2,'Shock Angle After')),get_data(type2,'Fractional upV'),name2)

%% Bow Shock Parameters
%Bow Shock Model Parameters


plot_histogram({'MMS Distance to BS'},0:2:12,...
    MMS_DistanceToBS,'MMS')

plot_histogram({'Distance to CS-BS Connection Point','[Re]'},0:4:8,...
    get_data(type1,'DistancetoCSBS'),name1,...
    get_data(type2,'DistancetoCSBS'),name2)

plot_histogram({'Normalized Distance to CS-BS Connection Point','[Re]'},0:4:8,...
    get_data(type1,'DistancetoCSBS'),name1,...
    get_data(type2,'DistancetoCSBS'),name2,...
    MMS_DistanceToBS)
% plot_histogram({'HFA Distance to CS-BS Connection Point','[Re]'},0:2:16,...
%     get_data('HFA SS','DistancetoCSBS'),'Substructured HFA',...
%     get_data('HFA NS','DistancetoCSBS'),'Non-Substructured HFA')


% % % % % plot_histogram({'MMS Minutely Normalized Closest Distance to Bow Shock','[Re]'},0:2:16,...
% % % % %     MMS_DistanceToBS,'MMS')

plot_histogram({'Closest Distance to Bow Shock','[Re]'},-2:2:6,...
    get_data(type1,'ClosestDistance'),name1,...
    get_data(type2,'ClosestDistance'),name2)

% % % % % plot_histogram({'MMS Minutely Normalized Closest Distance to Bow Shock','[Re]'},0:2:8,...
% % % % %     get_data(type1,'ClosestDistance'),name1,...
% % % % %     get_data(type2,'ClosestDistance'),name2,...
% % % % %     MMS_DistanceToBS)

% plot_histogram({'Closest Distance to Bow Shock','[Re]'},-2:2:8,...
%     get_data('HFA','ClosestDistance'),'HFA',...
%     get_data('SHFA','ClosestDistance'),'SHFA')
%
% plot_histogram({'MMS Minutely Normalized Closest Distance to Bow Shock','[Re]'},0:2:8,...
%     get_data('HFA','ClosestDistance'),'HFA',...
%     get_data('SHFA','ClosestDistance'),'SHFA',...
%     MMS_DistanceToBS)

%% Event Dimensions
%Event Width, Length, Area
plot_histogram({'Event Edges Angle','[\Theta]'},0:30:180,...
    get_data(type1,'Opening Angle'),name1,...
    get_data(type2,'Opening Angle'),name2)

plot_histogram({'Event Angle','[Re]'},0:30:180,...
    get_data(type1,'Event Edges Angle'),name1,...
    get_data(type2,'Event Edges Angle'),name2)



plot_histogram({'Event Width','[Re]'},0:4:16,...
    get_data(type1,'Event Width'),name1,...
    get_data(type2,'Event Width'),name2)

plot_histogram({'Event Height','[Re]'},0:3:12,...
    get_data(type1,'Event Height'),name1,...
    get_data(type2,'Event Height'),name2)

plot_histogram({'Event Area','[Re^2]'},0:3:21,...
    get_data(type1,'Event Area'),name1,...
    get_data(type2,'Event Area'),name2)

plot_scatter({'Event Width','Event Height'},...
    get_data(type1,'Event Width'),get_data(type1,'Event Height'),'Substructure',...
    get_data(type2,'Event Width'),get_data(type2,'Event Height'),'No Substructure')

%% Shock Angles at Event Edges Connection Points
%Shock Angles
plot_histogram({'Shock Angle at Leading Edge Connection Point','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'ShockAngleLeadingPoint')),name1,...
    angle90(get_data(type2,'ShockAngleLeadingPoint')),name2)

plot_histogram({'MMS Normalized Shock Angle Leading Edge','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'ShockAngleLeadingPoint')),name1,...
    angle90(get_data(type2,'ShockAngleLeadingPoint')),name2,...
    angle90(MMS_Shock_Down_Angle))

plot_histogram({'Shock Angle at Trailing Edge Connection Point','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'ShockAngleTrailingPoint')),name1,...
    angle90(get_data(type2,'ShockAngleTrailingPoint')),name2)

plot_histogram({'MMS Normalized Shock Angle Trailing Edge','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'ShockAngleTrailingPoint')),name1,...
    angle90(get_data(type2,'ShockAngleTrailingPoint')),name2,...
    angle90(MMS_Shock_Up_Angle))

% plot_scatter({'Shock Angle Downstream','Event Width'},...
%     angle90(get_data(type1,'ShockAngleLeadingPoint')),get_data(type1,'Event Width'),name1,...
%     angle90(get_data(type2,'ShockAngleLeadingPoint')),get_data(type2,'Event Width'),name2)
%
% plot_scatter({'Shock Angle Downstream','Event Height'},...
%     angle90(get_data(type1,'ShockAngleLeadingPoint')),get_data(type1,'Event Height'),name1,...
%     angle90(get_data(type2,'ShockAngleLeadingPoint')),get_data(type2,'Event Height'),name2)
%
% plot_scatter({'Shock Angle Downstream','Event Area'},...
%     angle90(get_data(type1,'ShockAngleLeadingPoint')),get_data(type1,'Event Area'),name1,...
%     angle90(get_data(type2,'ShockAngleLeadingPoint')),get_data(type2,'Event Area'),name2)
%
%
% plot_scatter({'Shock Angle Upstream','Event Width'},...
%     angle90(get_data(type1,'ShockAngleTrailingPoint')),get_data(type1,'Event Width'),name1,...
%     angle90(get_data(type2,'ShockAngleTrailingPoint')),get_data(type2,'Event Width'),name2)
%
% plot_scatter({'Shock Angle Upstream','Event Height'},...
%     angle90(get_data(type1,'ShockAngleTrailingPoint')),get_data(type1,'Event Height'),name1,...
%     angle90(get_data(type2,'ShockAngleTrailingPoint')),get_data(type2,'Event Height'),name2)
%
% plot_scatter({'Shock Angle Upstream','Event Area'},...
%     angle90(get_data(type1,'ShockAngleTrailingPoint')),get_data(type1,'Event Area'),name1,...
%     angle90(get_data(type2,'ShockAngleTrailingPoint')),get_data(type2,'Event Area'),name2)

% plot_scatter({'Shock Angle at Leading Edge Connection Point','Shock Angle at Trailing Edge Connection Point','Event Width'},...
%     angle90(get_data(type1,'ShockAngleLeadingPoint')),angle90(get_data(type1,'ShockAngleTrailingPoint')),name1,...
%     angle90(get_data(type2,'ShockAngleLeadingPoint')),angle90(get_data(type2,'ShockAngleTrailingPoint')),name2,...
%     60*get_data(type1,'Event Width'),60*get_data(type2,'Event Width'),'Event Width')
%
% plot_scatter({'Shock Angle at Leading Edge Connection Point','Shock Angle at Trailing Edge Connection Point','Event Height'},...
%     angle90(get_data(type1,'ShockAngleLeadingPoint')),angle90(get_data(type1,'ShockAngleTrailingPoint')),name1,...
%     angle90(get_data(type2,'ShockAngleLeadingPoint')),angle90(get_data(type2,'ShockAngleTrailingPoint')),name2,...
%     70*get_data(type1,'Event Height'),70*get_data(type2,'Event Height'),'Event Height')

plot_scatter({'Shock Angle at Leading Edge Connection Point','Shock Angle at Trailing Edge Connection Point','Event Area'},...
    angle90(get_data(type1,'ShockAngleLeadingPoint')),angle90(get_data(type1,'ShockAngleTrailingPoint')),name1,...
    angle90(get_data(type2,'ShockAngleLeadingPoint')),angle90(get_data(type2,'ShockAngleTrailingPoint')),name2,...
    30*get_data(type1,'Event Area'),30*get_data(type2,'Event Area'),'Event Area')

plot_scatter({'Shock Angle at Leading Edge Connection Point','Shock Angle at Trailing Edge Connection Point','Size'},...
    angle90(get_data(type1,'ShockAngleLeadingPoint')),angle90(get_data(type1,'ShockAngleTrailingPoint')),name1,...
    angle90(get_data(type2,'ShockAngleLeadingPoint')),angle90(get_data(type2,'ShockAngleTrailingPoint')),name2,...
    90*get_data(type1,'Size'),90*get_data(type2,'Size'),'Size')

%% CS-BS Angle
% %CS-BS Angle
plot_histogram({'MMS CS-BS Angle','[\Theta]'},15,...
    MMS_CSBS_Angle,'MMS')

plot_histogram({'CS-BS Angle','[\Theta]'},20,...
    get_data(type1,'CS-BS Angle'),name1,...
    get_data(type2,'CS-BS Angle'),name2)

plot_histogram({'MMS Normalized CS-BS Angle','[\Theta]'},30,...
    get_data(type1,'CS-BS Angle'),name1,...
    get_data(type2,'CS-BS Angle'),name2,...
    MMS_CSBS_Angle)

% % % % % plot_stackedHistogram({'CS-BS Angle','[\Theta]','Size'},...
% % % % %     get_data('all','CS-BS Angle'),15,'CS-BS Angle',...
% % % % %     get_data('all','Size'),3,'Size')

% plot_histogram({'HFA CS-BS Angle','[\Theta]'},0:15:180,...
%     get_data('HFA SS','CS-BS Angle'),'Substructured HFA',...
%     get_data('HFA NS','CS-BS Angle'),'Non-Substructured HFA')
%
% plot_histogram({'MMS Normalized HFA CS-BS','[\Theta]'},0:15:180,...
%     get_data('HFA SS','CS-BS Angle'),'Substructured HFA',...
%     get_data('HFA NS','CS-BS Angle'),'Non-Substructured HFA',...
%     MMS_CSBS_Angle)
%
% plot_histogram({'HFA Distance to CS-BS Point','[Re]'},0:2:10,...
%     get_data('HFA SS','DistancetoCSBS'),'Substructured HFA',...
%     get_data('HFA NS','DistancetoCSBS'),'Non-Substructured HFA')

% plot_histogram({'SHFA Closest Distance to BS','[Re]'},0:1:4,...
%     get_data('SHFA SS','ClosestDistance'),'Substructured SHFA',...
%     get_data('SHFA NS','ClosestDistance'),'Non-Substructured SHFA')
%
% plot_histogram({'HFA Closest Distance to BS','[Re]'},0:1:6,...
%     get_data('HFA SS','ClosestDistance'),'Substructured HFA',...
%     get_data('HFA NS','ClosestDistance'),'Non-Substructured HFA')
%
% plot_histogram({'Event Closest Distance to BS','[Re]'},0:1:6,...
%     get_data('HFA','ClosestDistance'),'HFA',...
%     get_data('SHFA','ClosestDistance'),'SHFA')

%% Shear Angle
%Shear Angle
plot_histogram({'Magnetic Shear Angle WIND','[\theta]'},15,...
    magneticShear,'WIND')

% plot_histogram({'Magnetic Shear Angle MMS','[\theta]'},15,...
%     MMS_magneticShear,'MMS')

plot_histogram({'Magnetic Shear Angle','[\theta]'},0:45:180,...
    get_data(type1,'Shear Angle'),name1,...
    get_data(type2,'Shear Angle'),name2)

plot_histogram({'WIND Normalized Magnetic Shear Angle','[\theta]'},0:45:180,...
    get_data(type1,'Shear Angle'),name1,...
    get_data(type2,'Shear Angle'),name2,...
    magneticShear)

plot_histogram({'MMS Normalized Magnetic Shear Angle','[\theta]'},0:45:180,...
    get_data(type1,'Shear Angle'),name1,...
    get_data(type2,'Shear Angle'),name2,...
    MMS_magneticShear)

%Magnetic Shear vs delta_V (SW Frame)
% plot_scatter({'Magnetic Shear','V_{diff}'},...
%     get_data(type1,'Shear Angle'),get_data(type1,'BoundaryExpansion'),name1,...
%     get_data(type2,'Shear Angle'),get_data(type2,'BoundaryExpansion'),name2)

%% Shock Angles
%Shock Angle Before using CS-BS Connection Point
% plot_histogram({'MMS TD Shock Angle Downstream','[\theta]'},0:15:90,...
%     angle90(MMS_Shock_Down_Angle),'MMS')

% plot_histogram({'Shock Angle Downstream','[\theta]'},0:22.5:90,...
%     angle90(get_data(type1,'Shock Angle Before')),name1,...
%     angle90(get_data(type2,'Shock Angle Before')),name2)
% 
% plot_histogram({'MMS Normalized Shock Angle Downstream','[\theta]'},0:22.5:90,...
%     angle90(get_data(type1,'Shock Angle Before')),name1,...
%     angle90(get_data(type2,'Shock Angle Before')),name2,...
%     angle90(MMS_Shock_Down_Angle))
% 
% plot_histogram({'Shock Angle Upstream','[\theta]'},0:22.5:90,...
%     angle90(get_data(type1,'Shock Angle After')),name1,...
%     angle90(get_data(type2,'Shock Angle After')),name2)
% 
% plot_histogram({'MMS Normalized Shock Angle Upstream','[\theta]'},0:22.5:90,...
%     angle90(get_data(type1,'Shock Angle After')),name1,...
%     angle90(get_data(type2,'Shock Angle After')),name2,...
%     angle90(MMS_Shock_Up_Angle))

plot_histogram({'Shock Angle Downstream','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'downstreamShockAngle')),name1,...
    angle90(get_data(type2,'downstreamShockAngle')),name2)

plot_histogram({'MMS Normalized Shock Angle Downstream','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'downstreamShockAngle')),name1,...
    angle90(get_data(type2,'downstreamShockAngle')),name2,...
    angle90(MMS_Shock_Down_Angle))

% plot_scatter({'Shock Angle Before','Duration'},...
%     angle90(get_data(type1,'downstreamShockAngle')),get_data(type1,'Duration'),name1,...
%     angle90(get_data(type2,'downstreamShockAngle')),get_data(type2,'Duration'),name2)
%
% plot_scatter({'Shock Angle Downstream','Size'},...
%     angle90(get_data(type1,'downstreamShockAngle')),get_data(type1,'Size'),name1,...
%     angle90(get_data(type2,'downstreamShockAngle')),get_data(type2,'Size'),name2)

%Shock Angle After using CS-BS Connection Point
% plot_histogram({'MMS TD Shock Angle Upstream','[\theta]'},0:15:90,...
%     angle90(MMS_Shock_Up_Angle),'MMS')
%


plot_histogram({'Shock Angle Upstream','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'upstreamShockAngle')),name1,...
    angle90(get_data(type2,'upstreamShockAngle')),name2)

plot_histogram({'MMS Normalized Shock Angle Upstream','[\theta]'},0:22.5:90,...
    angle90(get_data(type1,'upstreamShockAngle')),name1,...
    angle90(get_data(type2,'upstreamShockAngle')),name2,...
    angle90(MMS_Shock_Up_Angle))
%
% plot_scatter({'Shock Angle Upstream','Duration'},...
%     angle90(get_data(type1,'upstreamShockAngle')),get_data(type1,'Duration'),name1,...
%     angle90(get_data(type2,'upstreamShockAngle')),get_data(type2,'Duration'),name2)
%
% plot_scatter({'Shock Angle Upstream','Size'},...
%     angle90(get_data(type1,'upstreamShockAngle')),get_data(type1,'Size'),name1,...
%     angle90(get_data(type2,'upstreamShockAngle')),get_data(type2,'Size'),name2)
%
% plot_scatter({'Shock Angle Downstream','Shock Angle Upstream'},...
%     angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
%     angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2)

plot_scatter({'Shock Angle Downstream','Shock Angle Upstream','Event Area'},...
    angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
    angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2,...
    30*get_data(type1,'Event Area'),30*get_data(type2,'Event Area'),'Event Area')


plot_scatter({'Shock Angle Downstream','Shock Angle Upstream','Size'},...
    angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
    angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2,...
    150*get_data(type1,'Size'),150*get_data(type2,'Size'),'Size')


plot_scatter({'Shock Angle Downstream','Shock Angle Upstream','Duration'},...
    angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
    angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2,...
    5*get_data(type1,'Duration'),5*get_data(type2,'Duration'),'Duration')



plot_scatter({'Shock Angle Downstream','Shock Angle Upstream','Opening Angle'},...
    angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
    angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2,...
    5*get_data(type1,'Opening Angle'),5*get_data(type2,'Opening Angle'),'Opening Angle')

%% E-CS Orientation

%E-CS Before Mag
plot_histogram({'E-CS Magnitude Downstream WIND','[nT km/s]'},-400:50:400,...
    E_Down_Mag,'WIND')

% plot_histogram({'E-CS Magnitude Downstream MMS','[nT km/s]'},-2000:400:2000,...
%     MMS_E_Down_Mag,'MMS')

plot_histogram({'E-CS Magnitude Downstream','[nT km/s]'},-500:250:500,...
    get_data('HFA SS','E-CS Magnitude Before'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Magnitude Before'),'Non-Substructured HFAs')

plot_histogram({'WIND Normalized E-CS Magnitude Downstream','[nT km/s]'},-500:250:500,...
    get_data('HFA SS','E-CS Magnitude Before'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Magnitude Before'),'Non-Substructured HFAs',...
    E_Down_Mag)

% plot_histogram({'MMS Normalized E-CS Magnitude Downstream','[nT km/s]'},-2000:500:2000,...
%     get_data('HFA SS','E-CS Magnitude Before'),'Substructured HFAs',...
%     get_data('HFA NS','E-CS Magnitude Before'),'Non-Substructured HFAs',...
%     MMS_E_Down_Mag)

%E-CS Before Angle
plot_histogram({'E-CS Angle Downstream WIND','[\theta]'},15,...
    E_Down_Angle,'WIND')

% plot_histogram({'E-CS Angle Downstream MMS','[\theta]'},15,...
%     MMS_E_Down_Angle,'MMS')

plot_histogram({'E-CS Angle Downstream','[\theta]'},90,...
    get_data('HFA SS','E-CS Angle Before'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),'Non-Substructured HFAs')

plot_histogram({'WIND Normalized E-CS Angle Downstream','[\theta]'},90,...
    get_data('HFA SS','E-CS Angle Before'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),'Non-Substructured HFAs',...
    E_Down_Angle)

% plot_histogram({'MMS Normalized E-CS Angle Downstream','[\theta]'},90,...
%     get_data('HFA SS','E-CS Angle Before'),'Substructured HFAs',...
%     get_data('HFA NS','E-CS Angle Before'),'Non-Substructured HFAs',...
%     MMS_E_Down_Angle)

plot_scatter({'E-CS Angle Downstream','E-CS Magnitude Downstream','Event Area'},...
    get_data('HFA SS','E-CS Angle Before'),get_data('HFA SS','E-CS Magnitude Before'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),get_data('HFA NS','E-CS Magnitude Before'),'Non-Substructured HFAs',...
    10*get_data('HFA SS','Event Area'),10*get_data('HFA NS','Event Area'),'Event Area')

%E-CS After Mag
plot_histogram({'E-CS Magnitude Upstream WIND','[nT km/s]'},-400:50:400,...
    -E_Up_Mag,'WIND')

% plot_histogram({'E-CS Magnitude Upstream MMS','[nT km/s]'},-2000:400:2000,...
%     -MMS_E_Up_Mag,'MMS')

plot_histogram({'E-CS Magnitude Upstream','[nT km/s]'},-500:250:500,...
    -get_data('HFA SS','E-CS Magnitude After'),'Substructured HFAs',...
    -get_data('HFA NS','E-CS Magnitude After'),'Non-Substructured HFAs')

plot_histogram({'WIND Normalized E-CS Magnitude Upstream','[nT km/s]'},-500:250:500,...
    -get_data('HFA SS','E-CS Magnitude After'),'Substructured HFAs',...
    -get_data('HFA NS','E-CS Magnitude After'),'Non-Substructured HFAs',...
    -E_Up_Mag)
% 
% plot_histogram({'MMS Normalized E-CS Magnitude Upstream','[nT km/s]'},-2000:1000:2000,...
%     -get_data('HFA SS','E-CS Magnitude After'),'Substructured HFAs',...
%     -get_data('HFA NS','E-CS Magnitude After'),'Non-Substructured HFAs',...
%     -MMS_E_Up_Mag)

%E-CS After Angle
plot_histogram({'E-CS Angle Upstream WIND','[\theta]'},15,...
    180-E_Up_Angle,'WIND')
%
% plot_histogram({'E-CS Angle Upstream MMS','[\theta]'},15,...
%     180-MMS_E_Up_Angle,'MMS')

plot_histogram({'E-CS Angle Upstream','[\theta]'},90,...
    180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs')

plot_histogram({'WIND Normalized E-CS Angle Upstream','[\theta]'},90,...
    180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs',...
    180-E_Up_Angle)

% plot_histogram({'MMS Normalized E-CS Angle Upstream','[\theta]'},90,...
%     180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
%     180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs',...
%     180-MMS_E_Up_Angle)

plot_scatter({'E-CS Angle Downstream','E-CS Angle Upstream'},...
    get_data('HFA SS','E-CS Angle Before'),180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs')

plot_scatter({'E-CS Angle Upstream','E-CS Magnitude Upstream','Event Area'},...
    180-get_data('HFA SS','E-CS Angle After'),-get_data('HFA SS','E-CS Magnitude After'),'Substructured HFAs',...
    180-get_data('HFA NS','E-CS Angle After'),-get_data('HFA NS','E-CS Magnitude After'),'Non-Substructured HFAs',...
    20*get_data('HFA SS','Event Area'),20*get_data('HFA NS','Event Width'),'Event Area')


plot_scatter({'E-CS Angle Downstream','E-CS Angle Upstream','Size'},...
    get_data('HFA SS','E-CS Angle Before'),180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs',...
    150*get_data('HFA SS','Size'),150*get_data('HFA NS','Size'),'Size')

plot_scatter({'E-CS Angle Downstream','E-CS Angle Upstream','Duration'},...
    get_data('HFA SS','E-CS Angle Before'),180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs',...
    5*get_data('HFA SS','Duration'),5*get_data('HFA NS','Duration'),'Duration')


plot_scatter({'E-CS Angle Downstream','E-CS Angle Upstream','Event Area'},...
    get_data('HFA SS','E-CS Angle Before'),180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs',...
    45*get_data('HFA SS','Event Area'),45*get_data('HFA NS','Event Area'),'Event Area')


plot_scatter({'E-CS Angle Downstream','E-CS Angle Upstream','Opening Angle'},...
    get_data('HFA SS','E-CS Angle Before'),180-get_data('HFA SS','E-CS Angle After'),'Substructured HFAs',...
    get_data('HFA NS','E-CS Angle Before'),180-get_data('HFA NS','E-CS Angle After'),'Non-Substructured HFAs',...
    4*get_data('HFA SS','Opening Angle'),4*get_data('HFA NS','Opening Angle'),'Opening Angle')

%% Cone Angle
%Choosing Cone Angle that is closest to 0 or 180.
plot_histogram({'Cone Angle','[\theta]'},30,...
    coneAngle,'WIND')

plot_histogram({'Cone Angle','[\theta]'},30,...
    calculate_radialConeAngle(type1)',name1,...
    calculate_radialConeAngle(type2)',name2)

plot_histogram({'Normalized Cone Angle','[\theta]'},45,...
    calculate_radialConeAngle(type1)',name1,...
    calculate_radialConeAngle(type2)',name2,...
    coneAngle)
%%%%
% % % % plot_scatter({'Cone Angle','Duration'},...
% % % %     calculate_radialConeAngle(type1)',get_data(type1,'Duration'),name1,...
% % % %     calculate_radialConeAngle(type2)',get_data(type2,'Duration'),name2)
% % % %
% % % % plot_scatter({'Cone Angle','Event Width'},...
% % % %     calculate_radialConeAngle(type1)',get_data(type1,'Event Width'),name1,...
% % % %     calculate_radialConeAngle(type2)',get_data(type2,'Event Width'),name2)
% % % %
% % % % plot_scatter({'Cone Angle','Event Height'},...
% % % %     calculate_radialConeAngle(type1)',get_data(type1,'Event Height'),name1,...
% % % %     calculate_radialConeAngle(type2)',get_data(type2,'Event Height'),name2)
% % % %
% % % % plot_scatter({'Cone Angle','Event Area'},...
% % % %     calculate_radialConeAngle(type1)',get_data(type1,'Event Area'),name1,...
% % % %     calculate_radialConeAngle(type2)',get_data(type2,'Event Area'),name2)



% % % % plot_stackedHistogram({'Cone Angle','[\Theta]','Event Area'},...
% % % %     get_data('all','Cone Angle'),30,'Cone Angle',...
% % % %     get_data('all','Event Area'),0:5:30,'Event Area')
%%%%
% plot_histogram({'HFA Cone Angle','[\theta]'},30,...
%     calculate_radialConeAngle('HFA SS')','Substructured HFA',...
%     calculate_radialConeAngle('HFA NS')','Non-Substructured HFA')
%
% plot_histogram({'Normalized HFA Cone Angle','[\theta]'},30,...
%     calculate_radialConeAngle('HFA SS')','Substructured HFA',...
%     calculate_radialConeAngle('HFA NS')','Non-Substructured HFA',...
%     coneAngle)
%
% plot_scatter({'HFA Cone Angle','Size'},...
%     calculate_radialConeAngle('HFA SS'),get_data('HFA SS','Size'),'HFA SS',...
%     calculate_radialConeAngle('HFA NS'),get_data('HFA NS','Size'),'HFA NS')
% %%%%
% plot_histogram({'SHFA Cone Angle','[\theta]'},30,...
%     calculate_radialConeAngle('SHFA SS')','Substructured SHFA',...
%     calculate_radialConeAngle('SHFA NS')','Non-Substructured SHFA')
%
% plot_histogram({'Normalized SHFA Cone Angle','[\theta]'},30,...
%     calculate_radialConeAngle('SHFA SS')','Substructured SHFA',...
%     calculate_radialConeAngle('SHFA NS')','Non-Substructured SHFA',...
%     coneAngle)

%% Alfven V and Mach Numbers

% %Alfven V
% plot_histogram({'Alfven V Downstream','[km/s]'},0:20:80,...
%     get_data(type1,'Alfven V Before'),name1,...
%     get_data(type2,'Alfven V Before'),name2)
%
% plot_histogram({'Alfven V Upstream','[km/s]'},0:20:80,...
%     get_data(type1,'Alfven V After'),name1,...
%     get_data(type2,'Alfven V After'),name2)
%
% plot_histogram({'Alfven V'},0:10:110,...
%     V_Alfven,'WIND')
%
% plot_histogram({'Normalized Alfven V Downstream','[km/s]'},0:20:80,...
%     get_data(type1,'Alfven V Before'),name1,...
%     get_data(type2,'Alfven V Before'),name2, V_Alfven)
%
% plot_histogram({'Normalized Alfven V Upstream','[km/s]'},0:20:80,...
%     get_data(type1,'Alfven V After'),name1,...
%     get_data(type2,'Alfven V After'),name2, V_Alfven)


%Mach Number
plot_histogram({'Mach Number Downstream',''},6:3:18,...
    get_data(type1,'Mach Number Before'),name1,...
    get_data(type2,'Mach Number Before'),name2)

plot_histogram({'Mach Number Upstream',''},6:3:18,...
    get_data(type1,'Mach Number After'),name1,...
    get_data(type2,'Mach Number After'),name2)

% plot_scatter({'Mach Number Downstream','Mach Number Upstream'},...
%     get_data(type1,'Mach Number Before'),get_data(type1,'Mach Number After'),name1,...
%     get_data(type2,'Mach Number Before'),get_data(type2,'Mach Number After'),name2)
%
% plot_histogram({'Mach Number Upstream - Mach Number Downstream',''},-20:3:20,...
%     get_data(type1,'Mach Number After')- get_data(type1,'Mach Number Before'),name1,...
%     get_data(type2,'Mach Number After')- get_data(type2,'Mach Number Before'),name2)

plot_histogram({'Mach Number',''},0:2:20,...
    MachNumber,'WIND')

plot_histogram({'Normalized Mach Number Downstream',''},9:4:17,...
    get_data(type1,'Mach Number Before'),name1,...
    get_data(type2,'Mach Number Before'),name2, MachNumber)

plot_histogram({'Normalized Mach Number Upstream',''},9:4:17,...
    get_data(type1,'Mach Number After'),name1,...
    get_data(type2,'Mach Number After'),name2, MachNumber)

%% Solar Wind Speed
plot_histogram({'Solar Wind Speed','[km/s]'},[200:50:700],...
    vecnorm(vdata,2,2),'WIND')

plot_histogram({'Solar Wind Speed Downstream','[km/s]'},[400:50:600],...
    get_data(type1,'downstreamSpeed'),name1,get_data(type2,'downstreamSpeed'),name2)
%
plot_histogram({'Solar Wind Speed Upstream','[km/s]'},[400:50:600],...
    get_data(type1,'upstreamSpeed'),name1,get_data(type2,'upstreamSpeed'),name2)

plot_histogram({'Normalized Solar Wind Speed Downstream','[km/s]'},[400:100:600],...
    get_data(type1,'downstreamSpeed'),name1,get_data(type2,'downstreamSpeed'),name2,vecnorm(vdata,2,2))
%
plot_histogram({'Normalized Solar Wind Speed Upstream','[km/s]'},[400:100:600],...
    get_data(type1,'upstreamSpeed'),name1,get_data(type2,'upstreamSpeed'),name2,vecnorm(vdata,2,2))

%% Magnetic Field Magnitude
plot_histogram({'Solar Wind Magnetic Field Magnitude','[nT]'},0:3:12,...
    Bdatamag,'WIND')

% type_1_B_Before = [get_data(type1,'Bpre_x'),get_data(type1,'Bpre_y'),get_data(type1,'Bpre_z')];
% type_1_B_After = [get_data(type1,'Bpost_x'),get_data(type1,'Bpost_y'),get_data(type1,'Bpost_z')];
%
% type_2_B_Before = [get_data(type2,'Bpre_x'),get_data(type2,'Bpre_y'),get_data(type2,'Bpre_z')];
% type_2_B_After = [get_data(type2,'Bpost_x'),get_data(type2,'Bpost_y'),get_data(type2,'Bpost_z')];
%
% type_1_B_Before = vecnorm(type_1_B_Before,2,2);
% type_1_B_After = vecnorm(type_1_B_After,2,2);
% type_2_B_Before = vecnorm(type_2_B_Before,2,2);
% type_2_B_After = vecnorm(type_2_B_After,2,2);
%
%
% plot_histogram({'Magnetic Field Downstream','[km/s]'},0:2:12,...
%     type_1_B_Before,name1,type_1_B_After,name2)
%
% plot_histogram({'Magnetic Field Upstream','[km/s]'},0:2:12,...
%     type_1_B_After,name1,type_2_B_After,name2)
%
% plot_histogram({'Normalized Magnetic Field Downstream','[km/s]'},0:2:12,...
%    type_1_B_Before,name1,type_1_B_After,name2,Bdatamag)
%
% plot_histogram({'Normalized Solar Magnetic Field Upstream','[km/s]'},0:2:12,...
%    type_1_B_After,name1,type_2_B_After,name2,Bdatamag)


plot_histogram({'Magnetic Field Downstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpre'),name1,get_data(type2,'Bpre'),name2)

plot_histogram({'Magnetic Field Upstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpost'),name1,get_data(type2,'Bpost'),name2)

plot_histogram({'Normalized Magnetic Field Downstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpre'),name1,get_data(type2,'Bpre'),name2,Bdatamag)

plot_histogram({'Normalized Magnetic Field Upstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpost'),name1,get_data(type2,'Bpost'),name2,Bdatamag)


plot_histogram({'Magnetic Field Downstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpre_cs'),name1,get_data(type2,'Bpre_cs'),name2)

plot_histogram({'Magnetic Field Upstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpost_cs'),name1,get_data(type2,'Bpost_cs'),name2)

plot_histogram({'Normalized Magnetic Field Downstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpre_cs'),name1,get_data(type2,'Bpre_cs'),name2,Bdatamag)

plot_histogram({'Normalized Magnetic Field Upstream','[nT]'},3:1.5:6,...
    get_data(type1,'Bpost_cs'),name1,get_data(type2,'Bpost_cs'),name2,Bdatamag)

%% DD parameters
% plot_histogram({'Delta B over B WIND'},0.1,...
%     deltaBOverBmax,'WIND')
%
% plot_histogram({'Delta B over B MMS'},0.1,...
%     MMS_deltaBOverBmax,'MMS')
%
%
% plot_histogram({'Bn over B WIND'},0.1,...
%     BnmaxOverBmax,'WIND')
%
% plot_histogram({'Bn over B MMS'},0.1,...
%     MMS_BnmaxOverBmax,'MMS')
%
%
% plot_scatter({'Bn over B','Delta B over B'},...
%     BnmaxOverBmax,deltaBOverBmax,'WIND')

%% Size, Expansion Speed, Transversal Speed, Duration Distance Traveled
%Size
type1Size = rmoutliers(get_data(type1,'Size'));
type2Size = rmoutliers(get_data(type2,'Size'));
plot_histogram({'Size','[Re]'},0:1:5,...
    type1Size,name1,...
    type2Size,name2)

%Expansion Speed
plot_histogram({'Expansion Speed','[km/s]'},-200:100:200,...
    get_data(type1,'Expansion Speed'),name1,...
    get_data(type2,'Expansion Speed'),name2)

%Boundary Expansion in SW frame, v2 projected into v1
plot_histogram({'BoundaryExpansion','[km/s]'},-150:150:300,...
    get_data(type1,'BoundaryExpansion'),name1,...
    get_data(type2,'BoundaryExpansion'),name2)

plot_histogram({'BoundaryExpansionC','[km/s]'},-0:150:600,...
    get_data(type1,'BoundaryExpansionC'),name1,...
    get_data(type2,'BoundaryExpansionC'),name2)



plot_scatter({'BoundaryExpansionC','Size'},...
    get_data(type1,'BoundaryExpansionC'),get_data(type1,'Size'),name1,...
    get_data(type2,'BoundaryExpansionC'),get_data(type2,'Size'),name2)



plot_scatter({'Expansion Speed (SW Frame)','Size'},...
    get_data(type1,'HFA Expansion Speed SW'),get_data(type1,'Size'),name1,...
    get_data(type2,'HFA Expansion Speed SW'),get_data(type2,'Size'),name2)

%Transversal Speed
plot_histogram({'MMS Transversal Speed','[km/s]'},0:150:2000,...
    MMS_Transversal_Speed,'MMS')

plot_histogram({'Transversal Speed','[km/s]'},0:150:600,...
    get_data(type1,'Transversal Speed'),name1,...
    get_data(type2,'Transversal Speed'),name2)

plot_histogram({'MMS Normalized Transversal Speed','[km/s]'},0:150:600,...
    get_data(type1,'Transversal Speed'),name1,...
    get_data(type2,'Transversal Speed'),name2,...
    MMS_Transversal_Speed)

%Duration
type1Duration = rmoutliers(get_data(type1,'Duration'));
type2Duration = rmoutliers(get_data(type2,'Duration'));

plot_histogram({'Duration','[s]'},0:30:90,...
    type1Duration,name1,...
    type2Duration,name2)

% %Distance Traveled
% type1Distance =  rmoutliers(get_data(type1,'Distance Traveled'));
% type2Distance =  rmoutliers(get_data(type2,'Distance Traveled'));
% type1Distance = type1Distance(type1Distance < 15);
% type2Distance = type2Distance(type2Distance < 15);
%
% plot_histogram({'Distance Traveled','[Re]'},0:1:12,...
%     type1Distance,name1,...
%     type2Distance,name2)

%% Leading And Trailing Orientations &
%Leading and Trailing Velocities (SW Frame)
plot_scatter({'V_{leading}','V_{trailing}'},...
    get_data(type1,'V1'),get_data(type1,'V2'),name1,...
    get_data(type2,'V1'),get_data(type2,'V2'),name2)

%CS angle (with sun Earth) vs delta_V
currentSheetNormal_HFA = [get_data('HFA','CS_x'),get_data('HFA','CS_y'),get_data('HFA','CS_z')];
currentSheetNormal_type1 = currentSheetNormal_HFA(substructure == 1,:);
currentSheetNormal_type2 = currentSheetNormal_HFA(substructure == 0,:);

leadingEdgeNormal_HFA = [get_data('HFA','Leading Timing Normal_x'),get_data('HFA','Leading Timing Normal_y'),get_data('HFA','Leading Timing Normal_z')];
leadingEdgeNormal_type1 = leadingEdgeNormal_HFA(substructure == 1,:);
leadingEdgeNormal_type2 = leadingEdgeNormal_HFA(substructure == 0,:);


trailingEdgeNormal_HFA = [get_data('HFA','Trailing Timing Normal_x'),get_data('HFA', 'Trailing Timing Normal_y'),get_data('HFA','Trailing Timing Normal_z')];
trailingEdgeNormal_type1 = trailingEdgeNormal_HFA(substructure == 1,:);
trailingEdgeNormal_type2 = trailingEdgeNormal_HFA(substructure == 0,:);

% for i=1:length(currentSheetNormal_type1)
%     CSSE_Angle_type1(i) = angle(currentSheetNormal_type1(i,:),[1,0,0]);
%     CSLead_Angle_type1(i) = min(angle(currentSheetNormal_type1(i,:),leadingEdgeNormal_type1(i,:)),angle(-currentSheetNormal_type1(i,:),leadingEdgeNormal_type1(i,:)));
%     CSTrail_Angle_type1(i) = min(angle(currentSheetNormal_type1(i,:),trailingEdgeNormal_type1(i,:)),angle(-currentSheetNormal_type1(i,:),trailingEdgeNormal_type1(i,:)));
% end
% %
%
% for i=1:length(currentSheetNormal_type2)
%     CSSE_Angle_type2(i) = angle(currentSheetNormal_type2(i,:),[1,0,0]);
%     CSLead_Angle_type2(i) = min(angle(currentSheetNormal_type2(i,:),leadingEdgeNormal_type2(i,:)),angle(-currentSheetNormal_type2(i,:),leadingEdgeNormal_type2(i,:)));
%     CSTrail_Angle_type2(i) = min(angle(currentSheetNormal_type2(i,:),trailingEdgeNormal_type2(i,:)),angle(-currentSheetNormal_type2(i,:),trailingEdgeNormal_type2(i,:)));
% end

%% Misc Scatters
%Angles
plot_histogram({'Event Edges Angle','[\theta]'},45,...
    get_data(type1,'Event Edges Angle'),name1,...
    get_data(type2,'Event Edges Angle'),name2)

plot_scatter({'Event Edges Angle','Size'},...
    get_data(type1,'Event Edges Angle'),(get_data(type1,'Size')),name1,...
    get_data(type2,'Event Edges Angle'),(get_data(type2,'Size')),name2)

plot_scatter({'Event Edges Angle','Rx'},...
    get_data(type1,'Event Edges Angle'), get_data(type1,'Rx'),name1,...
    get_data(type2,'Event Edges Angle'), get_data(type2,'Rx'),name2)

R1 = (get_data(type1,'Rx').^2 + get_data(type1,'Ry').^2 + get_data(type1,'Rz').^2).^(1/2);
R2 = (get_data(type2,'Rx').^2 + get_data(type2,'Ry').^2 + get_data(type2,'Rz').^2).^(1/2);

plot_scatter({'Event Edges Angle','R'},...
    get_data(type1,'Event Edges Angle'), R1,name1,...
    get_data(type2,'Event Edges Angle'), R2,name2)

plot_scatter({'Event Edges Angle','Cone Angle'},...
    get_data(type1,'Event Edges Angle'), type1_coneAngle',name1,...
    get_data(type2,'Event Edges Angle'), type2_coneAngle',name2)


plot_scatter({'Event Edges Angle','Duration'},...
    get_data(type1,'Event Edges Angle'), get_data(type1,'Duration'),name1,...
    get_data(type2,'Event Edges Angle'), get_data(type2,'Duration'),name2)

plot_scatter({'Event Edges Angle','Transversal Speed'},...
    get_data(type1,'Event Edges Angle'), get_data(type1,'Transversal Speed'),name1,...
    get_data(type2,'Event Edges Angle'), get_data(type2,'Transversal Speed'),name2)



%CS Leading/Trailing Edge Angle
plot_histogram({'CS-Leading Edge Angle','[\theta]'},0:15:90,...
    CSLead_Angle_type1',name1,...
    CSLead_Angle_type2',name2)

plot_histogram({'CS-Trailing Edge Angle','[\theta]'},0:15:90,...
    CSTrail_Angle_type1',name1,...
    CSTrail_Angle_type2',name2)

plot_scatter({'CS-Leading Edge Angle','CS-Trailing Edge Angle'},...
    CSLead_Angle_type1',CSTrail_Angle_type1',name1,...
    CSLead_Angle_type2',CSTrail_Angle_type2',name2)

% %CS angle (with sun Earth) vs delta_V
currentSheetNormal_type1 = [get_data(type1,'CS_x'),get_data(type1,'CS_y'),get_data(type1,'CS_z')];
currentSheetNormal_type2 = [get_data(type2,'CS_x'),get_data(type2,'CS_y'),get_data(type2,'CS_z')];

for i=1:length(currentSheetNormal_type1)
    CSSE_Angle_type1(i) = angle(currentSheetNormal_type1(i,:),[1,0,0]);
end
%

for i=1:length(currentSheetNormal_type2)
    CSSE_Angle_type2(i) = angle(currentSheetNormal_type2(i,:),[1,0,0]);
end

% %CS angle vs delta_V (SW Frame)
plot_scatter({'CS Angle','V_{diff}'},...
    CSSE_Angle_type1,get_data(type1,'BoundaryExpansion2'),name1,...
    CSSE_Angle_type2,get_data(type2,'BoundaryExpansion2'),name2)

%CS-BS vs delta_V (SW Frame)
plot_scatter({'CS-BS Angle','V_{diff}'},...
    get_data('HFA','CS-BS Angle'),get_data('HFA','BoundaryExpansion2'),'HFA')

% %Event Size with delta_v (sw frame)
plot_scatter({'Event Size','V_{diff}'},...
    log10(get_data(type1,'Size')),get_data(type1,'BoundaryExpansion2'),name1,...
    log10(get_data(type2,'Size')),get_data(type2,'BoundaryExpansion2'),name2)

% %Boundary Edges Angle with Delta_v (swFrame)
plot_scatter({'Edges Angle','V_{diff}'},...
    get_data(type1,'Event Edges Angle'),get_data(type1,'BoundaryExpansion2'),name1,...
    get_data(type2,'Event Edges Angle'),get_data(type2,'BoundaryExpansion2'),name2)

% %Theta_Bn Angle with Delta_v (swFrame)
plot_scatter({'Theta_{Bn} Before','V_{diff}'},...
    get_data(type1,'Shock Angle Before'),get_data(type1,'BoundaryExpansion2'),name1,...
    get_data(type2,'Shock Angle Before'),get_data(type2,'BoundaryExpansion2'),name2)
%
plot_scatter({'Theta_{Bn} After','V_{diff}'},...
    get_data(type1,'Shock Angle After'),get_data(type1,'BoundaryExpansion2'),name1,...
    get_data(type2,'Shock Angle After'),get_data(type2,'BoundaryExpansion2'),name2)
%
%
plot_scatter({'Edges Angle','V_{diff}'},...
    [get_data(type1,'Event Edges Angle');get_data(type2,'Event Edges Angle')],[get_data(type1,'BoundaryExpansion2');get_data(type2,'BoundaryExpansion2')],'All')
%
% %Vdiff histograms
plot_histogram({'V_{diff}','[km/s]'},-300:50:300,...
    get_data(type1,'BoundaryExpansion'),name1,...
    get_data(type2,'BoundaryExpansion'),name2)

%difference in denssity
type1_diffN = abs(get_data(type1,'postN') - get_data(type1,'preN'))./mean([get_data(type1,'postN');get_data(type1,'preN')]);
type2_diffN = abs(get_data(type2,'postN') - get_data(type2,'preN'))./mean([get_data(type2,'postN');get_data(type2,'preN')]);

plot_histogram({'Difference in Density Across Event','[cm^-3]'},0:0.2:2,...
    type1_diffN,name1,type2_diffN,name2)





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions
%Plot the histogram, up to 3 data points
function [] = plot_histogram(parameter,bin_size,data1,label1,data2,label2,data3)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 350 350])
    %     co = [0 0 1;
    %         0 1 0;
    %         1 0 0;
    %         0 0 0];
    %     set(gcf,'defaultAxesColorOrder',co)
    set(gcf,'color','w');
    
    if nargin == 4
        data2=[];
        label2=[];
    end
    
    if size(bin_size) == 1
        
        %     data1=rmoutliers(data1,'ThresholdFactor',25);
        %     data2=rmoutliers(data2,'ThresholdFactor',25);
        if min([data1;data2])  >= 0 && max([data1;data2]) <= 1
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        elseif min([data1;data2])  >= -1 && max([data1;data2]) <= 1
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        elseif min([data1;data2])  >= 0 && max([data1;data2]) <= 75
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        elseif min([data1;data2])  >= 0 && max([data1;data2]) <= 180
            binRange = 0:bin_size:180;
        elseif  min([data1;data2]) >= -90 && max([data1;data2]) <= 90
            binRange = -90:bin_size:90;
        else
            binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
        end
        custom_bins=0;
    else
        binRange=bin_size;
        custom_bins = 1;
        bin_size = binRange(3)-binRange(2);
    end
    
    %Normalization
    if nargin == 7
        %         C = histcounts(data3,[-Inf binRange max([data1;data2])]);
        C = histcounts(data3,[-Inf binRange Inf]);
        
        if sum(contains(parameter,'ClosestDistance')) > 0
            %2 minutes per data point, per hour is 30 times.
            C = C./30;
        end
        %         C = C./sum(C);
        Cerror = C.^(-1/2);
        Cerror(Cerror == inf) = NaN;
        Cerror(Cerror == 0) = NaN;
        if sum(contains(parameter,'MMS')) > 0
            yAxisLabel = '# of Events normalized to MMS observations';
        else
            yAxisLabel = '# of Events normalized to WIND observations';
        end
    else
        C=1;
        Cerror = 1;
        yAxisLabel = '# of Events / Total # Events';
    end
    
    %weight = C./sum(C);
    if nargin ~= 4 %Wewighting from Spacecraft Data Over All Time Range
        A = histcounts(data1,[-Inf binRange Inf]);
        B = histcounts(data2,[-Inf binRange Inf]);
        
        A_normalized = A./C;
        B_normalized =B./C;
        
        if length(C) == 1
            A_normalized = A./sum(A);
            B_normalized =B./sum(B);
            Aerror = (A).^(-1/2)./sum(A);
            Berror = (B).^(-1/2)./sum(B);
            
            %             Aerror = A_normalized.*(A.^(-1) + sum(A).^(-1) ).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) + sum(B).^(-1) ).^(1/2);
        else
            Aerror = A_normalized.*(A.^(-1) +      C.^(-1)).^(1/2);
            Berror = B_normalized.*(B.^(-1) +      C.^(-1)).^(1/2);
        end
        
        if A(1) == 0 && B(1) == 0
            A = A(2:end);
            A_normalized = A_normalized(2:end);
            
            B = B(2:end);
            B_normalized = B_normalized(2:end);
            
            if length(C) ~= 1
                C = C(2:end);
                Cerror = Cerror(2:end);
                Aerror = Aerror(2:end);
                Berror = Berror(2:end);
            else
                Aerror = Aerror(2:end);
                Berror = Berror(2:end);
            end
        else
            binRange = [binRange(1) - bin_size, binRange];
        end
        xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        
        
        %         max(A_normalized)
        %         max(Aerror(Aerror~=Inf))
        %         if  max(Aerror(Aerror~=Inf))  > 0.25*max(A_normalized) ||  max(Berror(Berror~=Inf))  > 0.25*max(B_normalized)
        %             A_normalized = log10(A_normalized);
        %             B_normalized = log10(B_normalized);
        %             Aerror = log10(Aerror);
        %             Berror = log10(Berror);
        %
        %         end
        
        bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
        [XB,YB] = stairs(xlabels,B_normalized');
        XB = [XB;XB(end)+bin_size];
        YB = [YB;YB(end)];
        plot(XB,YB,'linewidth',3.5);
        
        %         stairs(xlabels,B_normalized','linewidth',3.5);
        
        if A(end) ~= 0 || B(end) ~= 0
            xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
            xlim([binRange(1) binRange(end)+(binRange(3)-binRange(2))])
        else
            xticks(binRange)
            xlim([binRange(1) binRange(end)])
        end
        
        %Labels
        if custom_bins == 1 && xlabels(1) ~= 0 && (max(data1) >= binRange(end) || max(data2) >= binRange(end) )
            xlabelCell = xticklabels;
            xlabelCell(1) = strcat('\leq',xlabelCell(1));
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        elseif custom_bins == 1 && (max(data1) >= binRange(end) || max(data2) >= binRange(end) )
            xlabelCell = xticklabels;
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        end
        
        %Find best placement for legend
        [maxA,Ia] = max(A_normalized(A_normalized~=inf));
        [maxB,Ib] = max(B_normalized(B_normalized~=inf));
        if maxA > maxB
            I=Ia;
        else
            I = Ib;
        end
        
        if I > length(binRange)/2
            legend(label1,label2,'FontSize',9,'Location','Northwest')
        else
            legend(label1,label2,'FontSize',9,'Location','Northeast')
        end
        
        hold on
        
        
        errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
        
        errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        
        
        hold off
    else
        A = histcounts(data1,[binRange Inf]);
        A_normalized = A./sum(A);
        Aerror = A.^(-1/2)./sum(A);
        xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        bar(xlabels+bin_size/2,A_normalized',1,'edgecolor','none'); hold on
        errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','black','LineWidth',1.5,'HandleVisibility','off')
        xticks(binRange)
        
        %Labels
        if custom_bins == 1 && xlabels(1) ~= 0
            xlabelCell = xticklabels;
            xlabelCell(1) = strcat('\leq',xlabelCell(1));
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        elseif custom_bins == 1
            xlabelCell = xticklabels;
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        end
        
        %Legend
        %Location of Max
        [~,I] = max(A(A~=inf));
        if I > length(binRange)/2
            legend(label1,'FontSize',9,'Location','Northwest')
        else
            legend(label1,'FontSize',9,'Location','Northeast')
        end
        hold on
        hold off
        
    end
    if contains(parameter(1),'E-CS Magnitude') & length(C) ~=1
        %%%%%%set(gca,'YScale','log')
    elseif contains(parameter(1),'Solar Wind Speed') & length(C) ~=1
        set(gca,'YScale','log')
    elseif contains(parameter(1),'Shock Angle') & length(C) ~=1
        set(gca,'YScale','log')
        
    end
    
    ylimits = ylim;
    if ylimits(1) < 0
        ylim([0 ylimits(2)])
        %         xlim([binRange(1) binRange(end)])
    end
    
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1);
    colormap('winter');
    xlabel(parameter,'FontSize',11)
    ylabel(yAxisLabel,'FontSize',11)
    title(strcat('f vs.', {' '}, parameter(1)),'FontSize',12,'FontWeight', 'normal')
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    
    print(gcf,'-dpng','-r300', '-loose', strcat(fileName));
    savefig(strcat(fileName));
    
end

function [] = plot_stackedHistogram(parameter,data1,binSize1,label1,data2,binSize2,label2)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 350 350])
    set(gcf,'color','w');
    if length(binSize1)==1
        binRange1 = [floor(min(data1)):binSize1:ceil(max(data1)),inf];
    else
        binRange1=[binSize1, inf];
    end
    if length(binSize2) == 1
        binRange2 = floor(min(data2)):binSize2:ceil(max(data2));
    else
        binRange2 = binSize2;
    end
    
    x = binRange1(1:end-1);
    y = [];
    for i=1:length(binRange1)-1
        indicesInBin = find(data1 >= binRange1(i) & data1 < binRange1(i+1));
        y = [y; histcounts(data2(indicesInBin),[binRange2 inf])];
        
    end
    
    bar(x+binSize1/2,y,'stacked')
    xticks(binRange1)
    set(gca,'YScale','log')
    ybinLabels = {};
    for i=1:length(binRange2)
        
        if i==length(binRange2)
            ybinLabels{i} = strcat(num2str(binRange2(i)),{'+ '},'Re');
        else
            ybinLabels{i} = strcat(num2str(binRange2(i)),{' '},'to',{' '},num2str(binRange2(i+1)),{' '},'Re');
        end
    end
    legend(string(ybinLabels),'FontSize',9,'Location','Northeast')
    
    fileName = strcat('stacked_',cell2mat(parameter(1)),'_',label1);
    colormap('winter');
    xlabel(parameter(1:2),'FontSize',11)
    ylabel(strcat('# of Events per',{' '},parameter(3),{ ' '},'Bin'),'FontSize',11)
    title(strcat('f vs.', {' '}, parameter(1),{' '},'and',{' '},parameter(3)),'FontSize',12,'FontWeight', 'normal')
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    
    print(gcf,'-dpng','-r300', '-loose', strcat(fileName));
    savefig(strcat(fileName));
    
end

%Plot the scatter, up to 3 data points
function [] = plot_scatter(parameter,data1a,data1b,label1,data2a,data2b,label2,data31,data32,label3,data3cutoff)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 600 400])
    %     co = [0 0 1;
    %         0 1 0;
    %         1 0 0;
    %         0 0 0];
    %     set(gcf,'defaultAxesColorOrder',co)
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
        scatter(data1a,data1b,data31,'filled','Displayname',label1)
        hold on
        scatter(data2a,data2b,data32,'filled','Displayname',label2)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1,'_',label2,'_',label3);
    elseif nargin > 10
        type1_indices_below = find(data31 < data3cutoff);
        type1_indices_above = find(data31 >= data3cutoff);
        type2_indices_below = find(data32 < data3cutoff);
        type2_indices_above = find(data32 >= data3cutoff);
        
        scatter(data1a(type1_indices_below),data1b(type1_indices_below),'o','filled','Displayname',string(strcat(label1, {' '},'below',{' ' },num2str(data3cutoff))))
        ax=gca;
        ax.ColorOrderIndex =1;
        hold on
        scatter(data1a(type1_indices_above),data1b(type1_indices_above),80,'s','linewidth',2,'Displayname',string(strcat(label1, {' '},'above',{' ' },num2str(data3cutoff))))
        ax.ColorOrderIndex =2;
        scatter(data2a(type2_indices_below),data2b(type2_indices_below),80,'o','filled','Displayname',string(strcat(label2, {' '},'below',{' ' },num2str(data3cutoff))))
        ax.ColorOrderIndex =2;
        scatter(data2a(type2_indices_above),data2b(type2_indices_above),80,'s','linewidth',2,'Displayname',string(strcat(label2, {' '},'above',{' ' },num2str(data3cutoff))))
        ax.ColorOrderIndex =2;
        
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1,'_',label2,'_',label3);
        
    end
    
    xlabel(parameter(1),'FontSize',14)
    ylabel(parameter(2),'FontSize',14)
    %title(strcat('Observation Position', { ' ' }, parameter(1),'-',parameter(2)),'FontSize',16,'FontWeight', 'normal')
    if nargin > 7
        title(strcat(parameter(1), { ' ' },'vs.',{ ' ' },parameter(2),{ ' '}, 'vs.',{ ' ' },parameter(3)),'FontSize',16,'FontWeight', 'normal')
    else
        title(strcat(parameter(1), { ' ' },'vs.',{ ' ' },parameter(2)),'FontSize',16,'FontWeight', 'normal')
    end
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
    savefig(strcat(fileName));
end

%Get the data from the Spreadsheet
function [data] = get_data(Event_Type,Event_Parameter)
    Page_number = 1;
    Column_number = 'A1';
    %     Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabaseFebruary26.xlsx';
    Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabase_Typical2.xlsx';
    
    switch Event_Type
        case 'SS'
            Page_number = 2;
        case 'NS'
            Page_number = 3;
        case 'HFA'
            Page_number = 4;
        case 'HFA SS'
            Page_number = 5;
        case 'HFA NS'
            Page_number = 6;
        case 'SHFA'
            Page_number = 7;
        case 'SHFA SS'
            Page_number = 8;
        case 'SHFA NS'
            Page_number = 9;
        case 'FB'
            Page_number = 10;
    end
    
    
    switch Event_Parameter
        case 'Event Date'
            Column_number = 'D';
        case 'Substructure'
            Column_number = 'C';
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
            
            
        case 'upstreamSpeed'
            Column_number = 'CX';
            
        case 'downstreamSpeed'
            Column_number = 'CY';
            
        case 'DistancetoCSBS'
            Column_number = 'CZ';
            
        case 'ClosestDistance'
            Column_number = 'DA';
            
        case 'ShockAngleLeadingPoint'
            Column_number = 'DN';
            
        case 'ShockAngleTrailingPoint'
            Column_number = 'DO';
            
        case 'Event Width'
            Column_number = 'DR';
            
        case 'Event Height'
            Column_number = 'DS';
            
        case 'Event Area'
            Column_number = 'DT';
            
        case 'Opening Angle'
            Column_number = 'DU';
            
        case 'Vpre Mag'
            Column_number = 'DV';
            
        case 'Vpost Mag'
            Column_number = 'DW';
            
            %         case 'Fractional downV'
            %             Column_number = 'EB';
            %
            %         case 'Fractional upV'
            %             Column_number = 'EC';
            
        case 'Bpre'
            Column_number = 'DX';
            
        case 'Bpost'
            Column_number = 'DY';
            
            
            
        case 'Rmag'
            Column_number = 'EB';
            
        case 'Bpre_cs'
            Column_number = 'EI';
            
        case 'Bpost_cs'
            Column_number = 'EJ';
            
            
        case 'downstreamShockAngle'
            Column_number = 'EK';
            
        case 'upstreamShockAngle'
            Column_number = 'EL';
            
            
        case 'Omni B'
            Column_number = 'EM';
        case 'Omni Speed'
            Column_number = 'EN';
        case 'Omni Density'
            Column_number = 'EO';
        case 'Omni Beta'
            Column_number = 'EP';
        case 'Omni Mach Number'
            Column_number = 'EQ';
            
        case 'Boundaries Expansion Size'
            Column_number = 'ER';
            
            
        case 'Edges Speed Difference'
            Column_number = 'EZ';
            
        case 'HFA Expansion Speed SW'
            Column_number = 'FA';
            
        case 'BoundaryExpansionMinus'
            Column_number = 'FB';
            
            
            
        case 'BoundaryExpansionB'
            Column_number = 'FL';
                        
        case 'BoundaryExpansionBminus'
            Column_number = 'FM';
            
        case 'BoundaryExpansionC'
            Column_number = 'FY';
            
        case 'BoundaryExpansionCminus'
            Column_number = 'FZ';
    end
    
    Column_range = strcat(Column_number,'2:',Column_number,'200');
    
    [data,txt] = xlsread(Database_Directory,Page_number,Column_range);
    data = data(~isnan(data));
    if strcmp(Event_Parameter,'Event Date')
        data = txt(txt~="");
    end
    
    if strcmp(Event_Parameter,'ClosestDistance')
        data = data(data ~= -1);
    end
    %A =datetime(type1_date,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS')
end

function [angle090] = angle90(angle)
    x = abs(cosd(angle));
    y = abs(sind(angle));
    angle090 = atand(y./x);
    
    if angle > 90
       angle090 =  180-angle;
    end
end

function [type_coneAngle] = calculate_radialConeAngle(type)
    type_coneAngle = [get_data(type,'Cone Angle Before'),get_data(type,'Cone Angle After')]';
    [~,I1] = max(abs(type_coneAngle - 90),[],1);
    
    type_coneAngle = type_coneAngle(I1+(0:2:2*(length(type_coneAngle)-1)));
end
