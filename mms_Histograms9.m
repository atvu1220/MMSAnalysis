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

type1b = 'HFA SS';
type2b = 'HFA NS';

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
% load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/WIND_2017_2019.mat')
% load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/WIND_2017_2019_TD4.mat')
% load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/MMS_2017_2019_TD8.mat')
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/MMS_2017_2019_TD12.mat')
%load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/OMNI_2017_2019.mat')
%load('/Users/andrewvu/iCloud Drive (Archive) - 1/Research/MATLAB Analysis/MMS_2017_2019_TD10.mat')
A=5;
%% Plots
% %% Substructures (SS) & No Substructures (NS)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plots for Paper
plot_histogram({'MMS Normalized Closest Distance to Bow Shock','[Re]'},1:2:5,...
    get_data(type1,'ClosestDistance'),'Substructures',...
    get_data(type2,'ClosestDistance'),'No Substructures',...
    MMS_TDClosestDistance)

plot_histogram({'MMS Normalized Distance to TD-BS Connection Point','[Re]'},1:2:5,...
    get_data(type1,'DistancetoCSBS'),'Substructures',...
    get_data(type2,'DistancetoCSBS'),'No Substructures',...
    MMS_BSCSDistance(MMS_magneticShear > 30))

plot_histogram({'MMS Normalized TD-BS Angle','[\Theta]'},0:15:90,...
    angle90(get_data(type1,'CS-BS Angle')),'Substructures',...
    angle90(get_data(type2,'CS-BS Angle')),'No Substructures',...
    angle90(MMS_CSBS_Angle(MMS_magneticShear > 30)))

plot_histogram({'Inner Edges Angle','[\theta]'},0:30:180,...
    get_data(type1,'Event Edges Angle'),'Substructures',...
    get_data(type2,'Event Edges Angle'),'No Substructures')

plot_histogram({'Width','[Re]'},0:0.5:3,...
    rmoutliers(get_data(type1,'Event Width'),'thresholdfactor',100),'Substructures',...
    rmoutliers(get_data(type2,'Event Width'),'thresholdfactor',100),'No Substructures') %113

plot_histogram({'Height','[Re]'},0:2:8,...
    rmoutliers(get_data(type1,'Event Height'),'thresholdfactor',100),'Substructures',...
    rmoutliers(get_data(type2,'Event Height'),'thresholdfactor',100),'No Substructures') %112

plot_histogram({'Area','[Re^2]'},0:1:5,...
    rmoutliers(get_data(type1,'Event Area'),'thresholdfactor',100),'Substructures',...
    rmoutliers(get_data(type2,'Event Area'),'thresholdfactor',100),'No Substructures') %101

plot_histogram({'MMS Normalized Magnetic Shear Angle','[\theta]'},0:20:90,...
    get_data(type1,'Shear Angle'),'Substructures',...
    get_data(type2,'Shear Angle'),'No Substructures',...
    MMS_magneticShear)

plot_histogram({'MMS Normalized Shock Geometry','[\theta]'},0:15:90,...
    angle90(get_data(type1,'downstreamShockAngle')) average_downup()  ,'Substructures',...
    angle90(get_data(type2,'upstreamShockAngle')),'Upstream',...
    angle90(MMS_Shock_Down_Angle))

plot_histogram({'MMS Normalized E-TD Angle','[\theta]'},90,...
    get_data('HFA','E-CS Angle Before'),'Downstream',...
    180-get_data('HFA','E-CS Angle After'),'Upstream',...
    [MMS_E_Down_Angle(MMS_magneticShear > 30);180-MMS_E_Up_Angle(MMS_magneticShear > 30)])

plot_histogram({'MMS Normalized Cone Angle','[\theta]'},30,...
    get_data('all','Cone Angle Before'),'Downstream',...
    get_data('all','Cone Angle After'),'Upstream',...
    MMS_ConeAngle)

plot_histogram({'MMS Normalized Alfvén Mach Number','[M_A]'},8:4:16,...
    get_data('all','Mach Number Before'),'Downstream',...
    get_data('all','Mach Number After'),'Upstream',...
    MMS_MachNumber)

plot_histogram({'MMS Normalized Solar Wind Speed','[km/s]'},[350:50:550],...
    get_data('all','downstreamSpeed'),'Downstream',...
    get_data('all','upstreamSpeed'),'Upstream',...
    MMS_SWSpeed)

plot_histogram({'Velocity Deflection','[\theta]'},0:15:90,...
    get_data('all','CoreSWVelocityDeflection'),'',[],'')

plot_histogram({'MMS Normalized Solar Wind Density','[cm^{-3}]'},[4:2:8],...
    get_data('all','preN'),'Downstream',...
    get_data('all','postN'),'Upstream',...
    (MMS_Nion))

plot_histogram({'MMS Normalized Magnetic Field Strength','[nT]'},3:1.5:6,...
    get_data('all','Bpre_cs'),'Downstream',...
    get_data('all','Bpost_cs'),'Upstream',...
    MMS_BMag)

BE1 = get_data('all','BoundaryExpansionC');
BE2 = get_data(type2,'BoundaryExpansionC');

%Expansion Speed
plot_histogram({'Inner Edges Expansion','[km/s]'},-150:150:600,...
    BE1,'',[],'') %98 Above 0
plot_histogram({'MMS Normalized Transit Speed','[km/s]'},0:150:450,...
    get_data('HFA','Transversal Speed'),'',[],'',...
    MMS_Transversal_Speed(MMS_magneticShear > 30))

%Durations
plot_histogram({'Core Duration','[s]'},0:15:90,...
    get_data('all','Core Duration'),'',[],'')

plot_histogram({'Total Duration','[s]'},0:15:90,...
    get_data('all','Duration'),'',[],'')

%Age
plot_histogram({'Age','[s]'},0:20:120,...
    rmoutliers(get_data('all','Age'),'thresholdfactor',100),'',[],'')

%Size
plot_histogram({'Size','[Re]'},0:1:4,...
    (get_data('HFA','Size')),'',[],'')

%Substructures
plot_histogram({'Substructure Density Ratio','[#]'},-0.25:0.5:5,...
    (get_data(type1,'SSCoreNRatio')),'Core',...
    (get_data(type1,'SSSWNRatio')),'SW')

plot_histogram({'Substructure Velocity Deflection','[\theta]'},0:15:90,...
    get_data(type1,'SSCoreDeflection'),'Core',...
    get_data(type1,'SSSWDeflection'),'SW')

plot_histogram({'Substructure Vx Ratio','[#]'},-0.25:0.5:4,...
    (get_data(type1,'SSCoreVxRatio')),'Core',...
    (get_data(type1,'SSSWVxRatio')),'SW')

plot_histogram({'Substructure Dynamic Pressure Ratio','[#]'},-0.5:1:6,...
    (get_data(type1,'SSCoreDPRatio')),'Core',...
    (get_data(type1,'SSSWDPRatio')),'SW')

plot_histogram({'Substructure Electron Temperature Ratio','[#]'},-0.25:0.5:4,...
    (get_data(type1,'maxTecoreratio')),'Core',...
    (get_data(type1,'maxTeSWRatio')),'SW')

plot_histogram({'Substructure n T Correlation','[#]'},-1.0:0.5:1.0,...
    get_data(type1,'SSenTRatio'),'Electron',...
    get_data(type1,'SSinTRatio'),'Ion')

plot_histogram({'Substructure Magnetic Field Strength Ratio','[#]'},-0.25:0.5:4,...
    (get_data(type1,'SSCoreBRatio')),'Core',...
    (get_data(type1,'SSSWBRatio')),'SW')

plot_histogram({'Substructure Correlation','[#]'},-1:0.5:1,...
    get_data(type1,'SSNVCorr'),'nV',...
    get_data(type1,'SSNBCorr'),'nB')

plot_histogram({'Substructure Duration','[s]'},0.0:1:10,...
    (get_data(type1,'SS Duration')),'Substructured Events',[],'')

plot_histogram({'Substructure Duration Ratio','[Total Core Duration]'},0:0.05:0.4,...
    get_data(type1,'SS Core Duration Ratio'),'Substructured Events',[],'')

plot_histogram({'Substructure Location in Event Core','[Position in Core]'},0:0.1:1,...
    get_data(type1,'SSStartFraction'),'Substructured Events',[],'')

plot_histogram({'Substructure Max','[\sigma]'},0:1:8,...
    get_data('all','MaxNSigma'),'Substructured Events',[],'')

plot_histogram({'Substructure Size','[km]'},0:400:2000,...
    (get_data(type1,'SSSizeBulkV')),'Substructured Events',[],'')

plot_histogram({'Substructure Size in Ion Inertial Lengths','[\lambda_{i}]'},0:0.1:1,...
    get_data(type1,'SSSizeinIonInertialLengths'),'Substructured Events',[],'')

pause
plot_histogram({'Substructure Size in Ion Gyroradii','[R_{g,i}]'},0:1:5,...
    get_data(type1,'SSSizeBulkVinIonGyro'),'Substructured Events',[],'')
plot_histogram({'Substructure Size in Ion Inertial Lengths','[\lambda_{i}]'},0:0.1:1,...
    get_data(type1,'SSSizeinIonInertialLengths'),'Substructured Events',[],'')

plot_histogram({'SS Ion Gyroradius','[km]'},0:100:1000,...
    get_data('all','ionGyroradius'),'Substructured Events',[],'')
plot_histogram({'SS Ion Inertial Length','[km]'},0:25:300,...
    get_data('all','ionInertiallength'),'Substructured Events',[],'')


plot_histogram({'Substructure Size in Electron Gyroradii','[R_{g,e}]'},0:100:1000,...
    get_data('all','SSSizeBulkVinElectronGyro'),'Substructured Events',[],'')
plot_histogram({'Substructure Size in Electron Inertial Lengths','[\lambda_e]'},0:4:40,...
    get_data('all','SSSizeinElectronInertialLengths'),'Substructured Events',[],'')

plot_histogram({'SS Electron Gyroradius','[km]'},0:1:10,...
    get_data('all','electronGyroradius'),'Substructured Events',[],'')
plot_histogram({'SS Electron Inertial Length','[km]'},0:1:10,...
    get_data('all','electronInertiallength'),'Substructured Events',[],'')








A = quantile(get_data('all','MaxNSigma'),[1/3 2/3])

Allsigma = get_data('all','MaxNSigma');
AllSigmaLower = (Allsigma<A(1));
AllSigmaUpper = (Allsigma>A(end));

HFAsigma = get_data('HFA','MaxNSigma');
HFASigmaLower = (HFAsigma<A(1));
HFASigmaUpper = (HFAsigma>A(end));

plot_histogram({'Width Compare','[Re]'},0:0.5:3,...
    rmoutliers(get_data('all','Event Width',AllSigmaLower),'thresholdfactor',100),'Lower',...
    rmoutliers(get_data('all','Event Width',AllSigmaUpper),'thresholdfactor',100),'Upper') %113

plot_histogram({'Core Duration Compare','[s]'},0:15:90,...
    get_data('all','Core Duration',AllSigmaLower),'Lower',...
    get_data('all','Core Duration',AllSigmaUpper),'Upper')  

plot_histogram({'Height','[Re]'},0:2:8,...
    rmoutliers(get_data('all','Event Height',AllSigmaLower),'thresholdfactor',100),'Lower',...
    rmoutliers(get_data('all','Event Height',AllSigmaUpper),'thresholdfactor',100),'Upper') %112

plot_histogram({'Size','[Re]'},0:1:4,...
    get_data('HFA','Size',HFASigmaLower),'Lower',...
    get_data('HFA','Size',HFASigmaUpper),'Upper')

%% Bow Shock Parameters

%Closest Distance to BS
% plot_histogram({'MMS Closest Distance to Bow Shock','[Re]'},-2:2:8,...
%     MMS_TDClosestDistance,'MMS')
%
% plot_histogram({'Closest Distance to Bow Shock','[Re]'},1:2:6,...
%     get_data(type1,'ClosestDistance'),name1,...
%     get_data(type2,'ClosestDistance'),name2)

% % plot_histogram({'MMS Normalized Closest Distance to Bow Shock','[Re]'},1:2:5,...
% %     get_data(type1,'ClosestDistance'),name1,...
% %     get_data(type2,'ClosestDistance'),name2,...
% %     MMS_TDClosestDistance)




% plot_histogram({'MMS Normalized Closest Distance to Bow Shock','[Re]'},0:2:8,...
%     get_data(type1,'ClosestDistance'),name1,...
%     get_data(type2,'ClosestDistance'),name2,...
%     MMS_DistanceToBS)


%Distance to CS-BS Connection Point
% plot_histogram({'MMS Distance to TD-BS Connection Point','[Re]'},-2:2:8,...
%     MMS_BSCSDistance(MMS_magneticShear > 30),'MMS')
%
% plot_histogram({'Distance to TD-BS Connection Point','[Re]'},1:2:6,...
%     get_data(type1b,'DistancetoCSBS'),name1,...
%     get_data(type2b,'DistancetoCSBS'),name2)

% % plot_histogram({'MMS Normalized Distance to TD-BS Connection Point','[Re]'},1:2:5,...
% %     get_data(type1b,'DistancetoCSBS'),name1,...
% %     get_data(type2b,'DistancetoCSBS'),name2,...
% %     MMS_BSCSDistance(MMS_magneticShear > 30))


% %CS-BS Angle
% plot_histogram({'MMS TD-BS Angle','[\theta]'},0:15:90,...
%     angle90(MMS_CSBS_Angle(MMS_magneticShear > 30)),'MMS')
%
% plot_histogram({'TD-BS Angle','[\Theta]'},0:15:90,...
%     angle90(get_data(type1b,'CS-BS Angle')),name1,...
%     angle90(get_data(type2b,'CS-BS Angle')),name2)

% % plot_histogram({'MMS Normalized TD-BS Angle','[\Theta]'},0:15:90,...
% %     angle90(get_data(type1b,'CS-BS Angle')),name1,...
% %     angle90(get_data(type2b,'CS-BS Angle')),name2,...
% %     angle90(MMS_CSBS_Angle(MMS_magneticShear > 30)))



%% Event Dimensions
%Event Width, Length, Area

%This is angle of the intersection point
% plot_histogram({'Event Edges Angle','[\Theta]'},0:30:180,...
%     get_data(type1,'Opening Angle'),name1,...
%     get_data(type2,'Opening Angle'),name2)

%This is the angle between event edges
% % plot_histogram({'Inner Edges Angle','[\theta]'},0:45:180,...
% %     get_data(type1,'Event Edges Angle'),name1,...
% %     get_data(type2,'Event Edges Angle'),name2)
% %
% % plot_histogram({'Width','[Re]'},0:1:4,...
% %     (get_data(type1,'Event Width')),name1,...
% %     (get_data(type2,'Event Width')),name2)



% plot_histogram({'Event Height','[Re]'},0:2:8,...
%     get_data(type1,'Event Height')-get_data(type1,'ClosestDistance'),name1,...
%     get_data(type2,'Event Height')-get_data(type2,'ClosestDistance'),name2)

% % plot_histogram({'Height','[Re]'},0:2:8,...
% %     (get_data(type1,'Event Height')),name1,...
% %     (get_data(type2,'Event Height')),name2)

% plot_histogram({'Closest Distance to Bow Shock','[Re]'},1:2:6,...
%     ,name1,...
%     ,name2)

% % plot_histogram({'Area','[Re^2]'},0:1:5,...
% %     rmoutliers(get_data(type1,'Event Area'),'thresholdfactor',5),name1,...
% %     rmoutliers(get_data(type2,'Event Area'),'thresholdfactor',5),name2)
% % %
% % plot_histogram({'Area','[Re^2]'},0:1:8,...
% %     get_data(type1,'Event Area'),name1,...
% %     get_data(type2,'Event Area'),name2)
% %
% % % plot_histogram({'Area','[Re^2]'},0:1:20,...
% % %     0.5.*(get_data(type1,'Event Width')).*(get_data(type1,'Event Height')),name1,...
% % %     0.5.*(get_data(type2,'Event Width')).*(get_data(type2,'Event Height')),name2)
% %
% % plot_scatter({'Event Width','Event Height'},...
% %     (get_data(type1,'Event Width')),(get_data(type1,'Event Height')),name1,...
% %     (get_data(type2,'Event Width')),(get_data(type2,'Event Height')),name2)
% %
% % A = get_data(type1,'V1');
% % Ad = get_data(type1,'Core Duration');
% % Asize = A.*Ad./6371;
% % B = get_data(type2,'V1');
% % Bd = get_data(type2,'Core Duration');
% % Bsize = B.*Bd./6371;
% % plot_histogram({'Size Test in SW Frame','[Re]'},0:1:5,...
% %     abs(Asize),name1,...
% %     abs(Bsize),name2)
% %
% % A = get_data(type1,'Vt1');
% % Ad = get_data(type1,'Core Duration');
% % Asize = A.*Ad./6371;
% % B = get_data(type2,'Vt1');
% % Bd = get_data(type2,'Core Duration');
% % Bsize = B.*Bd./6371;
% % plot_histogram({'Size Test in SC Frame','[Re]'},0:1:5,...
% %     abs(Asize),name1,...
% %     abs(Bsize),name2)
% %
% % A = get_data(type1,'BoundaryExpansionC');
% % Ad = get_data(type1,'Core Duration');
% % Asize = A.*Ad./6371;
% % B = get_data(type2,'BoundaryExpansionC');
% % Bd = get_data(type2,'Core Duration');
% % Bsize = B.*Bd./6371;
% % plot_histogram({'Size Test Using Boundaries Expansion','[Re]'},0:1:5,...
% %     abs(Asize),name1,...
% %     abs(Bsize),name2)



% plot_scatter({'Event Width','Event Height'},...
%     get_data(type1,'Event Width'),get_data(type1,'Event Height'),'Substructure',...
%     get_data(type2,'Event Width'),get_data(type2,'Event Height'),'No Substructure')

%% Shear Angle

% plot_histogram({'MMS Magnetic Shear Angle','[\theta]'},15,...
%     MMS_magneticShear,'MMS')
%
% plot_histogram({'Magnetic Shear Angle','[\theta]'},0:20:90,...
%     get_data(type1,'Shear Angle'),name1,...
%     get_data(type2,'Shear Angle'),name2)

% % plot_histogram({'MMS Normalized Magnetic Shear Angle','[\theta]'},0:20:90,...
% %     get_data(type1,'Shear Angle'),name1,...
% %     get_data(type2,'Shear Angle'),name2,...
% %     MMS_magneticShear)

%% Shock Angles
%
% plot_histogram({'MMS Shock Geometry Downstream','[\theta]'},0:15:90,...
%     angle90(MMS_Shock_Down_Angle),'MMS')
%
% plot_histogram({'Shock Geometry Downstream','[\theta]'},0:15:90,...
%     angle90(get_data(type1,'downstreamShockAngle')),name1,...
%     angle90(get_data(type2,'downstreamShockAngle')),name2)
%
% % plot_histogram({'MMS Normalized Shock Geometry Downstream','[\theta]'},0:22.5:90,...
% %     angle90(get_data(type1,'downstreamShockAngle')),name1,...
% %     angle90(get_data(type2,'downstreamShockAngle')),name2,...
% %     angle90(MMS_Shock_Down_Angle))
% % %
% % %
% % % plot_histogram({'MMS Shock Geometry Upstream','[\theta]'},0:15:90,...
% % %     angle90(MMS_Shock_Up_Angle),'MMS')
% % %
% % % plot_histogram({'Shock Geometry Upstream','[\theta]'},0:15:90,...
% % %     angle90(get_data(type1,'upstreamShockAngle')),name1,...
% % %     angle90(get_data(type2,'upstreamShockAngle')),name2)
% % %
% % plot_histogram({'MMS Normalized Shock Geometry Upstream','[\theta]'},0:22.5:90,...
% %     angle90(get_data(type1,'upstreamShockAngle')),name1,...
% %     angle90(get_data(type2,'upstreamShockAngle')),name2,...
% %     angle90(MMS_Shock_Up_Angle))
%
% plot_scatter({'Shock Geometry Downstream','Shock Geometry Upstream'},...
%     angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
%     angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2)
% xline(45,'--k')
% yline(45,'--k')
%
% plot_scatter({'Shock Geometry Downstream','Shock Geometry Upstream','Size'},...
%     angle90(get_data(type1,'Shock Angle Before')),angle90(get_data(type1,'Shock Angle After')),name1,...
%     angle90(get_data(type2,'Shock Angle Before')),angle90(get_data(type2,'Shock Angle After')),name2,...
%     150*get_data(type1,'Size'),150*get_data(type2,'Size'),'Size')

%% E-CS Orientation

%E-CS Before Mag
% plot_histogram({'E-CS Magnitude Downstream WIND','[nT km/s]'},-400:50:400,...
%     E_Down_Mag,'WIND')

% plot_histogram({'MMS E-TD Magnitude Downstream','[nT km/s]'},-1500:500:1500,...
%     MMS_E_Down_Mag(MMS_magneticShear > 30),'MMS')
%
% plot_histogram({'E-TD Magnitude Downstream','[nT km/s]'},-1000:500:1000,...
%     get_data(type1b,'E-CS Magnitude Before'),name1,...
%     get_data(type2b,'E-CS Magnitude Before'),name2)

% plot_histogram({'MMS Normalized E-TD Magnitude Downstream','[nT km/s]'},-1000:500:1000,...
%     get_data(type1b,'E-CS Magnitude Before'),name1,...
%     get_data(type2b,'E-CS Magnitude Before'),name2,...
%     MMS_E_Down_Mag(MMS_magneticShear > 30))


%E-CS Before Angle
% plot_histogram({'MMS E-TD Angle Downstream','[\theta]'},15,...
%     MMS_E_Down_Angle(MMS_magneticShear > 30),'MMS')
%
% plot_histogram({'E-TD Angle Downstream','[\theta]'},90,...
%     get_data(type1b,'E-CS Angle Before'),name1,...
%     get_data(type2b,'E-CS Angle Before'),name2)

% % plot_histogram({'MMS Normalized E-TD Angle Downstream','[\theta]'},90,...
% %     get_data(type1b,'E-CS Angle Before'),name1,...
% %     get_data(type2b,'E-CS Angle Before'),name2,...
% %     MMS_E_Down_Angle(MMS_magneticShear > 30))

%E-CS After Mag
% plot_histogram({'MMS E-TD Magnitude Upstream','[nT km/s]'},-1500:500:1500,...
%     -MMS_E_Up_Mag(MMS_magneticShear > 30),'MMS')
%
% plot_histogram({'E-TD Magnitude Upstream','[nT km/s]'},-1000:500:1000,...
%     -get_data(type1b,'E-CS Magnitude After'),name1,...
%     -get_data(type2b,'E-CS Magnitude After'),name2)
%
% plot_histogram({'MMS Normalized E-TD Magnitude Upstream','[nT km/s]'},-1000:500:1000,...
%     -get_data(type1b,'E-CS Magnitude After'),name1,...
%     -get_data(type2b,'E-CS Magnitude After'),name2,...
%     -MMS_E_Up_Mag(MMS_magneticShear > 30))

%E-CS After Angle
% plot_histogram({'MMS E-TD Angle Upstream','[\theta]'},15,...
%     180-MMS_E_Up_Angle,'MMS')
%
% plot_histogram({'E-TD Angle Upstream','[\theta]'},90,...
%     180-get_data(type1b,'E-CS Angle After'),name1,...
%     180-get_data(type2b,'E-CS Angle After'),name2)

% % plot_histogram({'MMS Normalized E-TD Angle Upstream','[\theta]'},90,...
% %     180-get_data(type1b,'E-CS Angle After'),name1,...
% %     180-get_data(type2b,'E-CS Angle After'),name2,...
% %     180-MMS_E_Up_Angle(MMS_magneticShear > 30))


%
% plot_scatter({'E-CS Angle Downstream','E-CS Angle Upstream','Event Area'},...
%     get_data(type1,'E-CS Angle Before'),180-get_data(type1,'E-CS Angle After'),name1,...
%     get_data(type2,'E-CS Angle Before'),180-get_data(type2,'E-CS Angle After'),name2,...
%     40*get_data(type1,'Event Area'),40*get_data(type2,'Event Area'),'Event Area')

% plot_scatter({'E-TD Angle Downstream','E-TD Angle Upstream','Size'},...
%     get_data(type1b,'E-CS Angle Before'),180-get_data(type1b,'E-CS Angle After'),name1,...
%     get_data(type2b,'E-CS Angle Before'),180-get_data(type2b,'E-CS Angle After'),name2,...
%     150*get_data(type1b,'Size'),150*get_data(type2b,'Size'),'Size')

% % plot_scatter({'E-TD Angle Downstream','E-TD Angle Upstream'},...
% %     get_data(type1b,'E-CS Angle Before'),180-get_data(type1b,'E-CS Angle After'),name1,...
% %     get_data(type2b,'E-CS Angle Before'),180-get_data(type2b,'E-CS Angle After'),name2)

%% Cone Angle
%Choosing Cone Angle that is closest to 0 or 180.
% plot_histogram({'MMS Cone Angle','[\theta]'},30,...
%     MMS_ConeAngle,'MMS')
%
% plot_histogram({'Cone Angle','[\theta]'},30,...
%     calculate_radialConeAngle(type1)',name1,...
%     calculate_radialConeAngle(type2)',name2)

% % plot_histogram({'MMS Normalized Cone Angle','[\theta]'},30,...
% %     calculate_radialConeAngle(type1)',name1,...
% %     calculate_radialConeAngle(type2)',name2,...
% %     MMS_ConeAngle)

% [get_data(type,'Cone Angle Before'),get_data(type,'Cone Angle After')]

% plot_histogram({'Cone Angle Downstream','[\theta]'},30,...
%     get_data(type1,'Cone Angle Before'),name1,...
%     get_data(type2,'Cone Angle Before'),name2)

% % plot_histogram({'MMS Normalized Cone Angle Downstream','[\theta]'},30,...
% %     get_data(type1,'Cone Angle Before'),name1,...
% %     get_data(type2,'Cone Angle Before'),name2,...
% %     MMS_ConeAngle)


% plot_histogram({'Cone Angle Upstream','[\theta]'},30,...
%     get_data(type1,'Cone Angle After'),name1,...
%     get_data(type2,'Cone Angle After'),name2)

% % plot_histogram({'MMS Normalized Cone Angle Upstream','[\theta]'},30,...
% %     get_data(type1,'Cone Angle After'),name1,...
% %     get_data(type2,'Cone Angle After'),name2,...
% %     MMS_ConeAngle)

%% Mach Numbers
% plot_histogram({'Alfvén Mach Number','[#]'},0:2:20,...
%     MMS_MachNumber,'MMS')
%
% plot_histogram({'Alfvén Mach Number Downstream','[#]'},8:4:16,...
%     get_data(type1,'Mach Number Before'),name1,...
%     get_data(type2,'Mach Number Before'),name2)
%
% plot_histogram({'Alfvén Mach Number Upstream','[#]'},8:4:16,...
%     get_data(type1,'Mach Number After'),name1,...
%     get_data(type2,'Mach Number After'),name2)

% % plot_histogram({'MMS Normalized Alfvén Mach Number Downstream','[M_A]'},8:4:16,...
% %     get_data(type1,'Mach Number Before'),name1,...
% %     get_data(type2,'Mach Number Before'),name2, MMS_MachNumber)
% %
% % plot_histogram({'MMS Normalized Alfvén Mach Number Upstream','[M_A]'},8:4:16,...
% %     get_data(type1,'Mach Number After'),name1,...
% %     get_data(type2,'Mach Number After'),name2, MMS_MachNumber)

%% Solar Wind Speed
% plot_histogram({'MMS Solar Wind Speed','[km/s]'},[200:50:700],...
%     MMS_SWSpeed,'MMS')
%
% plot_histogram({'Solar Wind Speed Downstream','[km/s]'},[400:50:600],...
%     get_data(type1,'downstreamSpeed'),name1,...
%     get_data(type2,'downstreamSpeed'),name2)
%
% plot_histogram({'Solar Wind Speed Upstream','[km/s]'},[400:50:600],...
%     get_data(type1,'upstreamSpeed'),name1,...
%     get_data(type2,'upstreamSpeed'),name2)

% % plot_histogram({'MMS Normalized Solar Wind Speed Downstream','[km/s]'},[350:50:550],...
% %     get_data(type1,'downstreamSpeed'),name1,...
% %     get_data(type2,'downstreamSpeed'),name2,...
% %     MMS_SWSpeed)
% %
% % plot_histogram({'MMS Normalized Solar Wind Speed Upstream','[km/s]'},[350:50:550],...
% %     get_data(type1,'upstreamSpeed'),name1,...
% %     get_data(type2,'upstreamSpeed'),name2,...
% %     MMS_SWSpeed)
%% Velocity Deflection Angle
% plot_histogram({'MMS Velocity Deflection','[\theta]'},[0:15:150],...
%     MMS_V_DeflectionAngle*180/pi,'MMS')

% % plot_histogram({'Velocity Deflection','[\theta]'},0:15:90,...
% %     get_data(type1,'CoreSWVelocityDeflection'),name1,...
% %     get_data(type2,'CoreSWVelocityDeflection'),name2)

% plot_histogram({'MMS Normalized Velocity Deflection','[\theta]'},[0:15:150],...
%     get_data(type1,'downstreamSpeed'),name1,...
%     get_data(type2,'downstreamSpeed'),name2,...
%     MMS_V_DeflectionAngle)

%% Solar Wind Ion Density

% plot_histogram({'MMS Solar Wind Density','[cm^{-3}]'},[0:3:12],...
%     MMS_Nion,'MMS')
%
% plot_histogram({'Ion Density Downstream','[cm^{-3}]'},[0:3:12],...
%     get_data(type1,'downstream N'),name1,...
%     get_data(type2,'downstream N'),name2)
%
% plot_histogram({'Ion Density Upstream','[cm^{-3}]'},[0:3:12],...
%     get_data(type1,'upstream N'),name1,...
%     get_data(type2,'upstream N'),name2)

% % plot_histogram({'MMS Normalized Ion Density Downstream','[cm^{-3}]'},[4:2:10],...
% %     get_data(type1,'downstream N'),name1,...
% %     get_data(type2,'downstream N'),name2,...
% %     (MMS_Nion))
% %
% % plot_histogram({'MMS Normalized Ion Density Upstream','[cm^{-3}]'},[4:2:10],...
% %     get_data(type1,'upstream N'),name1,...
% %     get_data(type2,'upstream N'),name2,...
% %     (MMS_Nion))
% %
% %
% % plot_histogram({'MMS Normalized Solar Wind Density Downstream','[cm^{-3}]'},[4:2:8],...
% %     get_data(type1,'preN'),name1,...
% %     get_data(type2,'preN'),name2,...
% %     (MMS_Nion))
% %
% % plot_histogram({'MMS Normalized Solar Wind Density Upstream','[cm^{-3}]'},[4:2:8],...
% %     get_data(type1,'postN'),name1,...
% %     get_data(type2,'postN'),name2,...
% %     (MMS_Nion))


% A = (get_data(type1,'preN')-get_data(type1,'downstream N'))./get_data(type1,'downstream N');
% B = (get_data(type1,'postN')-get_data(type1,'upstream N'))./get_data(type1,'upstream N');
% C = (get_data(type2,'preN')-get_data(type2,'downstream N'))./get_data(type2,'downstream N');
% D = (get_data(type2,'postN')-get_data(type2,'upstream N'))./get_data(type2,'upstream N');
%
% % JJ = [get_data(type1,'preN'),get_data(type1,'downstream N'),get_data(type1,'preN')-get_data(type1,'downstream N')]
% % KK = [get_data(type1,'postN'),get_data(type1,'upstream N'),get_data(type1,'postN')-get_data(type1,'upstream N')]
%
% plot_histogram({'Foreshock Density Fraction Downstream','[#]'},[0:0.025:0.2],...
%     A(~(A<0)),name1,...
%     C(~(C<0)),name2)
%
% plot_histogram({'Foreshock Density Fraction Upstream','[#]'},[0:0.025:0.2],...
%     B(~(B<0)),name1,...
%     D(~(D<0)),name2)

%% Magnetic Field Magnitude
% plot_histogram({'MMS Magnetic Field Strength','[nT]'},0:3:12,...
%     MMS_BMag,'MMS')

% plot_histogram({'Magnetic Field Downstream','[nT]'},3:1.5:6,...
%     get_data(type1,'Bpre'),name1,get_data(type2,'Bpre'),name2)
%
% plot_histogram({'Magnetic Field Upstream','[nT]'},3:1.5:6,...
%     get_data(type1,'Bpost'),name1,get_data(type2,'Bpost'),name2)
%
% plot_histogram({'Normalized Magnetic Field Downstream','[nT]'},3:1.5:6,...
%     get_data(type1,'Bpre'),name1,get_data(type2,'Bpre'),name2,Bdatamag)
%
% plot_histogram({'Normalized Magnetic Field Upstream','[nT]'},3:1.5:6,...
%     get_data(type1,'Bpost'),name1,get_data(type2,'Bpost'),name2,Bdatamag)


% plot_histogram({'Magnetic Field Strength Downstream','[nT]'},3:1.5:6,...
%     get_data(type1,'Bpre_cs'),name1,...
%     get_data(type2,'Bpre_cs'),name2)
%
% plot_histogram({'Magnetic Field Strength Upstream','[nT]'},3:1.5:6,...
%     get_data(type1,'Bpost_cs'),name1,...
%     get_data(type2,'Bpost_cs'),name2)

% % plot_histogram({'MMS Normalized Magnetic Field Strength Downstream','[nT]'},3:1.5:6,...
% %     get_data(type1,'Bpre_cs'),name1,...
% %     get_data(type2,'Bpre_cs'),name2,...
% %     MMS_BMag)
% %
% % plot_histogram({'MMS Normalized Magnetic Field Strength Upstream','[nT]'},3:1.5:6,...
% %     get_data(type1,'Bpost_cs'),name1,...
% %     get_data(type2,'Bpost_cs'),name2,...
% %     MMS_BMag)

%% Size, Expansion Speed, Transversal Speed

% % BE1 = get_data(type1,'BoundaryExpansionC');
% % BE2 = get_data(type2,'BoundaryExpansionC');
% %
% % %Expansion Speed
% % plot_histogram({'Inner Edges Expansion','[km/s]'},-150:200:600,...
% %     BE1(BE1>0),name1,...
% %     BE2(BE2>0),name2)
% %
% % plot_histogram({'Inner Edges Expansion','[km/s]'},-100:200:600,...
% %     BE1,name1,...
% %     BE2,name2)

% plot_histogram({'Boundaries Expansion','[km/s]'},-0:100:500,...
%     rmoutliers(get_data(type1,'BoundaryExpansionC')),name1,...
%     rmoutliers(get_data(type2,'BoundaryExpansionC')),name2)

%Transversal Speed
% plot_histogram({'MMS Transit Speed','[km/s]'},0:150:600,...
%     MMS_Transversal_Speed(MMS_magneticShear > 30),'MMS')
%
% plot_histogram({'Transit Speed','[km/s]'},0:150:450,...
%     get_data(type1b,'Transversal Speed'),name1,...
%     get_data(type2b,'Transversal Speed'),name2)

% % plot_histogram({'MMS Normalized Transit Speed','[km/s]'},0:150:450,...
% %     get_data(type1b,'Transversal Speed'),name1,...
% %     get_data(type2b,'Transversal Speed'),name2,...
% %     MMS_Transversal_Speed(MMS_magneticShear > 30))
% %
% % %Durations
% % plot_histogram({'Core Duration','[s]'},0:15:90,...
% %     get_data(type1,'Core Duration'),name1,...
% %     get_data(type2,'Core Duration'),name2)
% %
% % plot_histogram({'Total Duration','[s]'},0:15:90,...
% %     get_data(type1,'Duration'),name1,...
% %     get_data(type2,'Duration'),name2)
% %
% % %Age
% % % plot_histogram({'Event Age','[s]'},0:20:100,...
% % %     get_data(type1,'Age'),name1,...
% % %     get_data(type2,'Age'),name2)
% %
% % %Size
% % plot_histogram({'Size','[Re]'},0:1:4,...
% %     (get_data(type1,'Size')),name1,...
% %     (get_data(type2,'Size')),name2)

% plot_histogram({'CoreSize1IonLength','[Re]'},0:4:30,...
%     rmoutliers(get_data(type1,'CoreSize1IonLength')),name1)
%
% plot_histogram({'CoreSize2IonLength','[Re]'},0:4:30,...
%     rmoutliers(get_data(type1,'CoreSize2IonLength')),name1)


%Age
% plot_histogram({'Event Age','[s]'},0:15:75,...
%    abs( get_data(type1,'Size')./get_data(type1,'BoundaryExpansionC').*6371.2),name1,...
%    abs( get_data(type2,'Size')./get_data(type2,'BoundaryExpansionC').*6371.2),name2)

%Age
% plot_histogram({'Age','[s]'},0:30:150,...
%     rmoutliers(get_data(type1,'Age')),name1,...
%     rmoutliers(get_data(type2,'Age')),name2)
% 
% plot_histogram({'Event Age 2p','[s]'},0:5:50,...
%     (get_data(type1,'Age2p')),name1,...
%     (get_data(type2,'Age2p')),name2)

% plot_histogram({'Event Age3','[s]'},0:15:300,...
%    get_data(type1,'Age2'),name1,...
%    get_data(type2,'Age2'),name2)





%Size
% type1Size = rmoutliers(get_data(type1,'Size'));
% type2Size = rmoutliers(get_data(type2,'Size'));
% plot_histogram({'Size','[Re]'},0:1:5,...
%     type1Size,name1,...
%     type2Size,name2)

% type1_coreDuration = get_data(type1,'Core Duration');
% type2_coreDuration = get_data(type2,'Core Duration');
%
% type1_boundary1Speed = get_data(type1,'Boundary1Vmag');
% type1_boundary2Speed = get_data(type1,'Boundary2Vmag');
%
% type2_boundary1Speed = get_data(type2,'Boundary1Vmag');
% type2_boundary2Speed = get_data(type2,'Boundary2Vmag');
%
% type1_sizeBoundary1 = abs(type1_coreDuration.*type1_boundary1Speed./6371);
% type1_sizeBoundary2 = abs(type1_coreDuration.*type1_boundary2Speed./6371);
%
% type2_sizeBoundary1 = abs(type2_coreDuration.*type2_boundary1Speed./6371);
% type2_sizeBoundary2 = abs(type2_coreDuration.*type2_boundary2Speed./6371);
%
% plot_histogram({'Size using Boundary1 Speed','[Re]'},0:1:5,...
%     type1_sizeBoundary1,name1,...
%     type2_sizeBoundary1,name2)
%
% plot_histogram({'Size using Boundary2 Speed','[Re]'},0:1:5,...
%     type1_sizeBoundary2,name1,...
%     type2_sizeBoundary2,name2)

% plot_scatter({'Solar Wind Speed Upstream','Boundaries Expansion'},...
%     get_data(type1,'upstreamSpeed'),get_data(type1,'BoundaryExpansionC'),name1,...
%     get_data(type2,'upstreamSpeed'),get_data(type2,'BoundaryExpansionC'),name2)

%% Substructure Analysis

% % plot_histogram({'Substructure Density Ratio','[#]'},-0.25:0.5:5,...
% %     (get_data(type1,'SSCoreNRatio')),'Core',...
% %     (get_data(type1,'SSSWNRatio')),'SW')
% %
% % plot_histogram({'Substructure Velocity Deflection','[\theta]'},0:15:90,...
% %     get_data(type1,'SSCoreDeflection'),'Core',...
% %     get_data(type1,'SSSWDeflection'),'SW')
% %
% % plot_histogram({'Substructure Vx Ratio','[#]'},-0.25:0.5:4,...
% %     (get_data(type1,'SSCoreVxRatio')),'Core',...
% %     (get_data(type1,'SSSWVxRatio')),'SW')
% %
% % plot_scatter({'Substructure Density Ratio','Substructure Vx Ratio'},...
% %     get_data(type1b,'SSSWNRatio'),get_data(type1b,'SSSWVxRatio'),'SW')
% %
% % % plot_histogram({'Substructure Dynamic Pressure Ratio','[#]'},-0.25:0.5:5,...
% % %     (get_data(type1,'SSCoreDPRatio')),'Core',...
% % %     (get_data(type1,'SSSWDPRatio')),'SW')
% % plot_histogram({'Substructure Dynamic Pressure Ratio','[#]'},-0.5:1:6,...
% %     (get_data(type1,'SSCoreDPRatio')),'Core',...
% %     (get_data(type1,'SSSWDPRatio')),'SW')
% %
% % % A = get_data(type1,'SSSWDeflection');
% % %
% % %
% % % B =  get_data(type1,'Rx');
% % % B = B(logical((A > 30).*(A < 45)));
% % % C = get_data(type1,'Ry');
% % % C = C(logical((A > 30).*(A < 45)));
% % %
% % % plot_scatter({'Rx','Ry'},...
% % %     B,C,name1,...
% % %     get_data(type2,'Rx'),get_data(type2,'Ry'),name2)
% %
% %
% % plot_histogram({'Substructure Temperature Ratio','[#]'},0:0.5:4,...
% %     get_data(type1,'SSCoreTRatio'),'Core',...
% %     get_data(type1,'SSSWTRatio'),'SW')
% % plot_histogram({'Substructure Electron Temperature Ratio','[#]'},-0.25:0.5:5,...
% %     (get_data(type1,'maxTecoreratio')),'Core',...
% %     (get_data(type1,'maxTeSWRatio')),'SW')
% %
% % plot_histogram({'Substructure & Core Temperature Ratio','[#]'},-0.25:0.5:2,...
% %     get_data(type1,'SSCoreTRatio'),'Ion',...
% %     get_data(type1,'maxTecoreratio'),'Electron')
% %
% % plot_histogram({'Substructure n T Correlation','[#]'},-1.0:0.5:1.0,...
% %     get_data(type1,'SSenTRatio'),'Electron',...
% %     get_data(type1,'SSinTRatio'),'Ion')
% %
% %
% %
% % plot_histogram({'Substructure Magnetic Field Strength Ratio','[#]'},-0.25:0.5:5,...
% %     (get_data(type1,'SSCoreBRatio')),'Core',...
% %     (get_data(type1,'SSSWBRatio')),'SW')
% %
% % CoreB = get_data('all','SSCoreBRatio');
% % SWB = get_data('all','SSSWBRatio');
% % coreSWB = SWB./CoreB;
% %
% % plot_histogram({'Core Over Ambient Magnetic Field Ratio','[#]'},0:0.5:4,...
% %     coreSWB,'HFAs,FBs')
% %
% % plot_histogram({'Substructure Correlation','[#]'},-1:1/3:1,...
% %     get_data(type1,'SSNVCorr'),'nV',...
% %     get_data(type1,'SSNBCorr'),'nB')
% %
% %
% %
% % % plot_histogram({'Substructure Size','[km]'},0:5:30,...
% % %     (get_data(type1,'SS Size')),'Substructured Events') %Timing Method V
% %
% % % plot_histogram({'Substructure Size in Electron Inertial Lengths','[\lambda_e]'},0:1:9,...
% % %     rmoutliers(get_data(type1,'SS Electron Scale'),'thresholdfactor',1),'Substructured Events')
% %
% % % plot_histogram({'Substructure Speed','[km/s]'},0:5:30,...
% % %     rmoutliers(get_data(type1,'SS Speed'),'thresholdfactor',2),'Substructured Events')
% %
% % % plot_histogram({'Substructure Duration','[s]'},0.0:0.25:2,...
% % %     (get_data(type1,'SS Duration')),'Substructured Events')
% % plot_histogram({'Substructure Duration','[s]'},0.0:1:10,...
% %     (get_data(type1,'SS Duration')),'Substructured Events')
% %
% % %plot_histogram({'Substructure Duration Ratio','[Total Core Duration]'},0:0.01:.1,...
% % %   get_data(type1,'SS Core Duration Ratio'),'Substructured Events')
% % plot_histogram({'Substructure Duration Ratio','[Total Core Duration]'},0:0.05:0.5,...
% %     get_data(type1,'SS Core Duration Ratio'),'Substructured Events')
% %
% % plot_histogram({'Substructure Location in Event Core','[Position in Core]'},0:0.1:1,...
% %     get_data(type1,'SSStartFraction'),'Substructured Events')
% %
% % % plot_histogram({'CoreSSBeginIonLength2','[\lambda_i]'},0:0.5:10,...
% % %     rmoutliers(get_data(type1,'CoreSSBeginIonLength2')),'Substructured Events')
% % %
% % % plot_histogram({'CoreSSBeginIonLength1','[\lambda_i]'},0:0.5:10,...
% % %     rmoutliers(get_data(type1,'CoreSSBeginIonLength1')),'Substructured Events')
% %
% %
% % %
% % % plot_histogram({'SubstructureSW Density Ratio','[#]'},0:1:6,...
% % %     get_data(type1,'SSSWNRatio'),name1)
% %
% % % plot_histogram({'SubstructureSW Velocity Deflection','[\theta]'},0:15:90,...
% % %     get_data(type1,'SSSWDeflection'),name1)
% %
% % % plot_histogram({'SubstructureSW Vx Ratio','[#]'},0:0.5:3,...
% % %     get_data(type1,'SSSWVxRatio'),name1)
% % %
% % % plot_histogram({'SubstructureSW Dynamic Pressure Ratio','[#]'},0:1:8,...
% % %     get_data(type1,'SSSWDPRatio'),name1)
% % %
% % % plot_histogram({'SubstructureSW Temperature Ratio','[#]'},0:0.2:1,...
% % %     get_data(type1,'SSSWTRatio'),name1)
% %
% %
% %
% % plot_histogram({'Substructure Size','[km]'},0:400:2000,...
% %     (get_data(type1,'SSSizeBulkV')),'Substructured Events')
% %
% % plot_histogram({'Substructure Size in Ion Gyroradii','[R_{g,i}]'},0:0.5:8,...
% %     get_data(type1,'SSSizeBulkVinIonGyro'),'Substructured Events')

% plot_histogram({'Substructure Size with Duration And Event Size','[km]'},0:100:800,...
%     rmoutliers(get_data(type1,'SSSizewithDurationRatioandSize')*6371.2),'Substructured Events')
%Some Calculations

N = (get_data(type1,'SSSWNRatio'));
Vx = (get_data(type1,'SSSWVxRatio'));
Te = (get_data(type1,'maxTeSWRatio'));
ThetaV = (get_data(type1,'SSSWDeflection'));
Msh = sum((N>= 1.5).*(Vx <= 0.5).*(Te >= 1.5));
Sw = sum((N>=0.50 & N <= 1.5).*(Vx >= 0.5 & Vx <= 1.5).*(ThetaV <= 15));
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions
%Average sownstaream and upstream for each type
function [avg] = average_downup(downstream,upstream)
    avg = mean([downstream;upstream]);
end
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
    
    %     if nargin <= 4
    %         data2=[];
    %         label2=[];
    %         label1=[];
    %     end
%     if isempty(data2)
%         label1=[];
%     end
    
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
        offsetXlabels=0;
    else
        binRange=bin_size;
        custom_bins = 1;
        bin_size = binRange(3)-binRange(2);
        if binRange(1) == -(binRange(3)-binRange(2))/2;
            offsetXlabels = 1;
        else
            offsetXlabels = 0;
        end
    end
    
    %Normalization Initialization if MMS data in given
    if nargin == 7
        %         C = histcounts(data3,[-Inf binRange max([data1;data2])]);
        %Normalize to Hourly, MMS Survey is 4.5-s, 800 datapoints per hour
        
        
        
        C = histcounts(data3,[-Inf binRange Inf]);
        if sum(contains(parameter,'Solar Wind Speed')) > 0 ||...
                sum(contains(parameter,'Solar Wind Density')) > 0 ||...
                sum(contains(parameter,'Magnetic Field Strength')) > 0 ||...
                sum(contains(parameter,'Mach Number')) > 0 ||...
                sum(contains(parameter,'Cone Angle')) > 0
            C = C./800;
            
        end
        
        %         if sum(contains(parameter,'ClosestDistance')) > 0
        %             %2 minutes per data point, per hour is 30 times.
        %             C = C./30;
        %         end
        %         C = C./sum(C);
        Cerror = C.^(-1/2);
        Cerror(Cerror == inf) = NaN;
        Cerror(Cerror == 0) = NaN;
        if sum(contains(parameter,'MMS')) > 0
            yAxisLabel = '# of Events normalized to MMS observations';
            if sum(contains(parameter,'Shear')) > 0 ||...
                    sum(contains(parameter,'E-TD')) > 0 ||...
                    sum(contains(parameter,'TD-BS')) > 0 ||...
                    sum(contains(parameter,'Shock')) > 0 ||...
                    sum(contains(parameter,'Transit')) > 0 ||...
                    sum(contains(parameter,'Closest')) > 0
                yAxisLabel = '# of Events normalized to MMS solar wind discontinuities';
            else
                yAxisLabel = '# of Events per hour';
                %                 Ctotal = sum(C);
                %                 HoursOfC = Ctotal /(60*60);
                %                 C = C./HoursOfC;
                
            end
        else
            yAxisLabel = '# of Events normalized to WIND observations';
        end
    else
        C=1;
        Cerror = 1;
        yAxisLabel = '# of Events / Total # Events';
    end
    
    
    %Plotting W
    if nargin ~= 4 %&& nargin ~= 3%Wewighting from Spacecraft Data Over All Time Range
        A = histcounts(data1,[-Inf binRange Inf]);
        B = histcounts(data2,[-Inf binRange Inf]);
        
        A_normalized = A./C;
        B_normalized =B./C;
        
        if length(C) == 1
            A_normalized = A./sum(A);
            B_normalized =B./sum(B);
            Aerror = (A).^(-1/2)./sum(A);
            Berror = (B).^(-1/2)./sum(B);
            
            Aerror = (A_normalized).^(-1).*((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
            Berror = (B_normalized).^(-1).*((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);
            
            Aerror = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
            Berror = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);
            
            %             Aerror = Aerror./ (A./sum(A));
            %             Berror = Berror./ (B./sum(B));
            %             Aerror = A_normalized.*(A.^(-1) + sum(A).^(-1) ).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) + sum(B).^(-1) ).^(1/2);
        else
            
            Aerror = A_normalized.*(A.^(-1) +      C.^(-1)).^(1/2);
            Berror = B_normalized.*(B.^(-1) +      C.^(-1)).^(1/2);
            
            %Binomial Error Bars
            dA = (((A_normalized).^(-1).*(1-A_normalized))./sum(A)).^(1/2); %Real Error
            dB = (((B_normalized).^(-1).*(1-B_normalized))./sum(B)).^(1/2); %Real Error
            
            dA = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2); %Real Error
            dB = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2); %Real Error
            
            %             dA = dA ./ (A./sum(A));
            %             dB = dB ./ (B./sum(B));
            
            C_normalized = C./sum(C);
            dC = ((C_normalized.*(1-C_normalized))./sum(C)).^(1/2);
            
            Aerror = A_normalized.*(dA + dC);
            Berror = B_normalized.*(dB + dC);
            
            Aerror = (A./C).*(dA + dC);
            Berror = (B./C).*(dB + dC);
            
            Aerror = (sum(A)./C).*(dA + dC);
            Berror = (sum(B)./C).*(dB + dC);
            
            Aerror = (sum(A)./C).*dA + (sum(C).*A./C./C).*dC;
            Berror = (sum(B)./C).*dB + (sum(C).*B./C./C).*dC;
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
        
        if length(C)==1
            
            bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
            if isempty(data2)
                
            else
                [XB,YB] = stairs(xlabels,B_normalized');
                XB = [XB;XB(end)+bin_size];
                YB = [YB;YB(end)];
                
                plot(XB,YB,'linewidth',3.5);
            end
            
            
        else
            if isempty(data2)
                
            else
                yyaxis left
            end
            bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
            
            if isempty(data2)
                
            else
                [XB,YB] = stairs(xlabels,B_normalized');
                XB = [XB;XB(end)+bin_size];
                YB = [YB;YB(end)];
                
                yyaxis right
                plot(XB,YB,'linewidth',3.5); % plot(XB,55/45.*YB,'linewidth',3.5,'color','r');
                set(gca,'ycolor','r')
            end
            
        end
        
        %         stairs(xlabels,B_normalized','linewidth',3.5);
        
        if A(end) ~= 0 || B(end) ~= 0
            xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
            xlim([binRange(1) binRange(end)+(binRange(3)-binRange(2))])
        else
            xticks(binRange)
            xlim([binRange(1) binRange(end)])
        end
        
        %Labels
        if isempty(data2)
            if custom_bins == 1 && xlabels(1) ~= 0 && (max(data1) >= binRange(end)  )
                xlabelCell = xticklabels;
                xlabelCell(1) = strcat('\leq',xlabelCell(1));
                xlabelCell(end) = strcat('\geq',xlabelCell(end));
                xticklabels(xlabelCell);
            elseif custom_bins == 1 && (max(data1) >= binRange(end)  )
                xlabelCell = xticklabels;
                xlabelCell(end) = strcat('\geq',xlabelCell(end));
                xticklabels(xlabelCell);
            end
        else
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
        end
        if offsetXlabels == 1
            xlabelCell{1} = '0';
            xticklabels(xlabelCell);
        else
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
        
        if isempty(data2)
            legend off
        else
        end
        hold on
        
        if length(C) == 1
            errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
            
            errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        else
            %Find the type with the greater number of events
            maxAerror = max(A_normalized + Aerror);
            maxBerror = max(B_normalized + Berror);
            [~,maxTypeNormalized] = max([maxAerror,maxBerror]);
            
            if maxTypeNormalized == 1
                
                d2Overd1Ratio = length(data2)/length(data1);
                if isempty(data2)
                else
                    yyaxis left
                end
                errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                d1_ylim = get(gca,'ylim');
                if isempty(data2)
                    
                else
                    yyaxis right
                    errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                    ylim([0 d2Overd1Ratio*d1_ylim(2)])
                end
            else
                d1Overd2Ratio = length(data2)/length(data1);
                yyaxis right
                errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                d2_ylim = get(gca,'ylim');
                if isempty(data2)
                    
                else
                    yyaxis left
                    errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                    ylim([0 d1Overd2Ratio*d2_ylim(2)])
                end
                
            end
            
            %             %Final Check
            %             yyaxis left
            %             d1_ylim = get(gca,'ylim');
            %             if d1_ylim(2) < max(A_normalized + Aerror)
            %                 d2Overd1Ratio = length(data2)/length(data1);
            %                 ylim([0 d2Overd1Ratio*d1_ylim(2)])
            %             end
            
            
            
            
        end
        
        
        
        hold off
    else
        %Single Type of Event Type
        A = histcounts(data1,[-Inf binRange Inf]);
        A_normalized = A./sum(A);
        %Aerror = A.^(-1/2)./sum(A);
        Aerror = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
        
        
        
        
        if A(1) == 0
            A = A(2:end);
            A_normalized = A_normalized(2:end);
            Aerror = Aerror(2:end);
        else
            binRange = [binRange(1) - bin_size, binRange];
        end
        
        xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        
        bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
        
        if A(end) ~= 0
            xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
            xlim([binRange(1) binRange(end)+(binRange(3)-binRange(2))])
        else
            xticks(binRange)
            xlim([binRange(1) binRange(end)])
        end
        
        %Labels
        if custom_bins == 1 && xlabels(1) ~= 0 && (max(data1) >= binRange(end) )
            xlabelCell = xticklabels;
            xlabelCell(1) = strcat('\leq',xlabelCell(1));
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        elseif custom_bins == 1 && (max(data1) >= binRange(end) )
            xlabelCell = xticklabels;
            xlabelCell(end) = strcat('\geq',xlabelCell(end));
            xticklabels(xlabelCell);
        end
        
        
        
        errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','black','LineWidth',1.5,'HandleVisibility','off')
        
        
        
        
        
        
        
        %xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        %bar(xlabels+bin_size/2,A_normalized',1,'edgecolor','none'); hold on
        %errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','black','LineWidth',1.5,'HandleVisibility','off')
        
        
        
        
        %Legend
        %Location of Max
        [~,I] = max(A(A~=inf));
        if I > length(binRange)/2
            if isempty(label1)
            else
                legend(label1,'FontSize',9,'Location','Northwest')
            end
        else
            if isempty(label1)
            else
                legend(label1,'FontSize',9,'Location','Northeast')
            end
        end
        hold on
        hold off
        
    end
    
    
    
    % Plot Properties
    if contains(parameter(1),'E-TD Magnitude') & length(C) ~=1
        %%%%%%set(gca,'YScale','log')
    elseif contains(parameter(1),'Solar Wind Speed') & length(C) ~=1
        %         yyaxis left
        %         set(gca,'YScale','log')
        %         yyaxis right
        %         set(gca,'YScale','log')
        %         %     elseif contains(parameter(1),'Shear Angle') & length(C) ~=1
        %         %         set(gca,'YScale','log')
        %
    end
    if length(C) == 1
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
            %         xlim([binRange(1) binRange(end)])
        end
    else
        if isempty(data2)
        else
            yyaxis left
        end
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
        end
        if isempty(data2)
        else
            yyaxis right
        end
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
        end
    end
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1);
    colormap('winter');
    if sum(contains(parameter(2),'[Re]')) > 0
        parameter(2) = {'[R_E]'};
    elseif sum(contains(parameter(2),'[Re^2]')) > 0
        parameter(2) = {'[{R_E}^2]'};
    end
    
    
    if sum(contains(parameter,'Shock Geometry')) > 0
        xlabel({'\Theta_{Bn}'})
    else
        xlabel(erase(erase(erase(erase(parameter,'Normalized'), 'MMS'),'Upstream'),'Downstream'),'FontSize',11)
    end
    title(strcat('f vs.', {' '}, parameter(1)),'FontSize',12,'FontWeight', 'normal')
    if length(C)==1
        ylabel(yAxisLabel,'FontSize',11)
        
        set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    else
        if isempty(data2)
        else
        yyaxis left
        end
        ylabel(yAxisLabel,'FontSize',11)
        set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
        if isempty(data2)
        else
        yyaxis right
        end
        ylabel(yAxisLabel,'FontSize',11)
        set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    end
    
    
    print(gcf,'-dpng','-r300', '-loose', strcat(fileName));
    savefig(strcat(fileName));
    
end
function [] = plot_singleHistogram(parameter,bin_size,data1,data3)
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
        offsetXlabels=0;
    else
        binRange=bin_size;
        custom_bins = 1;
        bin_size = binRange(3)-binRange(2);
        if binRange(1) == -(binRange(3)-binRange(2))/2;
            offsetXlabels = 1;
        else
            offsetXlabels = 0;
        end
    end
    
    %Normalization Initialization if MMS data in given
    if nargin == 4
        %         C = histcounts(data3,[-Inf binRange max([data1;data2])]);
        %Normalize to Hourly, MMS Survey is 4.5-s, 800 datapoints per hour
        
        
        
        C = histcounts(data3,[-Inf binRange Inf]);
        if sum(contains(parameter,'Solar Wind Speed')) > 0 ||...
                sum(contains(parameter,'Solar Wind Density')) > 0 ||...
                sum(contains(parameter,'Magnetic Field Strength')) > 0 ||...
                sum(contains(parameter,'Mach Number')) > 0 ||...
                sum(contains(parameter,'Cone Angle')) > 0
            C = C./800;
            
        end
        
        %         if sum(contains(parameter,'ClosestDistance')) > 0
        %             %2 minutes per data point, per hour is 30 times.
        %             C = C./30;
        %         end
        %         C = C./sum(C);
        Cerror = C.^(-1/2);
        Cerror(Cerror == inf) = NaN;
        Cerror(Cerror == 0) = NaN;
        if sum(contains(parameter,'MMS')) > 0
            yAxisLabel = '# of Events normalized to MMS observations';
            if sum(contains(parameter,'Shear')) > 0 ||...
                    sum(contains(parameter,'E-TD')) > 0 ||...
                    sum(contains(parameter,'TD-BS')) > 0 ||...
                    sum(contains(parameter,'Shock')) > 0 ||...
                    sum(contains(parameter,'Transit')) > 0 ||...
                    sum(contains(parameter,'Closest')) > 0
                yAxisLabel = '# of Events normalized to MMS solar wind discontinuities';
            else
                yAxisLabel = '# of Events per hour';
                %                 Ctotal = sum(C);
                %                 HoursOfC = Ctotal /(60*60);
                %                 C = C./HoursOfC;
                
            end
        else
            yAxisLabel = '# of Events normalized to WIND observations';
        end
    else
        C=1;
        Cerror = 1;
        yAxisLabel = '# of Events / Total # Events';
    end
    
    
    %Plotting W
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
            
            Aerror = (A_normalized).^(-1).*((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
            Berror = (B_normalized).^(-1).*((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);
            
            Aerror = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
            Berror = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);
            
            %             Aerror = Aerror./ (A./sum(A));
            %             Berror = Berror./ (B./sum(B));
            %             Aerror = A_normalized.*(A.^(-1) + sum(A).^(-1) ).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) + sum(B).^(-1) ).^(1/2);
        else
            
            Aerror = A_normalized.*(A.^(-1) +      C.^(-1)).^(1/2);
            Berror = B_normalized.*(B.^(-1) +      C.^(-1)).^(1/2);
            
            %Binomial Error Bars
            dA = (((A_normalized).^(-1).*(1-A_normalized))./sum(A)).^(1/2); %Real Error
            dB = (((B_normalized).^(-1).*(1-B_normalized))./sum(B)).^(1/2); %Real Error
            
            dA = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2); %Real Error
            dB = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2); %Real Error
            
            %             dA = dA ./ (A./sum(A));
            %             dB = dB ./ (B./sum(B));
            
            C_normalized = C./sum(C);
            dC = ((C_normalized.*(1-C_normalized))./sum(C)).^(1/2);
            
            Aerror = A_normalized.*(dA + dC);
            Berror = B_normalized.*(dB + dC);
            
            Aerror = (A./C).*(dA + dC);
            Berror = (B./C).*(dB + dC);
            
            Aerror = (sum(A)./C).*(dA + dC);
            Berror = (sum(B)./C).*(dB + dC);
            
            Aerror = (sum(A)./C).*dA + (sum(C).*A./C./C).*dC;
            Berror = (sum(B)./C).*dB + (sum(C).*B./C./C).*dC;
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
        
        if length(C)==1
            
            bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
            
            [XB,YB] = stairs(xlabels,B_normalized');
            XB = [XB;XB(end)+bin_size];
            YB = [YB;YB(end)];
            
            plot(XB,YB,'linewidth',3.5);
            
            
            
        else
            
            yyaxis left
            bar(xlabels+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
            
            
            [XB,YB] = stairs(xlabels,B_normalized');
            XB = [XB;XB(end)+bin_size];
            YB = [YB;YB(end)];
            
            yyaxis right
            plot(XB,YB,'linewidth',3.5); % plot(XB,55/45.*YB,'linewidth',3.5,'color','r');
            set(gca,'ycolor','r')
            
        end
        
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
        if offsetXlabels == 1
            xlabelCell{1} = '0';
            xticklabels(xlabelCell);
        else
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
        
        if length(C) == 1
            errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
            
            errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        else
            %Find the type with the greater number of events
            maxAerror = max(A_normalized + Aerror);
            maxBerror = max(B_normalized + Berror);
            [~,maxTypeNormalized] = max([maxAerror,maxBerror]);
            
            if maxTypeNormalized == 1
                
                d2Overd1Ratio = length(data2)/length(data1);
                yyaxis left
                errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                d1_ylim = get(gca,'ylim');
                yyaxis right
                errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                ylim([0 d2Overd1Ratio*d1_ylim(2)])
            else
                d1Overd2Ratio = length(data2)/length(data1);
                yyaxis right
                errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                d2_ylim = get(gca,'ylim');
                yyaxis left
                errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                ylim([0 d1Overd2Ratio*d2_ylim(2)])
                
                
            end
            
        end
        
        hold off
        
    end
    
    
    
    % Plot Properties
    if contains(parameter(1),'E-TD Magnitude') & length(C) ~=1
        %%%%%%set(gca,'YScale','log')
    elseif contains(parameter(1),'Solar Wind Speed') & length(C) ~=1
        %         yyaxis left
        %         set(gca,'YScale','log')
        %         yyaxis right
        %         set(gca,'YScale','log')
        %         %     elseif contains(parameter(1),'Shear Angle') & length(C) ~=1
        %         %         set(gca,'YScale','log')
        %
    end
    if length(C) == 1
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
            %         xlim([binRange(1) binRange(end)])
        end
    else
        yyaxis left
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
        end
        yyaxis right
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
        end
    end
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1);
    colormap('winter');
    if sum(contains(parameter(2),'[Re]')) > 0
        parameter(2) = {'[R_E]'};
    elseif sum(contains(parameter(2),'[Re^2]')) > 0
        parameter(2) = {'[{R_E}^2]'};
    end
    
    
    if sum(contains(parameter,'Shock Geometry')) > 0
        xlabel({'\Theta_{Bn}'})
    else
        xlabel(erase(erase(erase(erase(parameter,'Normalized'), 'MMS'),'Upstream'),'Downstream'),'FontSize',11)
    end
    title(strcat('f vs.', {' '}, parameter(1)),'FontSize',12,'FontWeight', 'normal')
    if length(C)==1
        ylabel(yAxisLabel,'FontSize',11)
        
        set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    else
        yyaxis left
        ylabel(yAxisLabel,'FontSize',11)
        set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
        yyaxis right
        ylabel(yAxisLabel,'FontSize',11)
        set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    end
    
    
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
        scatter(data1a,data1b,85,'filled','Displayname',label1)
        fileName = strcat(cell2mat(parameter(1)),'-',cell2mat(parameter(2)),'_',label1);
    elseif nargin == 7
        scatter(data1a,data1b,85,'filled','Displayname',label1)
        hold on
        scatter(data2a,data2b,85,'filled','Displayname',label2)
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
function [data] = get_data(Event_Type,Event_Parameter,filter)
    Page_number = 1;
    Column_number = 'A1';
    %     Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabaseFebruary26.xlsx';
    Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabase_Typical_Final.xlsx';
    
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
        case 'all'
            Page_number = 1;
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
            
        case 'Boundary1Vmag'
            Column_number = 'FW';
            
        case 'Boundary2Vmag'
            Column_number = 'FX';
            
        case 'BoundaryExpansionC'
            Column_number = 'FY';
            
        case 'BoundaryExpansionCminus'
            Column_number = 'FZ';
            
        case 'SSCoreNRatio'
            Column_number = 'GA';
            
        case 'SSCoreDeflection'
            Column_number = 'GB';
            
        case 'SSCoreVxRatio'
            Column_number = 'GC';
            
        case 'SSCoreDPRatio'
            Column_number = 'GD';
            
        case 'SSCoreTRatio'
            Column_number = 'GE';
            
        case 'SSNVCorr'
            Column_number = 'GF';
            
        case 'SSSWNRatio'
            Column_number = 'GG';
            
        case 'SSSWDeflection'
            Column_number = 'GH';
            
        case 'SSSWVxRatio'
            Column_number = 'GI';
            
        case 'SSSWDPRatio'
            Column_number = 'GJ';
            
        case 'SSSWTRatio'
            Column_number = 'GK';
            
        case 'SSCoreBRatio'
            Column_number = 'GL';
            
        case 'SSSWBRatio'
            Column_number = 'GM';
            
        case 'SSNBCorr'
            Column_number = 'GN';
            
            
            
        case 'SS Size'
            Column_number = 'GO';
            
        case 'SSSizeinElectronInertialLengths'
            Column_number = 'GP';
            
            
        case 'SS Speed'
            Column_number = 'GQ';
            
        case 'SS Duration'
            Column_number = 'GR';
            
            
        case 'SS Core Duration Ratio'
            Column_number = 'GS';
            
            
            
            
        case 'CoreSize1IonLength'
            Column_number = 'GT';
            
        case 'CoreSize2IonLength'
            Column_number = 'GU';
            
        case 'SSStartFraction'
            Column_number = 'GV';
            
        case 'CoreSSBeginIonLength1'
            Column_number = 'GW';
            
        case 'CoreSSBeginIonLength2'
            Column_number = 'GX';
            
        case 'Age2'
            Column_number = 'GY';
            
        case 'SSSizeBulkV'
            Column_number = 'GZ';
            
        case 'SSSizeBulkVinIonGyro'
            Column_number = 'HA';
            
        case 'SSSizewithDurationRatioandSize'
            Column_number = 'HB';
            
        case 'Age2p'
            Column_number = 'HC';
            
        case 'downstream N'
            Column_number = 'HD';
            
        case 'upstream N'
            Column_number = 'HE';
            
        case 'CoreSWVelocityDeflection'
            Column_number = 'HF';
            
        case 'maxTecoreratio'
            Column_number = 'HG';
            
        case 'maxTeSWRatio'
            Column_number = 'HH';
            
        case 'SSenTRatio'
            Column_number = 'HI';
            
        case 'SSinTRatio'
            Column_number = 'HJ';
            
        case 'SSSizeinIonInertialLengths'
            Column_number = 'HK';
            
        case 'MaxNSigma'
            Column_number = 'HL';
            
        case 'ionGyroradius'
            Column_number = 'HM';
            
        case 'ionInertiallength'
            Column_number = 'HN';
            
        case 'SSSizeBulkVinElectronGyro'
            Column_number = 'HO';
            
        case 'electronGyroradius'
            Column_number = 'HP';
            
        case 'electronInertiallength'
            Column_number = 'HQ';
            
    end
    
    Column_range = strcat(Column_number,'2:',Column_number,'999');
    
    [data,txt] = xlsread(Database_Directory,Page_number,Column_range);
    data = data(~isnan(data));
    if strcmp(Event_Parameter,'Event Date')
        data = txt(txt~="");
    end
    
    if strcmp(Event_Parameter,'ClosestDistance')
        data = data(data ~= -1);
    end
    %A =datetime(type1_date,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS')
    
    if nargin == 3
    data = data(filter);
    end
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
