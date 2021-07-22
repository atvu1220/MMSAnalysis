close all
clear;
tic


cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Histograms'

%Strings for types Of events
type1 = 'SS';
name1 = 'Substructure';
type2 = 'NS';
name2 = 'No Substructure';

type1b = 'HFA SS';
type2b = 'HFA NS';

%2021 Edit First revision for Paper. No SHFAs are included.
type1 = 'FB SS';
name1 = 'FB SS';
type2 = 'FB NS';
name2 = 'FB NS';
type1b = type1;
name1b = name1;
type2b = type2;
name2b = name2;

FBnumber = size(get_data('FB','Shear Angle'))
HFAnumber = size(get_data('HFA','Shear Angle'))
SHFAnumber = size(get_data('SHFA','Shear Angle'))

SSnumber = size(get_data('SS','Shear Angle'))
NSnumber = size(get_data('NS','Shear Angle'))

FBSSnumber = size(get_data('FB SS','Shear Angle'))
FBNSnumber = size(get_data('FB NS','Shear Angle'))

HFASSnumber = size(get_data('HFA SS','Shear Angle'))
HFANSnumber = size(get_data('HFA NS','Shear Angle'))

SHFASSnumber = size(get_data('SHFA SS','Shear Angle'))
SHFANSnumber = size(get_data('SHFA NS','Shear Angle'))


%Load WIND Parameters
% load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/WIND_2017_2019.mat')
% load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/WIND_2017_2019_TD4.mat')
% load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/MMS_2017_2019_TD8.mat')
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/MMS_2017_2019_TD12.mat')
load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/OMNI_atMMStimes_2017_2019b.mat')
%load('/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis/OMNI_2017_2019.mat')
%load('/Users/andrewvu/iCloud Drive (Archive) - 1/Research/MATLAB Analysis/MMS_2017_2019_TD10.mat')
A=5;
%% Plots
% %% Substructures (SS) & No Substructures (NS)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Get the solar wind velocity
%Get the bow shock normal
%Get the velocity component along BS normal
%Get the density
%Get the magnetic field
% 
%             
% DownstreamVx_type1 = get_data(type1,'Vdownx');
% DownstreamVy_type1 = get_data(type1,'Vdowny');
% DownstreamVz_type1 = get_data(type1,'Vdownz');
% 
% UpstreamVx_type1 = get_data(type1,'Vupx');
% UpstreamVy_type1 = get_data(type1,'Vupy');
% UpstreamVz_type1 = get_data(type1,'Vupz');
% 
% DownstreamVx_type2 = get_data(type2,'Vdownx');
% DownstreamVy_type2 = get_data(type2,'Vdowny');
% DownstreamVz_type2 = get_data(type2,'Vdownz');
% 
% UpstreamVx_type2 = get_data(type2,'Vupx');
% UpstreamVy_type2 = get_data(type2,'Vupy');
% UpstreamVz_type2 = get_data(type2,'Vupz');
% 
% 
% 
% ClosestBSNx_type1 = get_data(type1,'ClosestBSNx');
% ClosestBSNy_type1 = get_data(type1,'ClosestBSNy');
% ClosestBSNz_type1 = get_data(type1,'ClosestBSNz');
% 
% ClosestBSNx_type2 = get_data(type2,'ClosestBSNx');
% ClosestBSNy_type2 = get_data(type2,'ClosestBSNy');
% ClosestBSNz_type2 = get_data(type2,'ClosestBSNz');
% 
% 
%        
% DownstreamVA_type1 = get_data(type1,'Alfven V Before');
% UpstreamVA_type1 = get_data(type1,'Alfven V After');
% 
% DownstreamVA_type2 = get_data(type2,'Alfven V Before');
% UpstreamVA_type2 = get_data(type2,'Alfven V After');
% 
% 
% 
% 
% downstreamNormalSpeed_type1 = DownstreamVx_type1.*ClosestBSNx_type1 + DownstreamVy_type1.*ClosestBSNy_type1 + DownstreamVz_type1.*ClosestBSNz_type1;
% upstreamNormalSpeed_type1   = UpstreamVx_type1.*  ClosestBSNx_type1 + UpstreamVy_type1.*  ClosestBSNy_type1 + UpstreamVz_type1.*  ClosestBSNz_type1;
% 
% downstreamNormalSpeed_type2 = DownstreamVx_type2.*ClosestBSNx_type2 + DownstreamVy_type2.*ClosestBSNy_type2 + DownstreamVz_type2.*ClosestBSNz_type2;
% upstreamNormalSpeed_type2   = UpstreamVx_type2.*  ClosestBSNx_type2 + UpstreamVy_type2.*  ClosestBSNy_type2 + UpstreamVz_type2.*  ClosestBSNz_type2;
% 
% downstreamMA_type1 = abs(downstreamNormalSpeed_type1 ./ DownstreamVA_type1);
% upstreamMA_type1 = abs(upstreamNormalSpeed_type1 ./ UpstreamVA_type1);
% 
% downstreamMA_type2 = abs(downstreamNormalSpeed_type2 ./ DownstreamVA_type2);
% upstreamMA_type2 = abs(upstreamNormalSpeed_type2 ./ UpstreamVA_type2);
% 
% 
% 
% 
% avgNormalSpeed_type1 = mean([DownstreamVx_type1,UpstreamVx_type1],2).*ClosestBSNx_type1 + mean([DownstreamVy_type1,UpstreamVy_type1],2).*ClosestBSNy_type1 + mean([DownstreamVz_type1,UpstreamVz_type1],2).*ClosestBSNz_type1;
% % upstreamNormalSpeed_type1   = UpstreamVx_type1.*  ClosestBSNx_type1 + UpstreamVy_type1.*  ClosestBSNy_type1 + UpstreamVz_type1.*  ClosestBSNz_type1;
% 
% avgNormalSpeed_type2 = mean([DownstreamVx_type2,UpstreamVx_type2],2).*ClosestBSNx_type2 + mean([DownstreamVy_type2,UpstreamVy_type2],2).*ClosestBSNy_type2 + mean([DownstreamVz_type2,UpstreamVz_type2],2).*ClosestBSNz_type2;
% 
% avgNormalSpeed_type1 = max([DownstreamVx_type1,UpstreamVx_type1],[],2).*ClosestBSNx_type1 + min([DownstreamVy_type1,UpstreamVy_type1],[],2).*ClosestBSNy_type1 + min([DownstreamVz_type1,UpstreamVz_type1],[],2).*ClosestBSNz_type1
% avgNormalSpeed_type2 = max([DownstreamVx_type2,UpstreamVx_type2],[],2).*ClosestBSNx_type2 + min([DownstreamVy_type2,UpstreamVy_type2],[],2).*ClosestBSNy_type2 + min([DownstreamVz_type2,UpstreamVz_type2],[],2).*ClosestBSNz_type2
% 
% 
% 
% avgMA_type1 = abs(avgNormalSpeed_type1 ./ DownstreamVA_type1);
% avgMA_type2 = abs(avgNormalSpeed_type2 ./ DownstreamVA_type2);
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},6:2:14,...
%     avgMA_type1,name1,...
%     avgMA_type2,name2,...
%     MMS_MachNumber)
% 
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},8:1.5:14,...
%     combine_downup(downstreamMA_type1,upstreamMA_type1),name1,...
%     combine_downup(downstreamMA_type2,upstreamMA_type2),name2,...
%     MMS_MachNumber)
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},8:2:14,...
%     combine_downup(downstreamMA_type1,upstreamMA_type1),name1,...
%     combine_downup(downstreamMA_type2,upstreamMA_type2),name2,...
%     MMS_MachNumber)
% 
% 
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},8:4:16,...
%     combine_downup(downstreamMA_type1,upstreamMA_type1),name1,...
%     combine_downup(downstreamMA_type2,upstreamMA_type2),name2,...
%     MMS_MachNumber)
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},6:2:13,...
%     downstreamMA_type1,name1,...
%     downstreamMA_type2,name2,...
%     MMS_MachNumber)
% 
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},6:2:13,...
%     upstreamMA_type1,name1,...
%     upstreamMA_type2,name2,...
%     MMS_MachNumber)
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},6:2:13,...
%     mean([downstreamMA_type1,upstreamMA_type1],2),name1,...
%     mean([downstreamMA_type2,upstreamMA_type2],2),name2,...
%     MMS_MachNumber)


%             
% plot_histogram({'MMS Normalized Alfvén Mach Number Downstream of BS Normal','[M_A]'},7:3:16,...
%     rmoutliers(get_data(type1,'Downstream BS MA'),'thresholdfactor',3),name1,...
%     rmoutliers(get_data(type2,'Downstream BS MA'),'thresholdfactor',3),name2,...
%     MMS_MachNumber)
% 
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number Upstream of BS Normal','[M_A]'},7:3:16,...
%     rmoutliers(get_data(type1,'Upstream BS MA'),'thresholdfactor',3),name1,...
%     rmoutliers(get_data(type2,'Upstream BS MA'),'thresholdfactor',3),name2,...
%     MMS_MachNumber)
% 
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},7:3:16,...
%     [rmoutliers(get_data(type1,'Downstream BS MA'),'thresholdfactor',3);rmoutliers(get_data(type1,'Upstream BS MA'),'thresholdfactor',3)],name1,...
%     [rmoutliers(get_data(type2,'Downstream BS MA'),'thresholdfactor',3);rmoutliers(get_data(type2,'Upstream BS MA'),'thresholdfactor',3)],name2,...
%     MMS_MachNumber)    
% 
% 
% plot_histogram({'OMNI Normalized Alfvén Mach Number at BS Normal','[M_A]'},7:3:16,...
%     get_data(type1,'OMNI MA'),name1,...
%     get_data(type2,'OMNI MA'),name2,...
%     Omni_MaBS_MMStimes)   
% 
% % plot_histogram({'MMS Normalized Alfvén Mach Number at BS Normal','[M_A]'},6:2:14,...
% %     mean([get_data(type1,'Downstream BS MA'),get_data(type1,'Upstream BS MA')],2),name1,...
% %     mean([get_data(type2,'Downstream BS MA'),get_data(type2,'Upstream BS MA')],2),name2,...
% %     MMS_MachNumber)    
% 
% 
% 
% plot_histogram({'MMS Normalized Msg Mach Number','[M_A]'},5:0.5:7.5,...
%     get_data(type1,'OMNI MGSMA'),name1,...
%     get_data(type2,'OMNI MGSMA'),name2,...
%     Omni_MGSM)
% 
% plot_histogram({'MMS Normalized Msg Mach Number Downstream of BS Normal','[M_A]'},5:1:7,...
%     rmoutliers(get_data(type1,'Downstream BS MGSMA'),'thresholdfactor',3),name1,...
%     rmoutliers(get_data(type2,'Downstream BS MGSMA'),'thresholdfactor',3),name2,...
%     Omni_MGSM)
% 
% plot_histogram({'MMS Normalized Msg Mach Number Upstream of BS Normal','[M_A]'},5:1:7,...
%     rmoutliers(get_data(type1,'Upstream BS MGSMA'),'thresholdfactor',3),name1,...
%     rmoutliers(get_data(type2,'Upstream BS MGSMA'),'thresholdfactor',3),name2,...
%     Omni_MGSM)
% 
% plot_histogram({'MMS Normalized Msg Mach Number at BS Normal','[M_A]'},5:1:7,...
%     [rmoutliers(get_data(type1,'Downstream BS MGSMA'),'thresholdfactor',3);rmoutliers(get_data(type1,'Upstream BS MGSMA'),'thresholdfactor',3)],name1,...
%     [rmoutliers(get_data(type2,'Downstream BS MGSMA'),'thresholdfactor',3);rmoutliers(get_data(type2,'Upstream BS MGSMA'),'thresholdfactor',3)],name2,...
%     Omni_MGSM)  
% 
% % plot_histogram({'MMS Normalized Magnetosonic Mach Number at BS Normal','[M_A]'},2:1:8,...
% %     mean([get_data(type1,'Downstream BS MGSMA'),get_data(type1,'Upstream BS MGSMA')],2),name1,...
% %     mean([get_data(type2,'Downstream BS MGSMA'),get_data(type2,'Upstream BS MGSMA')],2),name2)    



plot_histogram({'OMNI Magnetic Field Strength','[nT]'},2:2:12,...
    Omni_B_MMStimes,'Omni',...
    MMS_BMag,'MMS')

plot_histogram({'OMNI Magnetic Field Strength','[nT]'},2:2:12,...
    get_data(type1,'Omni B'),name1,...
    get_data(type2,'Omni B'),name2)

plot_histogram({'OMNI Normalized Magnetic Field Strength','[nT]'},2:2:12,...
    get_data(type1,'Omni B'),name1,...
    get_data(type2,'Omni B'),name2,...
    Omni_B_MMStimes)


plot_histogram({'OMNI Ion Bulk Speed','km/s'},300:100:600,...
    Omni_V_MMStimes,'Omni',...
    MMS_SWSpeed,'MMS')

plot_histogram({'OMNI Ion Bulk Speed','km/s'},200:100:600,...
    get_data(type1,'Omni Speed'),name1,...
    get_data(type2,'Omni Speed'),name2)

plot_histogram({'OMNI Normalized Ion Bulk Speed','km/s'},250:50:550,...
    get_data(type1,'Omni Speed'),name1,...
    get_data(type2,'Omni Speed'),name2,...
    Omni_V_MMStimes)

plot_histogram({'OMNI Ion Bulk Speed at BS','km/s'},250:100:550,...
    get_data(type1,'Omni Speed').* get_data(type1,'Omni Ma BS')./ get_data(type1,'Omni Ma'),name1,...
    get_data(type2,'Omni Speed').* get_data(type2,'Omni Ma BS')./ get_data(type2,'Omni Ma'),name2)

plot_histogram({'OMNI Normalized Ion Bulk Speed at BS','km/s'},250:100:550,...
    get_data(type1,'Omni Speed').* get_data(type1,'Omni Ma BS')./ get_data(type1,'Omni Ma'),name1,...
    get_data(type2,'Omni Speed').* get_data(type2,'Omni Ma BS')./ get_data(type2,'Omni Ma'),name2,...
    Omni_V_MMStimes.*Omni_BSNormalx_MMStimes)

plot_histogram({'OMNI Ion Density','cm^{-3}'},5:5:25,...
    Omni_n_MMStimes,'Omni',...
    MMS_Nion,'MMS')

plot_histogram({'OMNI Ion Density','cm^{-3}'},5:5:25,...
    get_data(type1,'Omni Density'),name1,...
    get_data(type2,'Omni Density'),name2)

plot_histogram({'OMNI Normalized Ion Density','cm^{-3}'},5:5:25,...
    get_data(type1,'Omni Density'),name1,...
    get_data(type2,'Omni Density'),name2,...
    Omni_n_MMStimes)


plot_histogram({'OMNI Ma Number','[M_A]'},5:3:16,...
    abs(Omni_Ma_MMStimes),'Omni',...
    abs(MMS_MachNumber),'MMS')

plot_histogram({'OMNI Ma Number','[M_A]'},5:3:16,...
    get_data(type1,'Omni Ma'),name1,...
    get_data(type2,'Omni Ma'),name2)

plot_histogram({'OMNI Normalized Ma Number','[M_A]'},5:3:16,...
    get_data(type1,'Omni Ma'),name1,...
    get_data(type2,'Omni Ma'),name2,...
    abs(Omni_Ma_MMStimes))


plot_histogram({'OMNI Ma Number at BS','[M_A]'},6:2:16,...
   abs(Omni_MaBS_MMStimes),'Omni',...
    abs(Omni_MaBS_MMStimes),'Omni')

plot_histogram({'OMNI Ma Number at BS','[M_A]'},6:3:16,...
    get_data(type1,'Omni Ma BS'),name1,...
    get_data(type2,'Omni Ma BS'),name2)

plot_histogram({'OMNI Normalized Ma Number at BS','[M_A]'},6:3:16,...
    get_data(type1,'Omni Ma BS'),name1,...
    get_data(type2,'Omni Ma BS'),name2,...
    abs(Omni_MaBS_MMStimes))


 plot_histogram({'OMNI Ms Number','[M_A]'},4:1:8,...
    Omni_Ms_MMStimes,'Omni',...
    Omni_Ms_MMStimes,'Omni')

 plot_histogram({'OMNI Ms Number','[M_A]'},4:1:8,...
    get_data(type1,'Omni Ms'),name1,...
    get_data(type2,'Omni Ms'),name2)

 plot_histogram({'OMNI Normalized Ms Number','[M_A]'},4:1:8,...
    get_data(type1,'Omni Ms'),name1,...
    get_data(type2,'Omni Ms'),name2,...
    Omni_Ms_MMStimes)


 plot_histogram({'OMNI Ms Number at BS','[M_A]'},0:1:7,...
    abs(Omni_MsBS_MMStimes),'Omni',...
    abs(Omni_MsBS_MMStimes),'Omni')

 plot_histogram({'OMNI Ms Number at BS','[M_A]'},0:1:7,...
    get_data(type1,'Omni Ms BS'),name1,...
    get_data(type2,'Omni Ms BS'),name2)

 plot_histogram({'OMNI Normalized Ms Number at BS','[M_A]'},3:1:7,...
    get_data(type1,'Omni Ms BS'),name1,...
    get_data(type2,'Omni Ms BS'),name2,...
    abs(Omni_MsBS_MMStimes))  


 plot_histogram({'OMNI Theta_{Bn}','[°]'},0:15:90,...
    angle90(Omni_thetaBn_MMStimes),'Omni',...
    combine_downup(angle90(MMS_Shock_Down_Angle),angle90(MMS_Shock_Up_Angle)),'MMS')

 plot_histogram({'OMNI Theta_{Bn}','[°]'},0:15:90,...
    angle90(get_data(type1,'Omni thetaBn')),name1,...
    angle90(get_data(type2,'Omni thetaBn')),name2)

 plot_histogram({'OMNI Normalized Theta_{Bn}','[°]'},0:20:90,...
    angle90(get_data(type1,'Omni thetaBn')),name1,...
    angle90(get_data(type2,'Omni thetaBn')),name2,...
    angle90(Omni_thetaBn_MMStimes))


 plot_histogram({'OMNI Cone Angle','[°]'},0:30:180,...
    Omni_ConeAngle_MMStimes,'Omni',...
    MMS_ConeAngle,'MMS') 

 plot_histogram({'OMNI Cone Angle','[°]'},0:30:180,...
    get_data(type1,'Omni coneAngle'),name1,...
    get_data(type2,'Omni coneAngle'),name2) 

 plot_histogram({'OMNI Normalized Cone Angle','[°]'},0:30:180,...
    get_data(type1,'Omni coneAngle'),name1,...
    get_data(type2,'Omni coneAngle'),name2,...
    Omni_ConeAngle_MMStimes) 
            
% HFA_M =  get_data(type1,'OMNI MGSMA').*get_data(type1,'BS_x');
% FB_M =  get_data(type2,'OMNI MGSMA').*get_data(type2,'BS_x');
% 
% plot_histogram({'MMS Normalized Msg Mach Number','[M_A]'},5:0.5:7.5,...
%     HFA_M,name1,...
%     FB_M,name2,...
%     Omni_MGSM)
% plot_histogram({'MMS Normalized Msg Mach Number','[M_A]'},5:0.5:7.5,...
%     get_data(type1,'OMNI MGSMA'),name1,...
%     get_data(type2,'OMNI MGSMA'),name2,...
%     Omni_MGSM)


% [~,Index] = min([angle90(get_data('SS','downstreamShockAngle')),angle90(get_data('SS','upstreamShockAngle'))],[],2);
% 
% thetaBN = [angle90(get_data('SS','downstreamShockAngle')),angle90(get_data('SS','upstreamShockAngle'))];
% 
% thetaBnmin = [];
% for i=1:length(thetaBN)
%     thetaBnmin(i) = min(thetaBN(i,:));
%     
% end
% scatter( get_data('SS','ClosestBSNx'), thetaBnmin)
% hold on;
% thetaBN = [angle90(get_data('NS','downstreamShockAngle')),angle90(get_data('NS','upstreamShockAngle'))];
% thetaBnmin = [];
% for i=1:length(thetaBN)
%     thetaBnmin(i) = min(thetaBN(i,:));
% end
% scatter( get_data('NS','ClosestBSNx'), thetaBnmin)



%Calculate alfven much number
%% Plots for Paper

plot_histogram({'MMS Closest Distance to Bow Shock','[Re]'},1:2:5,...
    MMS_TDClosestDistance,'MMS',...
    MMS_TDClosestDistance,'MMS')

plot_histogram({'MMS Closest Distance to Bow Shock','[Re]'},1:2:5,...
    rmoutliers(get_data(type1,'ClosestDistance'),'thresholdfactor',3),name1,...
    rmoutliers(get_data(type2,'ClosestDistance'),'thresholdfactor',3),name2)


plot_histogram({'MMS Normalized Closest Distance to Bow Shock','[Re]'},1:2:5,...
    rmoutliers(get_data(type1,'ClosestDistance'),'thresholdfactor',3),name1,...
    rmoutliers(get_data(type2,'ClosestDistance'),'thresholdfactor',3),name2,...
    MMS_TDClosestDistance)

% plot_histogram({'MMS Closest Distance to Bow Shock','[Re]'},1:2:5,...
%    MMS_TDClosestDistance,'MMS')
% 
% plot_histogram({'Closest Distance to Bow Shock','[Re]'},1:2:5,...
%     rmoutliers(get_data(type1,'ClosestDistance'),'thresholdfactor',3),name1,...
%     rmoutliers(get_data(type2,'ClosestDistance'),'thresholdfactor',3),name2)


plot_histogram({'MMS Normalized Distance to TD-BS Connection Point','[Re]'},1:2:5,...
    rmoutliers(get_data(type1b,'DistancetoCSBS'),'thresholdfactor',3),name1,...
    rmoutliers(get_data(type2b,'DistancetoCSBS'),'thresholdfactor',3),name2,...
    MMS_BSCSDistance(MMS_magneticShear > 30))

plot_histogram({'MMS Normalized TD-BS Angle','[\theta]'},0:15:90,...
    angle90(get_data(type1b,'CS-BS Angle')),name1,...
    angle90(get_data(type2b,'CS-BS Angle')),name2,...
    angle90(MMS_CSBS_Angle(MMS_magneticShear > 30)))

plot_histogram({'Inner Edges Angle','[\theta]'},0:30:180,...
    rmoutliers(get_data(type1,'Event Edges Angle'),'thresholdfactor',3),name1,...
    rmoutliers(get_data(type2,'Event Edges Angle'),'thresholdfactor',3),name2)

plot_histogram({'Width','[Re]'},0:1:3,...
    rmoutliers(get_data(type1,'Event Width'),'thresholdfactor',5),name1,...
    rmoutliers(get_data(type2,'Event Width'),'thresholdfactor',5),name2) %113

plot_histogram({'Height','[Re]'},0:2:8,...
    rmoutliers(get_data(type1,'Event Height'),'thresholdfactor',2),name1,...
    rmoutliers(get_data(type2,'Event Height'),'thresholdfactor',2),name2) %112

%Expansion Speed 2020
% plot_histogram({'Inner Edges Expansion','[km/s]'},-150:100:600,...
%     rmoutliers(get_data(type1,'BoundaryExpansion2020'),'thresholdfactor',10),name1,...
%     rmoutliers(get_data(type2,'BoundaryExpansion2020'),'thresholdfactor',10),name2) 

plot_histogram({'Inner Edges Expansion','[km/s]'},-150:100:600,...
    rmoutliers(get_data(type1,'BoundaryExpansion2020'),'thresholdfactor',8),name1,...
    rmoutliers(get_data(type2,'BoundaryExpansion2020'),'thresholdfactor',8),name2) 

% plot_histogram({'Area','[Re^2]'},0:0.5:4,...
%     rmoutliers(get_data(type1,'Event Area'),'thresholdfactor',10),name1,...
%     rmoutliers(get_data(type2,'Event Area'),'thresholdfactor',10),name2) %101


plot_histogram({'MMS Magnetic Shear Angle','[\theta]'},0:15:90,...
    MMS_magneticShear,'MMS',...
    MMS_magneticShear,'MMS')

plot_histogram({'MMS Magnetic Shear Angle','[\theta]'},0:15:90,...
    rmoutliers(get_data(type1,'Shear Angle'),'percentiles',[0 100]),name1,...
    rmoutliers(get_data(type2,'Shear Angle'),'percentiles',[0 100]),name2)


plot_histogram({'MMS Normalized Magnetic Shear Angle','[\theta]'},0:20:90,...
    rmoutliers(get_data(type1,'Shear Angle'),'percentiles',[0 100]),name1,...
    rmoutliers(get_data(type2,'Shear Angle'),'percentiles',[0 100]),name2,...
    MMS_magneticShear)
% 
% plot_histogram({'MMS Normalized Magnetic Shear Angle','[\theta]'},0:15:120,...
%     get_data(type1,'Shear Angle'),name1,...
%     get_data(type2,'Shear Angle'),name2)


% plot_histogram({'MMS Normalized Upstream Shock Geometry','[\theta]'},0:15:90,...
%     angle90(get_data(type1,'upstreamShockAngle')),name1,...
%     angle90(get_data(type2,'upstreamShockAngle')),name2,...
%     angle90(MMS_Shock_Up_Angle))
% 
% plot_histogram({'MMS Normalized Downstream Shock Geometry','[\theta]'},0:15:90,...
%     angle90(get_data(type1,'downstreamShockAngle')),name1,...
%     angle90(get_data(type2,'downstreamShockAngle')),name2,...
%     angle90(MMS_Shock_Down_Angle))


plot_histogram({'OMNI MMS ThetaBn','[\theta]'},0:15:90,...
    angle90(Omni_thetaBn_MMStimes),'Omni',...
    combine_downup(angle90(MMS_Shock_Down_Angle),angle90(MMS_Shock_Up_Angle)),'MMS')

plot_histogram({'MMS Shock Geometry','[\theta]'},0:15:90,...
    combine_downup(angle90(get_data(type1,'downstreamShockAngle')),angle90(get_data(type1,'upstreamShockAngle'))),name1,...
    combine_downup(angle90(get_data(type2,'downstreamShockAngle')),angle90(get_data(type2,'upstreamShockAngle'))),name2)

plot_histogram({'MMS Normalized Shock Geometry','[\theta]'},0:20:90,...
    combine_downup(angle90(get_data(type1,'downstreamShockAngle')),angle90(get_data(type1,'upstreamShockAngle'))),name1,...
    combine_downup(angle90(get_data(type2,'downstreamShockAngle')),angle90(get_data(type2,'upstreamShockAngle'))),name2,...
    combine_downup(angle90(MMS_Shock_Down_Angle),angle90(MMS_Shock_Up_Angle)))
    
% plot_scatter({'Downstream ThetaBn','Upstream ThetaBn'},...
%     angle90(get_data(type1,'downstreamShockAngle')),angle90(get_data(type1,'upstreamShockAngle')),name1,...
%     angle90(get_data(type2,'downstreamShockAngle')),angle90(get_data(type2,'upstreamShockAngle')),name2)

plot_histogram({'MMS Normalized E-TD Angle','[\theta]'},90,...
    combine_downup(get_data(type1b,'E-CS Angle Before'),180-get_data(type1b,'E-CS Angle After')),name1,...
    combine_downup(get_data(type2b,'E-CS Angle Before'),180-get_data(type2b,'E-CS Angle After')),name2,...
    [MMS_E_Down_Angle(MMS_magneticShear > 30);180-MMS_E_Up_Angle(MMS_magneticShear > 30)])


plot_histogram({'OMNI MMS Cone Angle','[\theta]'},30,...
    Omni_ConeAngle_MMStimes,'Omni',...
    MMS_ConeAngle,'MMS')

plot_histogram({'MMS Cone Angle','[\theta]'},30,...
    combine_downup(get_data(type1,'Cone Angle Before'),get_data(type1,'Cone Angle After')),name1,...
    combine_downup(get_data(type2,'Cone Angle Before'),get_data(type2,'Cone Angle After')),name2)

plot_histogram({'MMS Normalized Cone Angle','[\theta]'},30,...
    combine_downup(get_data(type1,'Cone Angle Before'),get_data(type1,'Cone Angle After')),name1,...
    combine_downup(get_data(type2,'Cone Angle Before'),get_data(type2,'Cone Angle After')),name2,...
    MMS_ConeAngle)

% plot_histogram({'MMS Normalized Alfvén Mach Number','[M_A]'},4:4:24,...
%     rmoutliers(combine_downup(get_data(type1,'Mach Number Before'),get_data(type1,'Mach Number After')),'percentiles',[0 93]),name1,...
%     rmoutliers(combine_downup(get_data(type2,'Mach Number Before'),get_data(type2,'Mach Number After')),'percentiles',[0 93]),name2,...
%     rmoutliers(MMS_MachNumber,'thresholdfactor',8))
% 
% 
% plot_histogram({'MMS Normalized Alfvén Mach Number at BS','[M_A]'},4:4:24,...
%     rmoutliers(combine_downup(get_data(type1,'Mach Number Before').*get_data(type1,'ClosestBSNx'),get_data(type1,'Mach Number After').*get_data(type1,'ClosestBSNx')),'percentiles',[0 93]),name1,...
%     rmoutliers(combine_downup(get_data(type2,'Mach Number Before').*get_data(type2,'ClosestBSNx'),get_data(type2,'Mach Number After').*get_data(type2,'ClosestBSNx')),'percentiles',[0 93]),name2,...
%     rmoutliers(MMS_MachNumber,'thresholdfactor',8))



% plot_histogram({'OMNI Normalized Alfvén Mach Number','[M_A]'},8:4:16,...
%     combine_downup(get_data(type1,'Omni Mach Number'),get_data(type1,'Omni Mach Number')),name1,...
%     combine_downup(get_data(type2,'Omni Mach Number'),get_data(type2,'Omni Mach Number')),name2,...
%     Omni_M)


% plot_histogram({'MMS Normalized Solar Wind Speed','[km/s]'},[300:100:650],...
%     rmoutliers(combine_downup(get_data(type1,'downstreamSpeed'),get_data(type1,'upstreamSpeed')),'percentiles',[0 97]),name1,...
%     rmoutliers(combine_downup(get_data(type2,'downstreamSpeed'),get_data(type2,'upstreamSpeed')),'percentiles',[0 97]),name2,...
%     rmoutliers(MMS_SWSpeed,'percentiles',[0 100]))
% 
% plot_histogram({'Velocity Deflection','[\theta]'},0:15:90,...
%     get_data(type1,'CoreSWVelocityDeflection'),name1,...
%     get_data(type2,'CoreSWVelocityDeflection'),name2)
% 
% plot_histogram({'MMS Normalized Solar Wind Density','[cm^{-3}]'},[5:5:20],...
%     rmoutliers(combine_downup(get_data(type1,'preN'),get_data(type1,'postN')),'percentiles',[0 100]),name1,...
%     rmoutliers(combine_downup(get_data(type2,'preN'),get_data(type2,'postN')),'percentiles',[0 96.5]),name2,...
%     rmoutliers((MMS_Nion),'percentiles',[0 100]))
% 
% plot_histogram({'MMS Normalized Magnetic Field Strength','[nT]'},2:2:10,...
%     rmoutliers(combine_downup(get_data(type1,'Bpre_cs'),get_data(type1,'Bpost_cs')),'thresholdfactor',3),name1,...
%     rmoutliers(combine_downup(get_data(type2,'Bpre_cs'),get_data(type2,'Bpost_cs')),'thresholdfactor',3),name2,...
%     rmoutliers(MMS_BMag,'thresholdfactor',3))

%Expansion Speed
% % plot_histogram({'Inner Edges Expansion','[km/s]'},-200:200:600,...
% %     rmoutliers(get_data(type1,'BoundaryExpansionC'),'thresholdfactor',3),name1,...
% %     rmoutliers(get_data(type2,'BoundaryExpansionC'),'thresholdfactor',3),name2) 


%Durations
plot_histogram({'Core Duration','[s]'},0:20:60,...
    rmoutliers(get_data(type1,'Core Duration')),name1,...
    rmoutliers(get_data(type2,'Core Duration')),name2)

plot_histogram({'Total Duration','[s]'},0:15:120,...
    get_data(type1,'Duration'),name1,...
    get_data(type2,'Duration'),name2)

%% Expansion Speed 8/24
%Current Sheet, -x for earthward
n_cs = -[get_data('HFA','CS_x'),get_data('HFA','CS_y'),get_data('HFA','CS_z')];
%TD-driven Events Leading
n_leading = [get_data('HFA','Leading Timing Normal_x'),get_data('HFA','Leading Timing Normal_y'),get_data('HFA','Leading Timing Normal_z')];
v_leading = get_data('HFA','Leading Timing Speed');
leading_velocity = n_leading.*v_leading;



%Calculate angle between n_Cs and leading edge
ncs_leading_angle = [];

for i=1:size(n_cs,1)
    ncs_leading_angle(i) = angle(n_cs(i,:),leading_velocity(i,:));
    
    if (ncs_leading_angle(i) < 90) && (n_leading(i,1) > 0 ) %If Less than 90 degrees, make nleading negative (same as n_cs)
        
        n_leading(i,:) = - n_leading(i,:);
        v_leading(i)   = - v_leading(i);
        
    elseif (ncs_leading_angle(i) > 90) && (n_leading(i,1) < 0 ) %If Greater than 90 degrees, making nleading positive (opposite of n_cs)
        n_leading(i,:) = - n_leading(i,:);
        v_leading(i)   = - v_leading(i);
        
    end
end

%Plot Angles between leading / core
plot_histogram({'n_{cs} and v_{leading} angle','[\theta]'},0:15:180,...
    ncs_leading_angle,'All')

%A = [ncs_leading_angle(ncs_leading_angle > 90)]
%Check
sum(ncs_leading_angle < 90)
sum(n_leading(:,1) < 0)

sum(ncs_leading_angle > 90)
sum(n_leading(:,1) > 0)

% %Trailing Edge
n_trailing = [get_data('HFA','Trailing Timing Normal_x'),get_data('HFA','Trailing Timing Normal_y'),get_data('HFA','Trailing Timing Normal_z')];
v_trailing = -abs(get_data('HFA','Trailing Timing Speed'));
trailing_velocity = n_trailing.*v_trailing;
% 
% % %Calculate angle between n_Cs and trailing edge
% ncs_trailing_angle = [];
% 
% for i=1:size(n_cs,1)
%     ncs_trailing_angle(i) = angle(n_cs(i,:),trailing_velocity(i,:));
%     
%     if (ncs_trailing_angle(i) < 90) && (n_trailing(i,1) < 0 ) %If Less than 90 degrees, make ntrailing positive (same as n_cs)
%         
%         n_trailing(i,:) = - n_trailing(i,:);
%         v_trailing(i)   = - v_trailing(i);
%         
%     elseif (ncs_trailing_angle(i) > 90) && (n_trailing(i,1) > 0 ) %If Greater than 90 degrees, making ntrailing negative (opposite of n_cs)
%         n_trailing(i,:) = - n_trailing(i,:);
%         v_trailing(i)   = - v_trailing(i);
%         
%     end
% end

%Check
% sum(ncs_trailing_angle < 90)
% sum(n_trailing(:,1) > 0)
% 
% sum(ncs_trailing_angle > 90)
% sum(n_trailing(:,1) < 0)
% 






%Solar Wind Speeds
Vdown = [get_data('HFA','Vdownx'),get_data('HFA','Vdowny'),get_data('HFA','Vdownz')];
Vup = [get_data('HFA','Vupx'),get_data('HFA','Vupy'),get_data('HFA','Vupz')];
% swVelocity = (Vdown+Vup)./2;

% %Choose the max  velocity from either side
% Vdown_norm = vecnorm(Vdown,1,2);
% Vup_norm = vecnorm(Vup,1,2);
% [~,index] = max([Vdown_norm,Vup_norm],[],2);
% swVelocity = [];
% for i=1:length(index)
%     if index(i) == 1 
%         swVelocity = [swVelocity;Vdown(i,:)];
%     else
%         swVelocity = [swVelocity;Vup(i,:)];
%     end
% end


%Expansion Speed Calculation
V1c = (v_leading - dot(Vdown,n_leading,2));
V2c = (v_trailing - dot(Vup,n_trailing,2));
BoundariesExpansionc = V1c + V2c; 
plot_histogram({'Inner Edges Expansion (Aug24)','[km/s]'},-350:100:600,...
    rmoutliers(BoundariesExpansionc,'thresholdfactor',3),'HFA,FB') 


%Event Numbers
EventNumber = get_data('HFA','Event Number');
SubstructurePresent = get_data('HFA','Substructure');
AllEvents = [EventNumber,...
    SubstructurePresent,...
    V1c,...
    BoundariesExpansionc,...
    n_leading,...
    v_leading,...
    Vdown,...
    n_trailing,...
    v_trailing];
ContractingEvents = [EventNumber(BoundariesExpansionc(:,1) < 0),...
    SubstructurePresent(BoundariesExpansionc(:,1) < 0),...
    BoundariesExpansionc(BoundariesExpansionc(:,1) < 0),...
    n_leading(BoundariesExpansionc(:,1) < 0,:),...
    v_leading(BoundariesExpansionc(:,1) < 0),...
    Vdown(BoundariesExpansionc(:,1) < 0,:)]; %19 Events


%% 
%% Histograms of Expansion Speeds
% 
% %swVelocity=[get_data('all','V_solarwindFramec_x'),get_data('all','V_solarwindFramec_y'),get_data('all','V_solarwindFramec_z')];
% LNx=get_data('all','Leading Timing Normal_x');
% LNy=get_data('all','Leading Timing Normal_y');
% LNz=get_data('all','Leading Timing Normal_z');
% LV =get_data('all','Leading Timing Speed');
% TNx=get_data('all','Trailing Timing Normal_x');
% TNy=get_data('all','Trailing Timing Normal_y');
% TNz=get_data('all','Trailing Timing Normal_z');
% TV =get_data('all','Trailing Timing Speed');
% 
% 
%             
% 
% Vdown = [get_data('all','Vdownx'),get_data('all','Vdowny'),get_data('all','Vdownz')];
% Vup = [get_data('allSS','Vupx'),get_data('all','Vupy'),get_data('all','Vupz')];
% % swVelocity = (Vdown+Vup)./2;
% 
% %Choose the max  velocity from either side
% Vdown_norm = vecnorm(Vdown,1,2);
% Vup_norm = vecnorm(Vup,1,2);
% [value,index] = max([Vdown_norm,Vup_norm],[],2);
% swVelocity = [];
% for i=1:length(index)
%     if index(i) == 1 
%         swVelocity = [swVelocity;Vdown(i,:)];
%     else
%         swVelocity = [swVelocity;Vup(i,:)];
%     end
% end
% 
% Leading_Normal = [LNx,LNy,LNz];
% Trailing_Normal = [TNx,TNy,TNz];
% Leading_Velocity = Leading_Normal.*LV;
% Trailing_Velocity = Trailing_Normal.*TV;
% 
% 
% 
% 
% 
% 
% 
% 
% % Diff_Velocity = Trailing_Velocity - Leading_Velocity;
% % 
% % swSpeed = swVelocity ./ vecnorm(swVelocity,1,2);
% % %Method 1
% % Leading_Minus_SW = Leading_Velocity - Vdown;
% % Trailing_Minus_SW = Trailing_Velocity - Vup;
% % Diff_Edges = dot((Leading_Minus_SW - Trailing_Minus_SW),-Leading_Normal,2)
% 
% 
% %Normal
% V1c = (LV - dot(swVelocity,[LNx,LNy,LNz],2));
% V2c = (TV - dot(swVelocity,[TNx,TNy,TNz],2));
% BoundariesExpansionc = V1c + V2c; %positive is contracting, negative is expanding
% plot_histogram({'Inner Edges Expansion (Normal)','[km/s]'},-400:200:600,...
%     rmoutliers(BoundariesExpansionc),'all') 
% 
% %Leading Normal negative
% for i=1:size(Leading_Normal,1)
%     if Leading_Normal(i,1) > 0
%         Leading_Normal = -Leading_Normal;
%         LV = -LV;
%     end
% end
% LV = -LV( Leading_Normal(:,1) > 0, :);
% Leading_Normal = -Leading_Normal( Leading_Normal(:,1) > 0 ,:);
% 
% V1c = (LV - dot(swVelocity,Leading_Normal,2));
% V2c = (TV - dot(swVelocity,Trailing_Normal,2));
% BoundariesExpansionc = V1c + V2c; %positive is contracting, negative is expanding
% plot_histogram({'Inner Edges Expansion (LNeg)','[km/s]'},-300:100:600,...
%     rmoutliers(BoundariesExpansionc),'all') 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % %Leading Substract
% % V1c = (LV - dot(swVelocity,[LNx,LNy,LNz],2));
% % V2c = (TV - dot(swVelocity,[TNx,TNy,TNz],2));
% % BoundariesExpansionc = -V1c + V2c; %positive is contracting, negative is expanding
% % plot_histogram({'Inner Edges Expansion (LSub)','[km/s]'},-600:200:600,...
% %     BoundariesExpansionc,'all') 
% % %Trailing Normal negative
% % V1c = (LV - dot(swVelocity,[LNx,LNy,LNz],2));
% % V2c = (TV - dot(swVelocity,[-TNx,-TNy,-TNz],2));
% % BoundariesExpansionc = V1c + V2c; %positive is contracting, negative is expanding
% % plot_histogram({'Inner Edges Expansion (TNeg)','[km/s]'},-600:200:600,...
% %     BoundariesExpansionc,'all') 
% % %Trailing Substract
% % V1c = (LV - dot(swVelocity,[LNx,LNy,LNz],2));
% % V2c = (TV - dot(swVelocity,[TNx,TNy,TNz],2));
% % BoundariesExpansionc = V1c + -V2c; %positive is contracting, negative is expanding
% % plot_histogram({'Inner Edges Expansion (TSub)','[km/s]'},-600:200:600,...
% %     BoundariesExpansionc,'all') 
% % 
% % 
% % % swSpeedSS=get_data(type1,'V_solarwindFramec');
% % % LNx_SS=get_data(type1,'Leading Timing Normal_x');
% % % LNy_SS=get_data(type1,'Leading Timing Normal_y');
% % % LNz_SS=get_data(type1,'Leading Timing Normal_z');
% % % LV_SS =get_data(type1,'Leading Timing Speed');
% % % TNx_SS=get_data(type1,'Trailing Timing Normal_x');
% % % TNy_SS=get_data(type1,'Trailing Timing Normal_y');
% % % TNz_SS=get_data(type1,'Trailing Timing Normal_z');
% % % TV_SS =get_data(type1,'Trailing Timing Speed');
% % % 
% % % 
% % % 
% % % swSpeedNS=get_data(type1,'V_solarwindFramec');
% % % LNx_NS=get_data(type1,'Leading Timing Normal_x');
% % % LNy_NS=get_data(type1,'Leading Timing Normal_y');
% % % LNz_NS=get_data(type1,'Leading Timing Normal_z');
% % % LV_NS =get_data(type1,'Leading Timing Speed');
% % % TNx_NS=get_data(type1,'Trailing Timing Normal_x');
% % % TNy_NS=get_data(type1,'Trailing Timing Normal_y');
% % % TNz_NS=get_data(type1,'Trailing Timing Normal_z');
% % % TV_NS =get_data(type1,'Trailing Timing Speed');
% % 
% % 
% % 
% % % 
% % % 
% % % 
% % % 
% % % 

%%
plot_histogram({'MMS Normalized Transit Speed','[km/s]'},0:150:500,...
    rmoutliers(get_data(type1b,'Transversal Speed'),'thresholdfactor',5),name1,...
    rmoutliers(get_data(type2b,'Transversal Speed'),'thresholdfactor',5),name2,...
    MMS_Transversal_Speed(MMS_magneticShear > 30 & MMS_Transversal_Speed<600))



%Age
plot_histogram({'Age','[s]'},0:30:150,...
    rmoutliers(get_data(type1b,'Age'),'thresholdfactor',3),name1,...
    rmoutliers(get_data(type2b,'Age'),'thresholdfactor',3),name2)

%Size
plot_histogram({'Length','[Re]'},0:1:4,...
    rmoutliers(get_data(type1b,'Size'),'thresholdfactor',3),name1,...
   rmoutliers(get_data(type2b,'Size'),'thresholdfactor',3),name2)

% MaxNsigmaAll = get_data('all','MaxNSigma');
% %Substructures
% plot_histogram({'Substructure Max Density','[\sigma]'},0:1:8,...
%     MaxNsigmaAll(MaxNsigmaAll>3.0),name1,...
%     MaxNsigmaAll(MaxNsigmaAll<3.0),name2)

plot_histogram({'Substructure Max Density','[\sigma]'},0:1:8,...
    get_data('all','MaxNSigma'),'')
%%
type1='SS';

plot_histogramSS({'Substructure Density Ratio','[#]'},-0.25:0.5:4.5,...
    rmoutliers(get_data(type1,'SSCoreNRatio'),'thresholdfactor',10),'Core',...
    rmoutliers(get_data(type1,'SSSWNRatio'),'thresholdfactor',10),'SW')

plot_histogramSS({'Substructure Magnetic Field Strength Ratio','[#]'},-0.25:0.5:4.75,...
    rmoutliers((get_data(type1,'SSCoreBRatio')),'thresholdfactor',3),'Core',...
    rmoutliers((get_data(type1,'SSSWBRatio')),'thresholdfactor',3),'SW')

plot_histogram({'Substructure Velocity Deflection','[\theta]'},0:15:90,...
    get_data(type1,'SSCoreDeflection'),'Core',...
    get_data(type1,'SSSWDeflection'),'SW')

plot_scatter({'Magnetic Field Strength Ratio','Velocity Deflection'},...
    (get_data(type1,'SSSWBRatio')),(get_data(type1,'SSCoreDeflection')),'All')

plot_scatter({'Substructure Duration','n corr |B|'},...
    (get_data(type1,'SS Duration')),(get_data(type1,'SSNBCorr')),'All')
%VxBins = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.5,3]
% VxBins = [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,2.0,2.5,3]

plot_histogramSS({'Substructure Vx Ratio','[#]'},0:0.2:3,...
    (get_data(type1,'SSCoreVxRatio')),'Core',...
    (get_data(type1,'SSSWVxRatio')),'SW')
% histogram(get_data(type1,'SSCoreVxRatio'),VxBins)
% hold on
% histogram(get_data(type1,'SSSWVxRatio'),VxBins)

plot_histogramSS({'Substructure Dynamic Pressure Ratio','[#]'},0:0.25:5,...
    rmoutliers(get_data(type1,'SSCoreDPRatio'),'thresholdfactor',10),'Core',...
    rmoutliers(get_data(type1,'SSSWDPRatio'),'thresholdfactor',10),'SW')

plot_histogram2Bins({'Substructure Dynamic Pressure Ratio','[#]'},0:0.5:5,0:0.25:2,...
    rmoutliers(get_data(type1,'SSCoreDPRatio'),'thresholdfactor',2),'Core',...
    rmoutliers(get_data(type1,'SSSWDPRatio'),'thresholdfactor',3),'SW')

plot_histogram({'Substructure Electron Temperature Ratio','[#]'},0:0.25:2.5,...
    rmoutliers(get_data(type1,'maxTecoreratio'),'thresholdfactor',3),'Core',...
    rmoutliers(get_data(type1,'maxTeSWRatio'),'thresholdfactor',3),'SW')

plot_histogram({'Substructure n T Correlation','[#]'},-1.0:0.5:1.0,...
    get_data(type1,'SSenTRatio'),'Electron',...
    get_data(type1,'SSinTRatio'),'Ion')


%%
plot_histogram({'n |V|, n |B| Substructure Correlation','[#]'},-1:0.5:1,...
    get_data(type1,'SSNVCorr'),'n corr |V|',...
    get_data(type1,'SSNBCorr'),'n corr |B|')

plot_histogram({'Substructure Duration','[s]'},0.0:1:10,...
    rmoutliers(get_data(type1,'SS Duration'),'thresholdfactor',12),'Substructured Events',[],'')

plot_histogram({'Substructure Duration Ratio','[Total Core Duration]'},0:0.1:0.5,...
    get_data(type1,'SS Core Duration Ratio'),'Substructured Events',[],'')

plot_histogram({'Substructure Location in Event Core','[Position in Core]'},0:0.1:1,...
    get_data(type1,'SSStartFraction'),'Substructured Events',[],'')

plot_histogram({'Substructure Size','[km]'},0:250:2000,...
    rmoutliers((get_data(type1,'SSSizeBulkV'))),'Substructured Events',[],'')

plot_histogram({'Substructure Size in Ion Inertial Lengths','[\lambda_{i}]'},0:3:20,...
    rmoutliers(get_data(type1,'SSSizeinIonInertialLengths')),'Substructured Events',[],'')





% SSsizeUsingBulkV = get_data(type1,'SSSizeBulkV');
% SSionInertialLength = get_data(type1,'ionInertiallength');
% SSSizeinInertialLengths = SSsizeUsingBulkV./SSionInertialLength;
% 
% plot_histogram({'Substructure Size in Ion Inertial Lengths','[\lambda_{i}]'},0:4:20,...
%     rmoutliers(SSSizeinInertialLengths,'thresholdfactor',3),'Substructured Events',[],'')
% 
% 
% 
% 
% 

% plot_scatter({'Inner Edges Expansion','Inner Edges Angle'},...
%     (get_data('all','BoundaryExpansionC')),(get_data('all','Event Edges Angle')),'All')


% plot_scatter({'Core Duration','Substructure Duration'},...
%     (get_data(type1,'Core Duration')),(get_data(type1,'SS Duration')),'Substructured Events')
 A=5
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions
%Average sownstaream and upstream for each type
function [combined] = combine_downup(downstream,upstream)
    combined = mean([downstream;upstream],2); %combine both vectors into one, (use ;) %for average use (,)
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
        if abs((binRange(1)) - (-(binRange(3)-binRange(2))/2)) < 0.001;
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
            yAxisLabel = '# of Events normalized to MMS observations';%sum(contains(parameter,'Shock')) > 0 ||...
            if sum(contains(parameter,'Shear')) > 0 ||...
                    sum(contains(parameter,'E-TD')) > 0 ||...
                    sum(contains(parameter,'TD-BS')) > 0 ||...
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
%             Aerror = (A).^(-1/2)./sum(A);
%             Berror = (B).^(-1/2)./sum(B);
%             
%             Aerror = (A_normalized).^(-1).*((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
%             Berror = (B_normalized).^(-1).*((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);
            
            Aerror = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
            Berror = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);
            
            %             Aerror = Aerror./ (A./sum(A));
            %             Berror = Berror./ (B./sum(B));
            %             Aerror = A_normalized.*(A.^(-1) + sum(A).^(-1) ).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) + sum(B).^(-1) ).^(1/2);
        else
            A_normalized = A./sum(A);
            B_normalized =B./sum(B);
            %             Aerror = A_normalized.*(A.^(-1) +      C.^(-1)).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) +      C.^(-1)).^(1/2);
%             
            %Binomial Error Bars
%             dA = (((A_normalized).^(-1).*(1-A_normalized))./sum(A)).^(1/2); %Real Error
%             dB = (((B_normalized).^(-1).*(1-B_normalized))./sum(B)).^(1/2); %Real Error
            
            dA = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2); %Real Error
            dB = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2); %Real Error
            
            %             dA = dA ./ (A./sum(A));
            %             dB = dB ./ (B./sum(B));
            
            C_normalized = C./sum(C);
            dC = ((C_normalized.*(1-C_normalized))./sum(C)).^(1/2);
            
%             Aerror = A_normalized.*(dA + dC);
%             Berror = B_normalized.*(dB + dC);
%             
%             Aerror = (A./C).*(dA + dC);
%             Berror = (B./C).*(dB + dC);
%             
%             Aerror = (sum(A)./C).*(dA + dC);
%             Berror = (sum(B)./C).*(dB + dC);
            
            Aerror = (sum(A)./C).*dA + (sum(C).*A./C./C).*dC;
            Berror = (sum(B)./C).*dB + (sum(C).*B./C./C).*dC;
            A_normalized = A./(C);
            B_normalized =B./(C);
            
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
            if length(binRange) > 12
               xlim([binRange(1) binRange(end)+1*bin_size])
               %xlim([binRange(1) binRange(end)+1*bin_size])
            else
                xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
                
                xlim([binRange(1) binRange(end)+(binRange(3)-binRange(2))])
            end
        else
            if length(binRange) > 12
               xlim([binRange(1) binRange(end)+bin_size])
            else
            xticks(binRange)
            xlim([binRange(1) binRange(end)])
            end
        end
        
        %Labels
        if length(binRange) > 12
            xlabelCell = xticklabels;
            if max([data1,data2]) > binRange(end)
                xlabelCell(end) = strcat('\geq',xlabelCell(end));
            end
            
            xticklabels(xlabelCell);
        else
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
                if custom_bins == 1 && xlabels(1) ~= 0 && (max([data1;data2],[],'all') >= binRange(end)+bin_size && min([data1;data2],[],'all') <= binRange(1) )
                    xlabelCell = xticklabels;
                    xlabelCell(1) = strcat('\leq',xlabelCell(1));
                    xlabelCell(end) = strcat('\geq',xlabelCell(end));
                    xticklabels(xlabelCell);
                elseif custom_bins == 1 && (max(data1) >= binRange(end)+bin_size || max(data2) >= binRange(end)+bin_size )
                    xlabelCell = xticklabels;
                    xlabelCell(end) = strcat('\geq',xlabelCell(end));
                    xticklabels(xlabelCell);
                    
                elseif custom_bins == 1 && (min(data1) < binRange(1) || min(data2) < binRange(1) )
                    xlabelCell = xticklabels;
                   xlabelCell(1) = strcat('\leq',xlabelCell(1));
                    xticklabels(xlabelCell);
                end
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
            
            if maxTypeNormalized == 1 %Data 1 is larger, 
                
                %d2Overd1Ratio = length(data2)/length(data1);
 
                if isempty(data2)
                else
                    yyaxis left
                end
                errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                d1_ylim = get(gca,'ylim');
                ylim([0 d1_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Dataset 2
                    yyaxis right
                    errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                    d2_ylim = get(gca,'ylim');
                    d2Overd1Ratio = nansum(B_normalized./d2_ylim(2))/nansum(A_normalized./d1_ylim(2)); %Number is less than 1

                    ylim([0 d2Overd1Ratio*d2_ylim(2)])
                end
            else %Data 2 is larger 
                %d1Overd2Ratio = length(data2)/length(data1);

                yyaxis right
                errorbar(xlabels+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                d2_ylim = get(gca,'ylim');
                ylim([0 d2_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Data 1
                    yyaxis left
                    errorbar(xlabels+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                    d1_ylim = get(gca,'ylim');
                    d1Overd2Ratio = nansum(A_normalized./d1_ylim(2))/nansum(B_normalized./d2_ylim(2)); %Number is less than 1
                    
                    ylim([0 d1Overd2Ratio*d1_ylim(2)])
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
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2);
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
%     savefig(strcat(fileName));
    
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
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2);
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
function [] = plot_irregularHistogram(parameter,bin_size,data1,label1,data2,label2)
    screenSize = get(0,'ScreenSize');
    figure('Position',[0 screenSize(4) 350 350])

    set(gcf,'color','w');

    if size(bin_size) == 1

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
%         bin_size = binRange(3)-binRange(2);
%         if binRange(1) == -(binRange(3)-binRange(2))/2;
%             offsetXlabels = 1;
%         else
%             offsetXlabels = 0;
%         end
    %calculate midpoint of each bin
    for i=1:length(binRange)-1
        
        midBin(i) = (binRange(i+1) + binRange(i))/2
    end
    midBin(length(binRange)) = (binRange(end) - binRange(end-1))/2 + binRange(end)
    
    end
    
    %Normalization Initialization if MMS data in given
    if nargin == 7

    else
        C=1;
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
    
            Aerror = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2);
            Berror = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2);

        else
            
            dA = ((A_normalized.*(1-A_normalized))./sum(A)).^(1/2); %Real Error
            dB = ((B_normalized.*(1-B_normalized))./sum(B)).^(1/2); %Real Error
         
            
            C_normalized = C./sum(C);
            dC = ((C_normalized.*(1-C_normalized))./sum(C)).^(1/2);

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
                Aerror = Aerror(2:end);
                Berror = Berror(2:end);
            else
                Aerror = Aerror(2:end);
                Berror = Berror(2:end);
            end
        else %% Bin Sizes
%             binRange = [binRange(1) - bin_size, binRange];
        end
        xlabels = [binRange(2)-(binRange(3)-binRange(2)),binRange(2:end-1),binRange(end-1)+(binRange(3)-binRange(2))];
        xlabels = binRange;
        
        if length(C)==1
            
            bar(midBin,A_normalized',1.0,'edgecolor','none'); hold on
            if isempty(data2)
                
            else
                [XB,YB] = stairs(xlabels,B_normalized');
                XB = [XB;XB(end)];
                YB = [YB;YB(end)];
                
                plot(XB,YB,'linewidth',3.5);
            end
       
        end
        
        
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
%         if offsetXlabels == 1
%             xlabelCell{1} = '0';
%             xticklabels(xlabelCell);
%         else
%         end
        
        
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
            errorbar(xlabels,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
            
            errorbar(xlabels,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        else
            %Find the type with the greater number of events
            maxAerror = max(A_normalized + Aerror);
            maxBerror = max(B_normalized + Berror);
            [~,maxTypeNormalized] = max([maxAerror,maxBerror]);
            
            if maxTypeNormalized == 1 %Data 1 is larger, 
                
                %d2Overd1Ratio = length(data2)/length(data1);
 
                if isempty(data2)
                else
                    yyaxis left
                end
                errorbar(xlabels,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                d1_ylim = get(gca,'ylim');
                ylim([0 d1_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Dataset 2
                    yyaxis right
                    errorbar(xlabels,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                    d2_ylim = get(gca,'ylim');
                    d2Overd1Ratio = nansum(B_normalized./d2_ylim(2))/nansum(A_normalized./d1_ylim(2)); %Number is less than 1

                    ylim([0 d2Overd1Ratio*d2_ylim(2)])
                end
            else %Data 2 is larger 
                yyaxis right
                errorbar(xlabels,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                d2_ylim = get(gca,'ylim');
                ylim([0 d2_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Data 1
                    yyaxis left
                    errorbar(xlabels,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                    d1_ylim = get(gca,'ylim');
                    d1Overd2Ratio = nansum(A_normalized./d1_ylim(2))/nansum(B_normalized./d2_ylim(2)); %Number is less than 1
                    
                    ylim([0 d1Overd2Ratio*d1_ylim(2)])
                end  
            end     
        end
        hold off
    else
    end

    if length(C) == 1
        ylimits = ylim;
        if ylimits(1) < 0
            ylim([0 ylimits(2)])
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
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2);
    colormap('winter');
    if sum(contains(parameter(2),'[Re]')) > 0
        parameter(2) = {'[R_E]'};
    elseif sum(contains(parameter(2),'[Re^2]')) > 0
        parameter(2) = {'[{R_E}^2]'};
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
    Database_Directory = '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/EventParameterDatabase_All_Events_2021_3.xlsx';%_BeforeChsngingLEadingEdgeV-.xlsx';
    
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
        case 'FB SS'
            Page_number = 11;
        case 'FB NS'
            Page_number = 12;            
        case 'all'
            Page_number = 1;
    end
    
    
    switch Event_Parameter
        case 'Event Number' 
            Column_number = 'A';
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
            
        case 'V_solarwindFramec_x'
            Column_number = 'FN';
        case 'V_solarwindFramec_y'
            Column_number = 'FO';
        case 'V_solarwindFramec_z'
            Column_number = 'FP';
            
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
            
        case 'Vdownx'
            Column_number = 'HR';
            
        case 'Vdowny'
            Column_number = 'HS';
            
        case 'Vdownz'
            Column_number = 'HT';
            
        case 'Vupx'
            Column_number = 'HU';
            
        case 'Vupy'
            Column_number = 'HV';
            
        case 'Vupz'
            Column_number = 'HW';
            
        case 'BoundaryExpansion2020'
            Column_number = 'HX';
            
        case 'LeadingSpeed2020'
            Column_number = 'HY';
            
        case 'TrailingSpeed2020'
            Column_number = 'HZ';
            
        case 'ClosestBSNx'
            Column_number = 'DE';
        
        case 'ClosestBSNy'
            Column_number = 'DF';
            
        case 'ClosestBSNz'
            Column_number = 'DG';
            
        case 'Downstream BS MA'
            Column_number = 'ID';
            
        case 'Upstream BS MA'
            Column_number = 'IE';
            
        case 'OMNI MGSMA'
            Column_number = 'IF';
            
        case 'Downstream BS MGSMA'
            Column_number = 'IG';
            
        case 'Upstream BS MGSMA'
            Column_number = 'IH';
            
        case 'Omni Ma'
            Column_number = 'II';
            
        case 'Omni Ms'
            Column_number = 'IJ';
            
        case 'Omni Ma BS'
            Column_number = 'IK';
            
        case 'Omni Ms BS'
            Column_number = 'IL';
            
        case 'Omni thetaBn'
            Column_number = 'IM';
            
        case 'Omni coneAngle'
            Column_number = 'IN';
            
            
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


function [] = plot_histogramSS(parameter,bin_size,data1,label1,data2,label2,data3)
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
        if abs((binRange(1)) - (-(binRange(3)-binRange(2))/2)) < 0.001;
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
            yAxisLabel = '# of Events normalized to MMS observations';%sum(contains(parameter,'Shock')) > 0 ||...
            if sum(contains(parameter,'Shear')) > 0 ||...
                    sum(contains(parameter,'E-TD')) > 0 ||...
                    sum(contains(parameter,'TD-BS')) > 0 ||...
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
        A = histcounts(data1,[-Inf binRange ]);
        B = histcounts(data2,[-Inf binRange ]);
        
        A_normalized = A./C;
        B_normalized =B./C;
        
        if length(C) == 1
            A_normalized = A./length(data1);
            B_normalized =B./length(data2);
            Aerror = (A).^(-1/2)./length(data1);
            Berror = (B).^(-1/2)./length(data2);
            
            Aerror = (A_normalized).^(-1).*((A_normalized.*(1-A_normalized))./length(data1)).^(1/2);
            Berror = (B_normalized).^(-1).*((B_normalized.*(1-B_normalized))./length(data2)).^(1/2);
            
            Aerror = ((A_normalized.*(1-A_normalized))./length(data1)).^(1/2);
            Berror = ((B_normalized.*(1-B_normalized))./length(data2)).^(1/2);
            
            %             Aerror = Aerror./ (A./sum(A));
            %             Berror = Berror./ (B./sum(B));
            %             Aerror = A_normalized.*(A.^(-1) + sum(A).^(-1) ).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) + sum(B).^(-1) ).^(1/2);
        else
            
            Aerror = A_normalized.*(A.^(-1) +      C.^(-1)).^(1/2);
            Berror = B_normalized.*(B.^(-1) +      C.^(-1)).^(1/2);
            
            %Binomial Error Bars
            dA = (((A_normalized).^(-1).*(1-A_normalized))./length(data1)).^(1/2); %Real Error
            dB = (((B_normalized).^(-1).*(1-B_normalized))./length(data2)).^(1/2); %Real Error
            
            dA = ((A_normalized.*(1-A_normalized))./length(data1)).^(1/2); %Real Error
            dB = ((B_normalized.*(1-B_normalized))./length(data2)).^(1/2); %Real Error
            
            %             dA = dA ./ (A./sum(A));
            %             dB = dB ./ (B./sum(B));
            
            C_normalized = C./sum(C);
            dC = ((C_normalized.*(1-C_normalized))./sum(C)).^(1/2);
            
%             Aerror = A_normalized.*(dA + dC);
%             Berror = B_normalized.*(dB + dC);
%             
%             Aerror = (A./C).*(dA + dC);
%             Berror = (B./C).*(dB + dC);
%             
%             Aerror = (sum(A)./C).*(dA + dC);
%             Berror = (sum(B)./C).*(dB + dC);
            
            Aerror = (length(data1)./C).*dA + (sum(C).*A./C./C).*dC;
            Berror = (length(data2)./C).*dB + (sum(C).*B./C./C).*dC;
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
            
            bar(xlabels(1:end-1)+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
            if isempty(data2)
                
            else
                [XB,YB] = stairs(xlabels(1:end-1),B_normalized');
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
            if length(binRange) > 12
               xlim([binRange(1) binRange(end)+1*bin_size])
               %xlim([binRange(1) binRange(end)+1*bin_size])
            else
                xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
                
                xlim([0 binRange(end)+(binRange(3)-binRange(2))])
            end
        else
            if length(binRange) > 12
               xlim([binRange(1) binRange(end)+bin_size])
            else
            xticks(binRange)
            xlim([binRange(1) binRange(end)])
            end
        end
        
        %Labels
% % % % % % % %         if length(binRange) > 12
% % % % % % % %             xlabelCell = xticklabels;
% % % % % % % %             
% % % % % % % %             xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %             xticklabels(xlabelCell);
% % % % % % % %         else
% % % % % % % %             if isempty(data2)
% % % % % % % %                 if custom_bins == 1 && xlabels(1) ~= 0 && (max(data1) >= binRange(end)  )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(1) = strcat('\leq',xlabelCell(1));
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 elseif custom_bins == 1 && (max(data1) >= binRange(end)  )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 end
% % % % % % % %             else
% % % % % % % %                 if custom_bins == 1 && xlabels(1) ~= 0 && (max([data1;data2],[],'all') >= binRange(end)+bin_size && min([data1;data2],[],'all') <= binRange(1) )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(1) = strcat('\leq',xlabelCell(1));
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 elseif custom_bins == 1 && (max(data1) >= binRange(end)+bin_size || max(data2) >= binRange(end)+bin_size )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                     
% % % % % % % %                 elseif custom_bins == 1 && (min(data1) < binRange(1) || min(data2) < binRange(1) )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                    xlabelCell(1) = strcat('\leq',xlabelCell(1));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %         end
    xlabelCell = xticklabels
%         if offsetXlabels == 1
%            % xlabelCell{1} = '0';
%             %xlabelCell{end} = [];
%             xticklabels(xlabelCell);
%         else
%         end
        xlim([0 binRange(end)])
        
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
            errorbar(xlabels(1:end-1)+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
            
            errorbar(xlabels(1:end-1)+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        else
            %Find the type with the greater number of events
            maxAerror = max(A_normalized + Aerror);
            maxBerror = max(B_normalized + Berror);
            [~,maxTypeNormalized] = max([maxAerror,maxBerror]);
            
            if maxTypeNormalized == 1 %Data 1 is larger, 
                
                %d2Overd1Ratio = length(data2)/length(data1);
 
                if isempty(data2)
                else
                    yyaxis left
                end
                errorbar(xlabels(1:end-1)+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                d1_ylim = get(gca,'ylim');
                ylim([0 d1_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Dataset 2
                    yyaxis right
                    errorbar(xlabels(1:end-1)+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                    d2_ylim = get(gca,'ylim');
                    d2Overd1Ratio = nansum(B_normalized./d2_ylim(2))/nansum(A_normalized./d1_ylim(2)); %Number is less than 1

                    ylim([0 d2Overd1Ratio*d2_ylim(2)])
                end
            else %Data 2 is larger 
                %d1Overd2Ratio = length(data2)/length(data1);

                yyaxis right
                errorbar(xlabels(1:end-1)+8*bin_size/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                d2_ylim = get(gca,'ylim');
                ylim([0 d2_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Data 1
                    yyaxis left
                    errorbar(xlabels(1:end-1)+8*bin_size/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                    d1_ylim = get(gca,'ylim');
                    d1Overd2Ratio = nansum(A_normalized./d1_ylim(2))/nansum(B_normalized./d2_ylim(2)); %Number is less than 1
                    
                    ylim([0 d1Overd2Ratio*d1_ylim(2)])
                end
                
            end
            xlim([0 binRange(end)])
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
        
        xlim([0 binRange(end)])
        
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

function [] = plot_histogram2Bins(parameter,bin_size1,bin_size2,data1,label1,data2,label2,data3)
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
    
%     if size(bin_size1) == 1
        
%         %     data1=rmoutliers(data1,'ThresholdFactor',25);
%         %     data2=rmoutliers(data2,'ThresholdFactor',25);
%         if min([data1;data2])  >= 0 && max([data1;data2]) <= 1
%             binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
%         elseif min([data1;data2])  >= -1 && max([data1;data2]) <= 1
%             binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
%         elseif min([data1;data2])  >= 0 && max([data1;data2]) <= 75
%             binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
%         elseif min([data1;data2])  >= 0 && max([data1;data2]) <= 180
%             binRange = 0:bin_size:180;
%         elseif  min([data1;data2]) >= -90 && max([data1;data2]) <= 90
%             binRange = -90:bin_size:90;
%         else
%             binRange = floor(min([data1;data2])):bin_size:ceil(max([data1;data2]));
%         end

%     else
        binRange1=bin_size1;
        custom_bins = 1;
        bin_size1 = binRange1(3)-binRange1(2);
        if abs((binRange1(1)) - (-(binRange1(3)-binRange1(2))/2)) < 0.001;
            offsetXlabels1 = 1;
        else
            offsetXlabels1 = 0;
        end
        

        binRange2=bin_size2;
        custom_bins = 1;
        bin_size2 = binRange2(3)-binRange2(2);
        if abs((binRange2(1)) - (-(binRange2(3)-binRange2(2))/2)) < 0.001;
            offsetXlabels2 = 1;
        else
            offsetXlabels2 = 0;
        end
        
%     end
    
    %Normalization Initialization if MMS data in given
%     if nargin == 7
%         %         C = histcounts(data3,[-Inf binRange max([data1;data2])]);
%         %Normalize to Hourly, MMS Survey is 4.5-s, 800 datapoints per hour
%         
%         
%         
%         C = histcounts(data3,[-Inf binRange Inf]);
%         if sum(contains(parameter,'Solar Wind Speed')) > 0 ||...
%                 sum(contains(parameter,'Solar Wind Density')) > 0 ||...
%                 sum(contains(parameter,'Magnetic Field Strength')) > 0 ||...
%                 sum(contains(parameter,'Mach Number')) > 0 ||...
%                 sum(contains(parameter,'Cone Angle')) > 0
%             C = C./800;
%             
%         end
%         
%         %         if sum(contains(parameter,'ClosestDistance')) > 0
%         %             %2 minutes per data point, per hour is 30 times.
%         %             C = C./30;
%         %         end
%         %         C = C./sum(C);
%         Cerror = C.^(-1/2);
%         Cerror(Cerror == inf) = NaN;
%         Cerror(Cerror == 0) = NaN;
%         if sum(contains(parameter,'MMS')) > 0
%             yAxisLabel = '# of Events normalized to MMS observations';%sum(contains(parameter,'Shock')) > 0 ||...
%             if sum(contains(parameter,'Shear')) > 0 ||...
%                     sum(contains(parameter,'E-TD')) > 0 ||...
%                     sum(contains(parameter,'TD-BS')) > 0 ||...
%                     sum(contains(parameter,'Transit')) > 0 ||...
%                     sum(contains(parameter,'Closest')) > 0
%                 yAxisLabel = '# of Events normalized to MMS solar wind discontinuities';
%             else
%                 yAxisLabel = '# of Events per hour';
%                 %                 Ctotal = sum(C);
%                 %                 HoursOfC = Ctotal /(60*60);
%                 %                 C = C./HoursOfC;
%                 
%             end
%         else
%             yAxisLabel = '# of Events normalized to WIND observations';
%         end
%     else
%         C=1;
%         Cerror = 1;
%         yAxisLabel = '# of Events / Total # Events';
%     end
    
C=1;
Cerror = 1;
yAxisLabel = '# of Events / Total # Events';
    %Plotting W
    if nargin ~= 4 %&& nargin ~= 3%Wewighting from Spacecraft Data Over All Time Range
        A = histcounts(data1,[-Inf binRange1 ]);
        B = histcounts(data2,[-Inf binRange2 ]);
        
        A_normalized = A./C;
        B_normalized =B./C;
        
        if length(C) == 1
            A_normalized = A./length(data1);
            B_normalized =B./length(data2);
            Aerror = (A).^(-1/2)./length(data1);
            Berror = (B).^(-1/2)./length(data2);
            
            Aerror = (A_normalized).^(-1).*((A_normalized.*(1-A_normalized))./length(data1)).^(1/2);
            Berror = (B_normalized).^(-1).*((B_normalized.*(1-B_normalized))./length(data2)).^(1/2);
            
            Aerror = ((A_normalized.*(1-A_normalized))./length(data1)).^(1/2);
            Berror = ((B_normalized.*(1-B_normalized))./length(data2)).^(1/2);
            
            %             Aerror = Aerror./ (A./sum(A));
            %             Berror = Berror./ (B./sum(B));
            %             Aerror = A_normalized.*(A.^(-1) + sum(A).^(-1) ).^(1/2);
            %             Berror = B_normalized.*(B.^(-1) + sum(B).^(-1) ).^(1/2);
        else
            
            Aerror = A_normalized.*(A.^(-1) +      C.^(-1)).^(1/2);
            Berror = B_normalized.*(B.^(-1) +      C.^(-1)).^(1/2);
            
            %Binomial Error Bars
            dA = (((A_normalized).^(-1).*(1-A_normalized))./length(data1)).^(1/2); %Real Error
            dB = (((B_normalized).^(-1).*(1-B_normalized))./length(data2)).^(1/2); %Real Error
            
            dA = ((A_normalized.*(1-A_normalized))./length(data1)).^(1/2); %Real Error
            dB = ((B_normalized.*(1-B_normalized))./length(data2)).^(1/2); %Real Error
            
            %             dA = dA ./ (A./sum(A));
            %             dB = dB ./ (B./sum(B));
            
            C_normalized = C./sum(C);
            dC = ((C_normalized.*(1-C_normalized))./sum(C)).^(1/2);
            
%             Aerror = A_normalized.*(dA + dC);
%             Berror = B_normalized.*(dB + dC);
%             
%             Aerror = (A./C).*(dA + dC);
%             Berror = (B./C).*(dB + dC);
%             
%             Aerror = (sum(A)./C).*(dA + dC);
%             Berror = (sum(B)./C).*(dB + dC);
            
            Aerror = (length(data1)./C).*dA + (sum(C).*A./C./C).*dC;
            Berror = (length(data2)./C).*dB + (sum(C).*B./C./C).*dC;
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
            binRange1 = [binRange1(1) - bin_size1, binRange1];
            binRange2 = [binRange2(1) - bin_size2, binRange2];
        end
        xlabels1 = [binRange1(2)-(binRange1(3)-binRange1(2)),binRange1(2:end-1),binRange1(end-1)+(binRange1(3)-binRange1(2))];
        xlabels2 = [binRange2(2)-(binRange2(3)-binRange2(2)),binRange2(2:end-1),binRange2(end-1)+(binRange2(3)-binRange2(2))];

        
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
            
            bar(xlabels1(1:end-1)+bin_size1/2,A_normalized',1.0,'edgecolor','none'); hold on
            if isempty(data2)
                
            else
                [XB,YB] = stairs(xlabels2(1:end-1),B_normalized');
                XB = [XB;XB(end)+bin_size2];
                YB = [YB;YB(end)];
                
                plot(XB,YB,'linewidth',3.5);
            end
            
            
        else
            if isempty(data2)
                
            else
                yyaxis left
            end
            bar(xlabels1+bin_size/2,A_normalized',1.0,'edgecolor','none'); hold on
            
            if isempty(data2)
                
            else
                [XB,YB] = stairs(xlabels2,B_normalized');
                XB = [XB;XB(end)+bin_size2];
                YB = [YB;YB(end)];
                
                yyaxis right
                plot(XB,YB,'linewidth',3.5); % plot(XB,55/45.*YB,'linewidth',3.5,'color','r');
                set(gca,'ycolor','r')
            end
            
        end
        
        %         stairs(xlabels,B_normalized','linewidth',3.5);
%         
%         if A(end) ~= 0 || B(end) ~= 0
%             if length(binRange) > 12
%                xlim([binRange(1) binRange(end)+1*bin_size])
%                %xlim([binRange(1) binRange(end)+1*bin_size])
%             else
%                 xticks([binRange binRange(end)+(binRange(3)-binRange(2))])
%                 
%                 xlim([0 binRange(end)+(binRange(3)-binRange(2))])
%             end
%         else
%             if length(binRange) > 12
%                xlim([binRange(1) binRange(end)+bin_size])
%             else
%             xticks(binRange)
%             xlim([binRange(1) binRange(end)])
%             end
%         end
    xlabelCell = xticklabels
if binRange2(end)> binRange1(end)
    xlim([binRange2(1) binRange2(end)])
    xlim([0 binRange2(end)])
else
    xlim([binRange1(1) binRange1(end)])
    xlim([0 binRange1(end)])
end
        
        %Labels
% % % % % % % %         if length(binRange) > 12
% % % % % % % %             xlabelCell = xticklabels;
% % % % % % % %             
% % % % % % % %             xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %             xticklabels(xlabelCell);
% % % % % % % %         else
% % % % % % % %             if isempty(data2)
% % % % % % % %                 if custom_bins == 1 && xlabels(1) ~= 0 && (max(data1) >= binRange(end)  )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(1) = strcat('\leq',xlabelCell(1));
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 elseif custom_bins == 1 && (max(data1) >= binRange(end)  )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 end
% % % % % % % %             else
% % % % % % % %                 if custom_bins == 1 && xlabels(1) ~= 0 && (max([data1;data2],[],'all') >= binRange(end)+bin_size && min([data1;data2],[],'all') <= binRange(1) )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(1) = strcat('\leq',xlabelCell(1));
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 elseif custom_bins == 1 && (max(data1) >= binRange(end)+bin_size || max(data2) >= binRange(end)+bin_size )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                     xlabelCell(end) = strcat('\geq',xlabelCell(end));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                     
% % % % % % % %                 elseif custom_bins == 1 && (min(data1) < binRange(1) || min(data2) < binRange(1) )
% % % % % % % %                     xlabelCell = xticklabels;
% % % % % % % %                    xlabelCell(1) = strcat('\leq',xlabelCell(1));
% % % % % % % %                     xticklabels(xlabelCell);
% % % % % % % %                 end
% % % % % % % %             end
% % % % % % % %         end

%         if offsetXlabels == 1
%            % xlabelCell{1} = '0';
%             %xlabelCell{end} = [];
%             xticklabels(xlabelCell);
%         else
%         end
        
        
        %Find best placement for legend
        [maxA,Ia] = max(A_normalized(A_normalized~=inf));
        [maxB,Ib] = max(B_normalized(B_normalized~=inf));
        if maxA > maxB
            I=Ia;
        else
            I = Ib;
        end
        
%         if I > length(binRange)/2
            legend(label1,label2,'FontSize',9,'Location','Northeast')
%         else
%             legend(label1,label2,'FontSize',9,'Location','Northeast')
%         end
        
        if isempty(data2)
            legend off
        else
        end
        hold on
        
        if length(C) == 1
            errorbar(xlabels1(1:end-1)+8*bin_size1/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
            
            errorbar(xlabels2(1:end-1)+8*bin_size2/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
        else
            %Find the type with the greater number of events
            maxAerror = max(A_normalized + Aerror);
            maxBerror = max(B_normalized + Berror);
            [~,maxTypeNormalized] = max([maxAerror,maxBerror]);
            
            if maxTypeNormalized == 1 %Data 1 is larger, 
                
                %d2Overd1Ratio = length(data2)/length(data1);
 
                if isempty(data2)
                else
                    yyaxis left
                end
                errorbar(xlabels1(1:end-1)+8*bin_size1/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                d1_ylim = get(gca,'ylim');
                ylim([0 d1_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Dataset 2
                    yyaxis right
                    errorbar(xlabels2(1:end-1)+8*bin_size2/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                    d2_ylim = get(gca,'ylim');
                    d2Overd1Ratio = nansum(B_normalized./d2_ylim(2))/nansum(A_normalized./d1_ylim(2)); %Number is less than 1

                    ylim([0 d2Overd1Ratio*d2_ylim(2)])
                end
            else %Data 2 is larger 
                %d1Overd2Ratio = length(data2)/length(data1);

                yyaxis right
                errorbar(xlabels2(1:end-1)+8*bin_size2/16,B_normalized,Berror,'LineStyle','none','Color','red','LineWidth',1.5,'HandleVisibility','off')
                d2_ylim = get(gca,'ylim');
                ylim([0 d2_ylim(2)])
                if isempty(data2)
                    
                else %Embiggen Data 1
                    yyaxis left
                    errorbar(xlabels1(1:end-1)+8*bin_size1/16,A_normalized,Aerror,'LineStyle','none','Color','blue','LineWidth',1.5,'HandleVisibility','off')
                    d1_ylim = get(gca,'ylim');
                    d1Overd2Ratio = nansum(A_normalized./d1_ylim(2))/nansum(B_normalized./d2_ylim(2)); %Number is less than 1
                    
                    ylim([0 d1Overd2Ratio*d1_ylim(2)])
                end
                
            end
            xlim([0 binRange(end)])
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
        
        xlim([0 binRange(end)])
        
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
    
    fileName = strcat(cell2mat(parameter(1)),'_',label1,'_',label2);
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
