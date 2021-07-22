% function [v_HT] = calculate_vHT(fgm_timedata, fgm_bdata,...
%         fpi_i_timedata,fpi_i_ndata,fpi_i_vdata)
%
%
%
%
% ['','']
% end
%Walen Test

clear
%Program Parameters
event_start =  '2018-01-29 03:36:15.000';
event_end = '2018-01-29 03:36:17.800';
plot='1';

event_start =  '2017-12-29 19:11:17.000'
event_end = '2017-12-29 19:11:40.000'

event_start =  '2017-12-29 19:11:18.000'
event_end = '2017-12-29 19:11:20.000'

if plot=='1'
    cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/WalenTests'
    eventFolderName = strcat(event_start(1:10),'/',event_start(12:13),'-',event_start(15:16));
    if ~exist(eventFolderName,'dir')
        mkdir(eventFolderName)
    end
    cd(eventFolderName)
end
%%%%MMS1
%load FGM data
[fgm_timedata_raw, fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'brst');
%Load FPI_i
[fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,~,~,~,~] = load_fpi(event_start,event_end,1,'brst','i');


% % %%%%MMS2
% % figure
% % %load FGM data
% % [fgm_timedata_raw, fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'brst');
% % %Load FPI_i
% % [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,~,~,~,~] = load_fpi(event_start,event_end,2,'brst','i');
% %
% % plot_fgm_magnetic(event_start,event_end,fgm_timedata_raw,fgm_bdata_raw,3,1); datetick('keeplimits')
% % plot_fpi_number(event_start,event_end,fpi_i_timedata,fpi_i_ndata,3,2,'i'); datetick('keeplimits')
% % plot_fpi_bulkv(event_start,event_end,fpi_i_timedata,fpi_i_vdata,3,3); datetick('keeplimits')
% %
% %
% % %%%%MMS3
% % %load FGM data
% % [fgm_timedata_raw, fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'brst');
% % %Load FPI_i
% % [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,~,~,~,~] = load_fpi(event_start,event_end,3,'brst','i');
% %
% % plot_fgm_magnetic(event_start,event_end,fgm_timedata_raw,fgm_bdata_raw,3,1); datetick('keeplimits')
% % plot_fpi_number(event_start,event_end,fpi_i_timedata,fpi_i_ndata,3,2,'i'); datetick('keeplimits')
% % plot_fpi_bulkv(event_start,event_end,fpi_i_timedata,fpi_i_vdata,3,3); datetick('keeplimits')
% %
% %
% % %%%%MMS4
% % %load FGM data
% % [fgm_timedata_raw, fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'brst');
% % %Load FPI_i
% % [fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,~,~,~,~] = load_fpi(event_start,event_end,4,'brst','i');
% %
% % plot_fgm_magnetic(event_start,event_end,fgm_timedata_raw,fgm_bdata_raw,3,1); datetick('keeplimits')
% % plot_fpi_number(event_start,event_end,fpi_i_timedata,fpi_i_ndata,3,2,'i'); datetick('keeplimits')
% % plot_fpi_bulkv(event_start,event_end,fpi_i_timedata,fpi_i_vdata,3,3); datetick('keeplimits')











% %Interpolate data with fpi instrument, FGM down to FPI
fgm_bdata = zeros(length(fpi_i_timedata),3);
for i=1:4
    fgm_bdata(:,i)=interp1(fgm_timedata_raw,fgm_bdata_raw(:,i),fpi_i_timedata,'pchip');
end

% % Crop Data to specified time period
[~,fpi_i_ndata,~,~] = crop(fpi_i_timedata,fpi_i_ndata,event_start,event_end); %ndata
[~,fpi_i_vdata,~,~] = crop(fpi_i_timedata,fpi_i_vdata,event_start,event_end); %vdata
[fpi_i_timedata,fgm_bdata,~,~] = crop(fpi_i_timedata,fgm_bdata,event_start,event_end); %Bdata

figure
plot_fgm_magnetic(event_start,event_end,fpi_i_timedata,fgm_bdata,3,1); datetick('keeplimits')
title(strcat('MMS1- ',event_start(1:19)),'FontSize',14)
% plot_fgm_magnetic(event_start,event_end,fpi_i_timedata_cropped,fgm_bdata_raw,3,1);  datetick('keeplimits');
plot_fpi_number(event_start,event_end,fpi_i_timedata,fpi_i_ndata,3,2,'i'); datetick('keeplimits')
plot_fpi_bulkv(event_start,event_end,fpi_i_timedata,fpi_i_vdata,3,3); datetick('keeplimits')

if plot=='1'
    plot_name =  strcat('WalenTest_','mms','1','_BnV_',event_start(1:19),'_',event_end(1:19),'.eps');
    print(gcf,'-depsc2', '-loose', plot_name)
end



%Test from page 227
% % % % %Calculate VHT for all data points
% % % % %Calculate V_A for each data point
% % % % B = [-13.6 -24.7 54.6;...
% % % %     -14.8 -24.9 58.7;...
% % % %     -13.4 -17.2 62.4;...
% % % %     -14.0 -25.0 43.8;...
% % % %     -7.1 -4.5 33.5;...
% % % %     -0.9 -5.0 44.4;...
% % % %     -10.0 -0.4 44.6;...
% % % %     -6.1 -4.8 21.1;...
% % % %     1.2 1.6 21.0;...
% % % %     -3.4 -3.9 4.1;...
% % % %     -0.9 1.2 5.0;...
% % % %     -1.0 -1.5 12.3;...
% % % %     11.0 13.2 29.7;...
% % % %     19.1 34.4 20.1;...
% % % %     24.9 50.1 1.9;...
% % % %     29.2 47.1 -10.6];
% % % %
% % % %
% % % % B(:,4) = sqrt(B(:,1).^2 + B(:,2).^2 + B(:,3).^2)
% % % %
% % % % V = [ -111.0 -211.0 57.0;...
% % % % -102.0 -213.0 41.0;...
% % % % -111.0 -196.0 34.0;...
% % % % -135.0 -229.0 54.0;...
% % % % -128.0 -252.0 54.0;...
% % % % -82.0 -237.0 51.0;...
% % % % -139.0 -228.0 77.0;...
% % % % -143.0 -241.0 57.0;...
% % % % -132.0 -226.0 80.0;...
% % % % -117.0 -217.0 79.0;...
% % % % -112.0 -210.0 93.0;...
% % % % -98.0 -212.0 92.0;...
% % % % -90.0 -186.0 104.0;...
% % % % -83.0 -168.0 121.0;...
% % % % -82.0 -129.0 88.0;...
% % % % -93.0 -123.0 53.0];
% % % %
% % % % N = [12.18,9.22,9.92,18.08,20.39,15.00,20.19,23.53,24.31,25.91,26.17,24.49,22.20,22.86,17.56,17.86]';
% % % %
% % % %
% % % % %VHT Algorithm
% % % %
% % % % KV = zeros(3,1);
% % % % K0= zeros(3,3);
% % % % for i=1:16
% % % %    K = calculate_K(B(i,:));
% % % %    K0 = K0 + K;
% % % %    KV = KV + K*V(i,1:3)';
% % % % end
% % % %
% % % % KVT = KV./length(16); %Average
% % % % K0T = K0./length(16); %Average
% % % %
% % % % V_HT = inv(K0)*KV
% % % % V_HT = [-122.8;-223.0;76];
% % % %
% % % % VminusVHT = V-V_HT';
% % % % mu_0 = 4*pi*10^-7;
% % % % hplus_mass = 1.6726219e-27;
% % % % V_Alfven = B(:,1:3).*(10^-9).*(mu_0*hplus_mass*N.*10^-6).^(-1/2)./10^9; %in km/s
% % % %
% % % %
% % % % %Quality Test
% % % % figure
% % % % E_c = cross(-V,B(:,1:3));
% % % % E_HT = cross(-repmat(V_HT,1,length(B))',B(:,1:3));
% % % % scatter(E_HT(:,1),E_c(:,1));hold on
% % % % scatter(E_HT(:,2),E_c(:,2))
% % % % scatter(E_HT(:,3),E_c(:,3))
% % % % xlabel({'E_{HT}';'[\muV/m]'})
% % % % ylabel({'E_{c}';'[\muV/m]'})
% % % % CC_E_HT_E_C = corrcoef(E_HT(:,1),E_c(:,1))
% % % % title({strcat('HT-Frame Determination: Correlation Coefficient:',num2str(CC_E_HT_E_C(1,2)))})
% % % %
% % % % %Second Quality Test
% % % % D = calculate_residualE(V_HT,V,B);
% % % % D0 = calculate_residualE([0;0;0],V,B);
% % % % D_ratio = D/D0
% % % %
% % % % D2 = calculate_residualE2(V_HT,V,B)
% % % % D02 = calculate_residualE2([0;0;0],V,B);
% % % % D_ratio2 = D2/D02



% % %VHT Algorithm

KV = zeros(3,1);
K0= zeros(3,3);
for i=1:length(fpi_i_timedata)
    K = calculate_K(fgm_bdata(i,:));
    K0 = K0 + K;
    KV = KV + K*fpi_i_vdata(i,1:3)';
end

KV = KV./length(fpi_i_timedata); %Average
K0 = K0./length(fpi_i_timedata); %Average

V_HT = inv(K0)*KV;

VminusVHT = fpi_i_vdata-V_HT';
mu_0 = 4*pi*10^-7;
hplus_mass = 1.67e-27;
V_Alfven = fgm_bdata(:,1:3).*(10^-9).*(mu_0*hplus_mass*fpi_i_ndata.*10^-6).^(-1/2)./10^9; %in km/s

%Calculate V_Alfven in terms of B, V_Alfven para and V_Alfven Perp
B_V_A_Angle = atan2(norm(cross(fgm_bdata(:,1:3)',V_Alfven')),dot(fgm_bdata(:,1:3)',V_Alfven'))';
V_A_parallel = vecnorm(V_Alfven,1,2).*cos(B_V_A_Angle);
V_A_perp = vecnorm(V_Alfven,1,2).*sin(B_V_A_Angle);


%Calculate V-VHT in terms of B, VminusVHT para and VminusHT Perp
B_VminusVHT_Angle = atan2(norm(cross(fgm_bdata(:,1:3)',VminusVHT')),dot(fgm_bdata(:,1:3)',VminusVHT'))';
VminusVHT_parallel = vecnorm(VminusVHT,1,2).*cos(B_VminusVHT_Angle);
VminusVHT_perp = vecnorm(VminusVHT,1,2).*sin(B_VminusVHT_Angle);




%Plotting Walen Test
figure
scatter(V_Alfven(:,1),VminusVHT(:,1),'filled');hold on
scatter(V_Alfven(:,2),VminusVHT(:,2),'filled')
scatter(V_Alfven(:,3),VminusVHT(:,3),'filled')
axis([-max(V_Alfven,[],'all'),max(V_Alfven,[],'all'),-max(VminusVHT,[],'all'),max(VminusVHT,[],'all')])
xlabel({'V_{A}';'[km/s]'})
ylabel({'V-V_{HT}';'[km/s]'})
legend('V_x','V_y','V_z')
title({'Walen Test'},'FontSize',14)
if plot=='1'
    plot_name =  strcat('WalenTest_','mms','1','_VA_VHT_',event_start(1:19),'_',event_end(1:19),'.eps');
    print(gcf,'-depsc2', '-loose', plot_name)
end


%Plotting Walen Test with perp and parallel
figure
scatter(V_A_parallel,VminusVHT_parallel,'filled');hold on
scatter(V_A_perp,VminusVHT_perp,'filled')
axis([-max([V_A_parallel;V_A_perp],[],'all'),max([V_A_parallel;V_A_perp],[],'all'),-max([VminusVHT_parallel;VminusVHT_perp],[],'all'),max([VminusVHT_parallel;VminusVHT_perp],[],'all')])
xlabel({'V_{A}';'[km/s]'})
ylabel({'V-V_{HT}';'[km/s]'})
legend('V_{parallel}','V_{perp}')
title({'Walen Test'},'FontSize',14)



%Quality Test
figure
E_c = cross(-fpi_i_vdata,fgm_bdata(:,1:3));
E_HT = cross(-repmat(V_HT,1,length(fgm_bdata))',fgm_bdata(:,1:3));
scatter(E_HT(:,1),E_c(:,1),'filled');hold on
scatter(E_HT(:,2),E_c(:,2),'filled')
scatter(E_HT(:,3),E_c(:,3),'filled')
xlabel({'E_{HT}';'[\muV/m]'})
ylabel({'E_{c}';'[\muV/m]'})
legend('E_x','E_y','E_z')
CC_E_HT_E_C = corrcoef(E_HT(:,1),E_c(:,1))
title({strcat('HT-Frame Determination: Correlation Coefficient:',num2str(CC_E_HT_E_C(1,2)))},'FontSize',14)
if plot=='1'
    plot_name =  strcat('WalenTest_','mms','1','_Ec_EHT_',event_start(1:19),'_',event_end(1:19),'.eps');
    print(gcf,'-depsc2', '-loose', plot_name)
end

%Second Quality Test
D = calculate_residualE(V_HT,fpi_i_vdata,fgm_bdata);
D0 = calculate_residualE([0;0;0],fpi_i_vdata,fgm_bdata);
D_ratio = D/D0






function [D] = calculate_residualE(V,v,B)
    if length(V) ~= length(v)
        V = repmat(V,1,length(v))';
    end
    D = (1/length(B)) * sum(norm(cross((v-V),B(:,1:3)))^2);
end

% function [D] = calculate_residualE2(V,v,B)
%     if length(V) ~= length(v)
%         V = repmat(V,1,length(v))';
%     end
%
%     sum = [0;0;0];
%     for i=1:length(B)
%         diff = v-V;
%         crossprod= cross(diff,B(:,1:3));
%         mag = norm(crossprod);
%         sum = sum + mag^2;
%     end
%     D = sum/length(B);
% end

function [K] = calculate_K(bdata)
    K = [bdata(4)^2 - bdata(1)^2, -bdata(1)*bdata(2), -bdata(1)*bdata(3);...
        -bdata(2)*bdata(1), bdata(4)^2 - bdata(2)^2, -bdata(2)*bdata(3);...
        -bdata(3)*bdata(1), -bdata(3)*bdata(2), bdata(4)^2-bdata(3)^2];
end