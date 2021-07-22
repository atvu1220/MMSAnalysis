%Test Particle Simulation of Foreshock Ion at a Tangential Discontinuity
clear
% close all

%% Simulation input
%Constant
b0 = 5e-9; %nT
q = 1.6e-19; %C
mion = 1.6605e-27;%kg
mu0 = 4*pi*10^-7;
n0 = 5e6; %m^-3
c = 2.9999*10^8; %km/s
epsilon0 = 8.85e-12;

%Simulation Parameters
nx = 200;%inertial lengths
ny = 400;%inertial lengths
nz = 200;%inertial lengths
dx_frac = 1 ; %Fraction of inertial lengths
dt = 0.02;
ddthickness = 16;
thermal = 4;
T = 10;
n=3;
v0 = 0;
%Plasma Calculations
%plasma frequency
%plasma period T=2piOmega

omega_p = q*b0/mion;
lambda_i = c/sqrt(n0*q^2/(epsilon0*mion)); %meters
va = sqrt(b0^2 / (mu0*mion*n0)); %m/s
% dt = dt/(2*pi*2*omega_p);
T = 30;%T/(2*pi*omega_p);

%% Grid Initialization
%Magnetic Field Grid
grid_x = 0:lambda_i*dx_frac:(nx*lambda_i);
grid_y = 0:lambda_i*dx_frac:(ny*lambda_i);
grid_z = 0:lambda_i*dx_frac:(nz*lambda_i);
% Bx_slice = 0.5*b0 + 0.5*b0*tanh( (grid_z(nz/2/dx_frac) - grid_z )     / ( lambda_i*dx_frac*ddthickness ) );
% By_slice = 0.5*b0 + 0.5*b0*tanh( (grid_z - grid_z(nz/2/dx_frac) )     / ( lambda_i*dx_frac*ddthickness ) );
% Bz_slice = zeros(1,length(Bx_slice));
B = zeros(length(grid_x),length(grid_y),length(grid_z),3);
for i=1:length(grid_x)
    for j=1:length(grid_y)
        for k=1:length(grid_z)
            %                  B(i,j,k,1) = b0;
            %                  B(i,j,k,2) = 0.0;
            %                  B(i,j,k,3) = 0.0;
            if k <= nz/2
                B(i,j,k,1) = 0.5*b0 + 0.5*b0*tanh( (grid_z(nz/2/dx_frac) - grid_z(k) )     / ( lambda_i*dx_frac*ddthickness ) );
                B(i,j,k,2) = 0.5*b0 - 0.5*b0*tanh( (grid_z(nz/2/dx_frac) - grid_z(k) )     / ( lambda_i*dx_frac*ddthickness ) );
                B(i,j,k,3) = 0.0;
            end
            if k > nz/2
                B(i,j,k,1) = 0.5*b0 - 0.5*b0*tanh( (grid_z(k) - grid_z(nz/2/dx_frac) )     / ( lambda_i*dx_frac*ddthickness ) );
                B(i,j,k,2) = 0.5*b0 + 0.5*b0*tanh( (grid_z(k) - grid_z(nz/2/dx_frac) )     / ( lambda_i*dx_frac*ddthickness ) );
                B(i,j,k,3) = 0.0;
            end
        end
    end
end
% plot(grid_z/(lambda_i),squeeze(B(nx/2,ny/2,:,1)));
% hold on
% plot(grid_z/(lambda_i),squeeze(B(nx/2,ny/2,:,2)));
% xlabel('Z')
% ylabel('B')
% title('B at x=nx/2, ny=ny/2')
% legend(['Bx';'By'])
% xlim([0 nz])

% vp = zeros(10000,2);
% for i=1:length(vp)
%     vp(i,1) = sqrt(-2*log(rand))*cos(2*pi*rand);
%     vp(i,2) = sqrt(-2*log(rand))*cos(2*pi*rand);
% end
% scatter(vp(:,1),vp(:,2))
% axis square
% figure

%% Particle Initialization
%Random Parameters
pos1 = [grid_x(nx/2) + 2*lambda_i*(rand(1,n)-0.5);...
    grid_y(ny/2)+ 2*lambda_i*(rand(1,n)-0.5);...
    grid_z(nz/2)+ -0.5*ddthickness*lambda_i + ddthickness*1*lambda_i*(rand(1,n)-0.5)];
vthx = va*sqrt(-2*log(rand(1,n))).*cos(2*pi*rand(1,n));
vthy = va*sqrt(-2*log(rand(1,n))).*cos(2*pi*rand(1,n));
vthz = va*sqrt(-2*log(rand(1,n))).*cos(2*pi*rand(1,n));


%Set Parameters
pos1 = [...
    grid_x(nx),grid_y(ny/2),grid_z(nz/4);...
    grid_x(nx),grid_y(ny/2),grid_z(nz/2)-0.5*16*lambda_i;...
    grid_x(nx),grid_y(ny/2),grid_z(nz/2)-0.5*16*lambda_i;...
    grid_x(nx),grid_y(ny/2),grid_z(nz/2)-0.5*16*lambda_i;...
    grid_x(nx),grid_y(ny/2),grid_z(nz/2)-0.5*16*lambda_i;...
    ]';
vthx = thermal.*[...
    va,...
    va,...
    -va,...
    va,...
    va,...
    ];
vthy = thermal.*[...
    va,... %Black
    va,... %Green
    -va,...%Red
    2*va,...%Blue
    4*va,...%Cyan
    ];
vthz = thermal.*[...
    va,...
    va,...
    -va,...
    va,...
    4*va,...
    ];
n=3;



for nPart=1:n
    [Bpos,~] = get_FieldfromPos(pos1(:,nPart),grid_x,grid_y,grid_z,B,va);
    vel1(:,nPart) = [-24*va*Bpos(1)/(norm(Bpos)) + vthx(:,nPart); -24*va*Bpos(2)/(norm(Bpos)) + vthy(:,nPart);vthz(:,nPart)];
    %vel1(:,nPart) = [0 + vthx(:,nPart); + vthy(:,nPart);vthz(:,nPart)];
end

if v0 == 0
    vel0 = va.*ones(size(vel1));
else
    vel0=vel1;
end

%ExB speed
Ex = 0;
Ey = 0;
Ez = -(12*va*B(:,:,:,2));

ExBx = -(Ez.*B(:,:,:,2))./( (B(:,:,:,1).^2 + B(:,:,:,2).^2)   );
ExBy = +(Ez.*B(:,:,:,1))./( (B(:,:,:,1).^2 + B(:,:,:,2).^2)   );
ExBz = zeros(size(ExBx));

ExB(:,:,:,1) = ExBx;
ExB(:,:,:,2) = ExBy;
ExB(:,:,:,3) = ExBz;

ExBpos0 = [];

for nPart=1:n
    [ExBpos,~] =  get_FieldfromPos(pos1(:,nPart),grid_x,grid_y,grid_z,ExB,va)
    vel1(:,nPart) = vel1(:,nPart) + ExBpos;
    ExBpos0 = [ExBpos0,ExBpos];
end

% ExBpos0 = ExBpos;
Bx = squeeze(B(nx/2,ny/2,:,1));
By = squeeze(B(nx/2,ny/2,:,2));
Btot = vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2);
%Plot References
figure('color','white','position',[300 1600 400 1200]);
subplot(3,1,1);
plot(-24*squeeze(B(nx/2,ny/2,:,1))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)));hold on
plot(-24*squeeze(B(nx/2,ny/2,:,2))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)));
xline(nz/2-0.5*16)
title('Vpara speed','fontsize',14);legend({'Vparax';'Vparay'},'fontsize',14,'location','northwest');
xlim([0 nz]);xlabel('Z');ylabel('v/va'); grid on

subplot(3,1,2)
plot(squeeze(ExB(nx/2,ny/2,:,1))./va);hold on
plot(squeeze(ExB(nx/2,ny/2,:,2))./va);
plot(squeeze(ExB(nx/2,ny/2,:,3))./va);
xline(nz/2-0.5*16)
title('ExB speed','fontsize',14);legend({'ExBx';'ExBy';'ExBz'},'fontsize',14,'location','northwest');
xlim([0 nz]);xlabel('Z');ylabel('v/va'); grid on

subplot(3,1,3);
plot(-24*squeeze(B(nx/2,ny/2,:,1))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)) + squeeze(ExB(nx/2,ny/2,:,1))./va);hold on
plot(-24*squeeze(B(nx/2,ny/2,:,2))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)) + squeeze(ExB(nx/2,ny/2,:,2))./va);
xline(nz/2-0.5*16)
title('Sum','fontsize',14);legend({'Sumx';'Sumy'},'fontsize',14,'location','northwest');
xlim([0 nz]);xlabel('Z');ylabel('v/va'); grid on

%% Loop
%Particle Initialization
%     subplot(1,2,1)
%     scatter(vthx,vthz);
%     xlabel('Vx');
%     ylabel('Vz');
%     hold on
% scatterFig = figure;


timeFig = figure('color','w','position',[1500 0 800 1600]);
zFig = figure('color','w','position',[700 0 800 1600]);

vminus = zeros(3,n);

partColor = rand(3,n);
partColor = [0 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 1]';

zpos = zeros(T/dt,n);
vperpx = zeros(T/dt,n);
vperpy = zeros(T/dt,n);
vx = zeros(T/dt,n);
vy = zeros(T/dt,n);
Bx = zeros(T/dt,n);
By = zeros(T/dt,n);
for i=0:dt:T
    for j=1:n
        if pos1(1,j) < 0
            continue
        end
        [Bpos,Epos] = get_FieldfromPos(pos1(:,j),grid_x,grid_y,grid_z,B,va);
        Bpos = Bpos';
        Epos = Epos';
        %Velocity Update
        vminus(:,j) = vel1(:,j) + q*Epos*dt/(mion*2);
        
        t = q*Bpos/mion * dt/2;
        vprime = vminus(:,j) + (cross(vminus(:,j),t))';
        s = 2*t / (1+norm(t)^2);
        vplus(:,j) = vminus(:,j) + (cross(vprime,s))';
        
        vel2 = vplus(:,j) + q*Epos/mion * dt/2;
        
        %Position Update
        pos2 = pos1(:,j) + vel2*dt;
        
        if pos2(2) < 0
            pos2(2) = pos2(2) + max(grid_y);
        elseif pos2(2) > max(grid_y)
            pos2(2) = pos2(2) - max(grid_y);
        end
        
        if mod(i,5*dt) == 0
            %Calculate velocity components relative to local B
            [vparax,vparay,vparaz,vperpx,vperpy,vperpz] = get_Bcomp(Bpos,vel2);
            if v0 == 0
                vperpx0 = va;vperpy0 = va;
            else
                [vparax0,vparay0,vparaz0,vperpx0,vperpy0,vperpz0] = get_Bcomp(Bpos,vel0(:,j));
            end
            %% Time
            figure(timeFig);
            subplot(4,2,[1:1]);hold on;
            plot3(([i-0.25,i]./(2*2*pi)),([pos1(2,j) pos2(2)])/lambda_i,([pos1(3,j) pos2(3)])/lambda_i,'o','color',partColor(:,j)','markersize',2,'markerfacecolor',partColor(:,j))     
             xlim([0 T]./(2*2*pi));ylim([0 ny]);%zlim([25 125]);
            xlabel('\Omega_{pi}');ylabel('Y');zlabel('Z');
            title('Time vs. FS Ion Z Position','fontsize',14);
            grid on; view([-0.09079754601227,0.05046728971962])
            
            subplot(4,2,[2:2]);hold on;
            plot3(([i-0.25,i]./(2*2*pi)),([pos1(2,j) pos2(2)])/lambda_i,([pos1(3,j) pos2(3)])/lambda_i,'o','color',partColor(:,j)','markersize',2,'markerfacecolor',partColor(:,j))
             xlim([0 T./(2*2*pi)]);ylim([0 ny]);%zlim([25 125]);
            xlabel('\Omega_{pi}');ylabel('Y');zlabel('Z');
            title('Time vs. FS Ion Z Position','fontsize',14);
            grid on; view([-0.09079754601227,0.05046728971962])
            
            %B
            subplot(4,2,3);hold on
            scatter(i./(2*2*pi),Bpos(1)./b0,[],partColor(:,j)','filled'); ylabel('Bx', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14);%z vs xcomp
            xlim([0 T./(2*2*pi)]); grid on
            %             scatter(pos2(3)./lambda_i,Bpos(1)./b0,[],partColor(:,j)','filled'); %z vs xcomp
            
            subplot(4,2,4);hold on
            scatter(i./(2*2*pi)./(2*pi),Bpos(2)./b0,[],partColor(:,j)','filled'); ylabel('By', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%z vs xcomp
            xlim([0 T./(2*2*pi)]); grid on
            %             scatter(pos2(3)./lambda_i,Bpos(2)./b0,[],partColor(:,j)','filled'); %z vs xcomp
           
            subplot(4,2,5);hold on
            scatter(i./(2*2*pi),(vel2(1)-ExBpos0(1,j))./vel0(1,j),[],partColor(:,j)','filled');ylabel('Vx/va', 'fontsize',14); xlabel('\Omega_{pi}', 'fontsize',14) %z vs xcomp
            xlim([0 T./(2*2*pi)]); grid on
            %             scatter(pos2(3)./lambda_i,vel2(1)./vel0(1,j),[],partColor(:,j)','filled'); ylabel('Vx/Vx0')%z vs xcomp
            
            subplot(4,2,6); hold on
            scatter(i./(2*2*pi),(vel2(2)-ExBpos0(2,j))./vel0(2,j),[],partColor(:,j)','filled');  ylabel('Vy/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%y vs ycomp
            xlim([0 T./(2*2*pi)]); grid on
            %             scatter(pos2(3)./lambda_i,vel2(2)./vel0(2,j),[],partColor(:,j)','filled');  ylabel('Vy/Vy0')%y vs ycomp
       
            subplot(4,2,7);hold on            
            scatter(i./(2*2*pi),(vperpx-ExBpos0(1,j))./vperpx0,[],partColor(:,j)','filled');  ylabel('Vperpx/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%z vs perpxcomp
            xlim([0 T./(2*2*pi)]); grid on
            %             scatter(pos2(3)./lambda_i,vperpx./vperpx0,[],partColor(:,j)','filled'); ylabel('Vperpx/Vperpx0')%z vs perpxcomp

            subplot(4,2,8); hold on            
            scatter(i./(2*2*pi),(vperpy-ExBpos0(2,j))./vperpy0,[],partColor(:,j)','filled');ylabel('Vperpy/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14) %z vs perpycomp
            xlim([0 T./(2*2*pi)]); grid on
            %             scatter(pos2(3)./lambda_i,vperpy./vperpy0,[],partColor(:,j)','filled'); ylabel('Vperpy/Vperpy0')%z vs perpycomp
     
            
            zpos(i,j) = pos2(3);
            vperpx(i,j) = (vperpx-ExBpos0(1,j))./va;
            vperpy(i,j) = (vperpy-ExBpos0(2,j))./va;
            vx(i,j) = (vel2(1)-ExBpos0(1,j))./va;
            vy(i,j) = (vel2(2)-ExBpos0(2,j))./va;
            Bx(i,j) = Bpos(1)./b0;
            By(i,j) = Bpos(2)./b0;
            %% Z
            figure(zFig);
            subplot(5,2,[1:4]);hold on
            plot3(([pos1(1,j) pos2(1)])/lambda_i,([pos1(2,j) pos2(2)])/lambda_i,([pos1(3,j) pos2(3)])/lambda_i,'o','color',partColor(:,j)','markersize',8,'markerfacecolor',partColor(:,j))      
%             xlim([0 nx]);ylim([0 ny]);zlim([25 125])       
            xlabel('X');  ylabel('Y');zlabel('Z')  
            title('FS Ion Trajectory', 'fontsize',14)
            grid on; view([-0.09079754601227,0.05046728971962])
           
            %B
            subplot(5,2,5); hold on
            scatter(pos2(3)./lambda_i,Bpos(1)./b0,[],partColor(:,j)','filled');  ylabel('Bx', 'fontsize',14);xlabel('Z', 'fontsize',14)%z vs xcomp
            xline(nz/2-0.5*16)
%             scatter(i,Bpos(1)./b0,[],partColor(:,j)','filled'); %z vs xcomp
           
            subplot(5,2,6); hold on
            scatter(pos2(3)./lambda_i,Bpos(2)./b0,[],partColor(:,j)','filled'); ylabel('By', 'fontsize',14);xlabel('Z', 'fontsize',14)%z vs xcomp
            xline(nz/2-0.5*16)
%             scatter(i,Bpos(2)./b0,[],partColor(:,j)','filled'); %z vs xcomp
            
            subplot(5,2,7); hold on
            scatter(pos2(3)./lambda_i,(vel2(1)-ExBpos0(1,j))./vel0(1,j),[],partColor(:,j)','filled'); ylabel('Vx/va', 'fontsize',14);xlabel('Z', 'fontsize',14)%z vs xcomp
            xline(nz/2-0.5*16)
%             scatter(i,vel2(1)./vel0(1,j),[],partColor(:,j)','filled'); ylabel('Vx/Vx0')%z vs xcomp       

            subplot(5,2,8); hold on
            scatter(pos2(3)./lambda_i,(vel2(2)-ExBpos0(2,j))./vel0(2,j),[],partColor(:,j)','filled');   ylabel('Vy/va', 'fontsize',14);xlabel('Z', 'fontsize',14)%y vs ycomp
            xline(nz/2-0.5*16)
%             scatter(i,vel2(2)./vel0(2,j),[],partColor(:,j)','filled'); ylabel('Vy/Vy0')%y vs ycomp           
            
            subplot(5,2,9); hold on
            scatter(pos2(3)./lambda_i,(vperpx-ExBpos0(1,j))./vperpx0,[],partColor(:,j)','filled'); ylabel('Vperpx/va', 'fontsize',14);xlabel('Z', 'fontsize',14) %z vs perpxcomp
            xline(nz/2-0.5*16)
%             scatter(i,vperpx./vperpx0,[],partColor(:,j)','filled');ylabel('Vperpx/Vperpx0') %z vs perpxcomp         
            
            subplot(5,2,10);hold on;
            scatter(pos2(3)./lambda_i,(vperpy-ExBpos0(2,j))./vperpy0,[],partColor(:,j)','filled'); ylabel('Vperpy/va', 'fontsize',14);xlabel('Z', 'fontsize',14) %z vs perpycomp
            xline(nz/2-0.5*16)
%             scatter(i,vperpy./vperpy0,[],partColor(:,j)','filled'); ylabel('Vperpy/Vperpy0') %z vs perpycomp         
           
            
        end

        %         elseif mod(i,0.5) == 0
        %           plot3(([pos1(1) pos2(1)])/lambda_i,([pos1(2) pos2(2)])/lambda_i,([pos1(3) pos2(3)])/lambda_i,'|','color',partColor(:,j)','markersize',10,'markerfacecolor',partColor(:,j),'linewidth',10)
        %
        %             xlim([0 nx])
        %             ylim([0 ny])
        %             zlim([0 nz])
        %             xlabel('X')
        %             ylabel('Y')
        %             zlabel('Z')
        %             title('Foreshock Ion Trajectory at TD')
        %             hold on
        %             grid on
        %         end
        
%         if pos2(1) < 0
%             break;
%         end
        
        
        %Push down Velocity
        vel1(:,j) = vel2;
        pos1(:,j) = pos2;
    end
%     drawnow
    hold off
end









%
% figure
% [X,Y,Z] = meshgrid(grid_y,grid_x,grid_z);
% quiver3(X,Y,Z,squeeze(B(:,:,:,1)),squeeze(B(:,:,:,2)),squeeze(B(:,:,:,3)))
% xlabel('X'); ylabel('Y'); zlabel('Z')


%Particle Update

function[Bpos,Epos] =  get_FieldfromPos(pos,grid_x,grid_y,grid_z,B,va)

Bpos = squeeze(B(find(grid_x < pos(1),1,'last'),find(grid_y < pos(2),1,'last'),find(grid_z < pos(3),1,'last'),:));
sw_vel = [12*va,0,0]';
if (size(Bpos)) ~= (size(sw_vel))
    pause
end

Epos = - cross(sw_vel,Bpos);
Epos = Epos';
end

function[vparax,vparay,vparaz,vperpx,vperpy,vperpz] = get_Bcomp(Bpos,vel)
% vparax = Bpos(1).*vel(1)./(norm(Bpos));
% vparay = Bpos(2).*vel(2)./(norm(Bpos));
% vparaz = Bpos(3).*vel(3)./(norm(Bpos));


vpara = Bpos(1).*vel(1) + Bpos(2).*vel(2) + Bpos(3).*vel(3);
vparax = vpara.*Bpos(1)./norm(Bpos).^2;
vparay = vpara.*Bpos(2)./norm(Bpos).^2;
vparaz = vpara.*Bpos(3)./norm(Bpos).^2;
vperpx = vel(1) - vparax;
vperpy = vel(2) - vparay;
vperpz = vel(3) - vparaz;
end
