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
dt = 0.05;
ddthickness = 16;
thermal = 4;
T = 155;
n=200;


omega_p = q*b0/mion;
lambda_i = c/sqrt(n0*q^2/(epsilon0*mion)); %meters
va = sqrt(b0^2 / (mu0*mion*n0)); %m/s

z0 = -0.25*ddthickness*lambda_i;
% %% Grid Initialization
% %Magnetic Field Grid
% grid_x = 0:lambda_i*dx_frac:(nx*lambda_i);
% grid_y = 0:lambda_i*dx_frac:(ny*lambda_i);
% grid_z = 0:lambda_i*dx_frac:(nz*lambda_i);
% % Bx_slice = 0.5*b0 + 0.5*b0*tanh( (grid_z(nz/2/dx_frac) - grid_z )     / ( lambda_i*dx_frac*ddthickness ) );
% % By_slice = 0.5*b0 + 0.5*b0*tanh( (grid_z - grid_z(nz/2/dx_frac) )     / ( lambda_i*dx_frac*ddthickness ) );
% % Bz_slice = zeros(1,length(Bx_slice));
% B = zeros(length(grid_x),length(grid_y),length(grid_z),3);
% for i=1:length(grid_x)
%     for j=1:length(grid_y)
%         for k=1:length(grid_z)
%             if k <= nz/2
%                 B(i,j,k,1) = 0.5*b0 + 0.5*b0*tanh( (grid_z(nz/2/dx_frac) - grid_z(k) )     / ( lambda_i*dx_frac*ddthickness ) );
%                 B(i,j,k,2) = 0.5*b0 - 0.5*b0*tanh( (grid_z(nz/2/dx_frac) - grid_z(k) )     / ( lambda_i*dx_frac*ddthickness ) );
%                 B(i,j,k,3) = 0.0;
%             end
%             if k > nz/2
%                 B(i,j,k,1) = 0.5*b0 - 0.5*b0*tanh( (grid_z(k) - grid_z(nz/2/dx_frac) )     / ( lambda_i*dx_frac*ddthickness ) );
%                 B(i,j,k,2) = 0.5*b0 + 0.5*b0*tanh( (grid_z(k) - grid_z(nz/2/dx_frac) )     / ( lambda_i*dx_frac*ddthickness ) );
%                 B(i,j,k,3) = 0.0;
%             end
%         end
%     end
% end
% 

%% Particle Initialization
x0 = 0;%nx*lambda_i;

%Random Parameters
% pos1 = [grid_x(nx/2) + 2*lambda_i*(rand(1,n)-0.5);...
%     grid_y(ny/2)+ 2*lambda_i*(rand(1,n)-0.5);...
%     grid_z(nz/2)+ -0.5*ddthickness*lambda_i + ddthickness*1*lambda_i*(rand(1,n)-0.5)];

pos1 = [zeros(1,n);zeros(1,n); z0.*ones(1,n)] ;


vthx = thermal.*va*sqrt(-2*log(rand(1,n))).*cos(2*pi*rand(1,n));
vthy = thermal.*va*sqrt(-2*log(rand(1,n))).*cos(2*pi*rand(1,n));
vthz = thermal.*va*sqrt(-2*log(rand(1,n))).*cos(2*pi*rand(1,n));


%Set Parameters
% pos1 = [...
%     x0,0,-4*ddthickness*lambda_i;...
%     x0,0,z0;...
%     x0,0,z0;...
%     x0,0,z0;...
%     x0,0,z0;...
%     ]';
% vthx = thermal.*[...
%     va,...
%     va,...
%     -va,...
%     va,...
%     va,...
%     ];
% vthy = thermal.*[...
%     va,... %Black
%     va,... %Green
%     -va,...%Red
%     2*va,...%Blue
%     4*va,...%Cyan
%     ];
% vthz = thermal.*[...
%     va,...
%     va,...
%     -va,...
%     va,...
%     4*va,...
%     ];



ExBpos0 = [];
for nPart=1:n
    [Bpos,Epos] = get_FieldfromPos(pos1(:,nPart),ddthickness,b0,va,lambda_i);
    ExBpos = cross(Epos,Bpos)./norm(Bpos).^2;
    vel1(:,nPart) = [-24*va*Bpos(1)/(norm(Bpos)) + vthx(:,nPart) + ExBpos(1); -24*va*Bpos(2)/(norm(Bpos)) + vthy(:,nPart) + ExBpos(2); vthz(:,nPart) + ExBpos(3) ];
%     vx = (-24*va*Bpos(1)/(norm(Bpos)) + ExBpos(1))./va;
%     vy = (-24*va*Bpos(2)/(norm(Bpos)) + ExBpos(2))./va;
    ExBpos0 = [ExBpos0,ExBpos];
end


%Check B Profile
testX = -8*ddthickness*lambda_i:lambda_i:8*ddthickness*lambda_i;
testBx = [];
testBy = [];

% for i=1:length(testX)
%     [Bpos,~] = get_FieldfromPos([0,0,testX(i)],ddthickness,b0,va,lambda_i);
%     testBx = [testBx,Bpos(1)];
%     testBy = [testBy,Bpos(2)];
% end
% figure; plot(testX/lambda_i,testBx); hold on; plot(testX/lambda_i,testBy);
% 
% Ez = -(12*va*testBy);
% ExBx = -Ez.*testBy./(testBx.^2 + testBy.^2);
% ExBy = +Ez.*testBx./(testBx.^2 + testBy.^2);
% 
% figure; plot(testX/lambda_i,ExBx./va); hold on; plot(testX/lambda_i,ExBy./va);

% %ExB speed
% Ex = 0;
% Ey = 0;
% Ez = -(12*va*B(:,:,:,2));
% 
% ExBx = -(Ez.*B(:,:,:,2))./( (B(:,:,:,1).^2 + B(:,:,:,2).^2)   );
% ExBy = +(Ez.*B(:,:,:,1))./( (B(:,:,:,1).^2 + B(:,:,:,2).^2)   );
% ExBz = zeros(size(ExBx));
% 
% ExB(:,:,:,1) = ExBx;
% ExB(:,:,:,2) = ExBy;
% ExB(:,:,:,3) = ExBz;

% ExBpos0 = [];
% 
% for nPart=1:n
% %     [Bpos,~] = get_FieldfromPos(pos1(:,nPart),ddthickness,b0,va);
% %     Ex = 0;
% %     Ey = 0;
% %     Ez = -(12*va*Bpos(2));
% %     
% %     ExBx = -(Ez.*Bpos(2))./( (Bpos(1).^2 + Bpos(2).^2)   );
% %     ExBy = +(Ez.*Bpos(1))./( (Bpos(1).^2 + Bpos(2).^2)   );
% %     ExBz = 0;
% %     
% %     ExBpos = [ExBx,ExBy,ExBz];
%     vel1(:,nPart) = vel1(:,nPart) + ExBpos;
%     ExBpos0 = [ExBpos0,ExBpos];
% end


%% Loop
timeFig = figure('color','w','position',[800 0 800 1600]);
zFig = figure('color','w','position',[1600 0 800 1600]);
avgFig = figure('color','w','position',[000 0 800 800]);
vminus = zeros(3,n);

partColor = rand(3,n);
% partColor = [0 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 1]';

time = zeros(T/dt,1);
xpos= zeros(T/dt,n);
zpos = zeros(T/dt,n);
vpx = zeros(T/dt,n);
vpy = zeros(T/dt,n);
vx = zeros(T/dt,n);
vy = zeros(T/dt,n);
Bx = zeros(T/dt,n);
By = zeros(T/dt,n);
for i=1:1:T/dt
    time(i) = i*dt/(4*pi);
    for j=1:n
%         if pos1(1,j) < 0
%             continue
%         end
        [Bpos,Epos] = get_FieldfromPos(pos1(:,j),ddthickness,b0,va,lambda_i);
%         Bpos = Bpos';
%         Epos = Epos';
        %Velocity Update
        vminus(:,j) = vel1(:,j) + q*Epos*dt/(mion*2);
        
        t = q*Bpos/mion * dt/2;
        vprime = vminus(:,j) + (cross(vminus(:,j),t));
        s = 2*t / (1+norm(t)^2);
        vplus(:,j) = vminus(:,j) + (cross(vprime,s));
        
        vel2 = vplus(:,j) + q*Epos/mion * dt/2;
        
        %Position Update
        pos2 = pos1(:,j) + vel2*dt;
        
%         if pos2(2) < 0
%             pos2(2) = pos2(2) + max(grid_y);
%         elseif pos2(2) > max(grid_y)
%             pos2(2) = pos2(2) - max(grid_y);
%         end

        [vparax,vparay,vparaz,vperpx,vperpy,vperpz] = get_Bcomp(Bpos,vel2);
        
        xpos(i,j) = pos2(1)/lambda_i;
        zpos(i,j) = pos2(3)/lambda_i;
        vpx(i,j) = (vperpx-ExBpos0(1,j))./va;
        vpy(i,j) = (vperpy-ExBpos0(2,j))./va;
        vx(i,j) = (vel2(1)-ExBpos0(1,j))./va;
        vy(i,j) = (vel2(2)-ExBpos0(2,j))./va;
        Bx(i,j) = Bpos(1)./b0;
        By(i,j) = Bpos(2)./b0;
        
        %Push down Velocity
        vel1(:,j) = vel2;
        pos1(:,j) = pos2;
    end
end
figure(avgFig);
avgvx = mean(vx,2);
avgvy = mean(vy,2);
avgvpx = mean(vpx,2);
avgvpy = mean(vpy,2);
sgtitle(sprintf('Average for particles originating at z=%2d',z0/lambda_i),'fontsize',14)
subplot(2,2,1);hold on
plot(time(2:end),avgvx(2:end),'color','k');ylabel('avg Vx/va', 'fontsize',14); xlabel('\Omega_{pi}', 'fontsize',14) %z vs xcomp
xlim([0 T./(2*2*pi)]); grid on

subplot(2,2,2); hold on
plot(time(2:end),avgvy(2:end),'color','k');  ylabel('avg Vy/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%y vs ycomp
xlim([0 T./(2*2*pi)]); grid on

subplot(2,2,3);hold on
plot(time(2:end),avgvpx(2:end),'color','k');  ylabel('avg Vperpx/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%z vs perpxcomp
xlim([0 T./(2*2*pi)]); grid on

subplot(2,2,4); hold on
plot(time(2:end),avgvpy(2:end),'color','k');ylabel('avg Vperpy/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14) %z vs perpycomp
xlim([0 T./(2*2*pi)]); grid on

for j=1:n
    %% Time
    figure(timeFig);
    subplot(4,2,[1:1]);hold on;
    plot(time(2:end),zpos(2:end,j),'o','color',partColor(:,j)','markersize',2,'markerfacecolor',partColor(:,j))
    xlim([0 T]./(2*2*pi));
    xlabel('\Omega_{pi}');ylabel('Z');
    title('Time vs. FS Ion Z Position','fontsize',14);
    grid on; 
    
    subplot(4,2,[2:2]);hold on;
    plot(time(2:end),zpos(2:end,j),'o','color',partColor(:,j)','markersize',2,'markerfacecolor',partColor(:,j))
    xlim([0 T./(2*2*pi)]);
    xlabel('\Omega_{pi}');ylabel('Z');
    title('Time vs. FS Ion Z Position','fontsize',14);
    grid on; 
    
    %B
    subplot(4,2,3);hold on
    plot(time(2:end),Bx(2:end,j),'color',partColor(:,j)'); ylabel('Bx', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14);%z vs xcomp
    xlim([0 T./(2*2*pi)]); grid on
    
    subplot(4,2,4);hold on
    plot(time(2:end),By(2:end,j),'color',partColor(:,j)'); ylabel('By', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%z vs xcomp
    xlim([0 T./(2*2*pi)]); grid on
    
    subplot(4,2,5);hold on
    plot(time(2:end),vx(2:end,j),'color',partColor(:,j)');ylabel('Vx/va', 'fontsize',14); xlabel('\Omega_{pi}', 'fontsize',14) %z vs xcomp
    xlim([0 T./(2*2*pi)]); grid on
    
    subplot(4,2,6); hold on
    plot(time(2:end),vy(2:end,j),'color',partColor(:,j)');  ylabel('Vy/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%y vs ycomp
    xlim([0 T./(2*2*pi)]); grid on
    
    subplot(4,2,7);hold on
    plot(time(2:end),vpx(2:end,j),'color',partColor(:,j)');  ylabel('Vperpx/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14)%z vs perpxcomp
    xlim([0 T./(2*2*pi)]); grid on
    
    subplot(4,2,8); hold on
    plot(time(2:end),vpy(2:end,j),'color',partColor(:,j)');ylabel('Vperpy/va', 'fontsize',14);xlabel('\Omega_{pi}', 'fontsize',14) %z vs perpycomp
    xlim([0 T./(2*2*pi)]); grid on
    
    %% Z
    figure(zFig);
    subplot(5,2,[1:4]);hold on
    plot(xpos(2:end,j),zpos(2:end,j),'o','color',partColor(:,j)','markersize',8,'markerfacecolor',partColor(:,j))
    xlabel('X');  ylabel('Z');
    title('FS Ion Trajectory', 'fontsize',14)
    grid on; 
    
    %B
    subplot(5,2,5); hold on
    plot(zpos(2:end,j),Bx(2:end,j),'color',partColor(:,j)');  ylabel('Bx', 'fontsize',14);xlabel('Z', 'fontsize',14)%z vs xcomp
    xline(z0./lambda_i)
    
    subplot(5,2,6); hold on
    plot(zpos(2:end,j),By(2:end,j),'color',partColor(:,j)'); ylabel('By', 'fontsize',14);xlabel('Z', 'fontsize',14)%z vs xcomp
    xline(z0./lambda_i)
    
    subplot(5,2,7); hold on
    plot(zpos(2:end,j),vx(2:end,j),'color',partColor(:,j)'); ylabel('Vx/va', 'fontsize',14);xlabel('Z', 'fontsize',14)%z vs xcomp
    xline(z0./lambda_i)
    
    subplot(5,2,8); hold on
    plot(zpos(2:end,j),vy(2:end,j),'color',partColor(:,j)');   ylabel('Vy/va', 'fontsize',14);xlabel('Z', 'fontsize',14)%y vs ycomp
    xline(z0./lambda_i)
    
    subplot(5,2,9); hold on
    plot(zpos(2:end,j),vpx(2:end,j),'color',partColor(:,j)'); ylabel('Vperpx/va', 'fontsize',14);xlabel('Z', 'fontsize',14) %z vs perpxcomp
    xline(z0./lambda_i)
    
    subplot(5,2,10);hold on;
    plot(zpos(2:end,j),vpy(2:end,j),'color',partColor(:,j)'); ylabel('Vperpy/va', 'fontsize',14);xlabel('Z', 'fontsize',14) %z vs perpycomp
    xline(z0./lambda_i)
end






function[Bpos,Epos] =  get_FieldfromPos(pos,ddthickness,b0,va,lambda_i)


if pos(3) <= 0
    Bpos(1) = 0.5*b0 - 0.5*b0*tanh( (pos(3) )     / ( lambda_i*ddthickness ) );
    Bpos(2) = 0.5*b0 + 0.5*b0*tanh( (pos(3))     / ( lambda_i*ddthickness ) );
    Bpos(3) = 0.0;
end
if pos(3) > 0
    Bpos(1) = 0.5*b0 - 0.5*b0*tanh( (pos(3) )     / ( lambda_i*ddthickness ) );
    Bpos(2) = 0.5*b0 + 0.5*b0*tanh( (pos(3) )     / ( lambda_i*ddthickness ) );
    Bpos(3) = 0.0;
end


% % Bpos(1) = b0/2;
% % Bpos(2) = b0/2;
% % if pos(3) < 4*ddthickness*lambda_i
% %     Bpos(1) = Bpos(1) + b0/2;
% %     Bpos(2) = Bpos(2) - b0/2;
% %     Bpos(3) = 0;
% % 
% % elseif pos(3) > 4*ddthickness*lambda_i
% %     Bpos(1) = Bpos(1) + b0/2;
% %     Bpos(2) = Bpos(2) - b0/2;
% %     Bpos(3) = 0;
% % else
% %     Bpos(1) = Bpos(1) - 0.5*(4*ddthickness*lambda_i)*pos(3);
% %     Bpos(2) = Bpos(2) + 0.5*(4*ddthickness*lambda_i)*pos(3);
% %     Bpos(3) = 0;
% % end

% Bpos = squeeze(B(find(grid_x < pos(1),1,'last'),find(grid_y < pos(2),1,'last'),find(grid_z < pos(3),1,'last'),:));


sw_vel = [12*va,0,0];
if (size(Bpos)) ~= (size(sw_vel))
    pause
end

Epos = - cross(sw_vel,Bpos);
Epos = Epos';
Bpos=Bpos';
end

function[vparax,vparay,vparaz,vperpx,vperpy,vperpz] = get_Bcomp(Bpos,vel)
vpara = Bpos(1).*vel(1) + Bpos(2).*vel(2) + Bpos(3).*vel(3);
vparax = vpara.*Bpos(1)./norm(Bpos).^2;
vparay = vpara.*Bpos(2)./norm(Bpos).^2;
vparaz = vpara.*Bpos(3)./norm(Bpos).^2;
vperpx = vel(1) - vparax;
vperpy = vel(2) - vparay;
vperpz = vel(3) - vparaz;
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


% %Plot References
% figure('color','white','position',[300 1600 400 1200]);
% subplot(3,1,1);
% plot(-24*squeeze(B(nx/2,ny/2,:,1))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)));hold on
% plot(-24*squeeze(B(nx/2,ny/2,:,2))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)));
% xline(nz/2-0.5*16)
% title('Vpara speed','fontsize',14);legend({'Vparax';'Vparay'},'fontsize',14,'location','northwest');
% xlim([0 nz/2]);xlabel('Z');ylabel('v/va'); grid on
%
% subplot(3,1,2)
% plot(squeeze(ExB(nx/2,ny/2,:,1))./va);hold on
% plot(squeeze(ExB(nx/2,ny/2,:,2))./va);
% plot(squeeze(ExB(nx/2,ny/2,:,3))./va);
% xline(nz/2-0.5*16)
% title('ExB speed','fontsize',14);legend({'ExBx';'ExBy';'ExBz'},'fontsize',14,'location','northwest');
% xlim([0 nz/2]);xlabel('Z');ylabel('v/va'); grid on
%
% subplot(3,1,3);
% plot(-24*squeeze(B(nx/2,ny/2,:,1))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)) + squeeze(ExB(nx/2,ny/2,:,1))./va);hold on
% plot(-24*squeeze(B(nx/2,ny/2,:,2))./(vecnorm([squeeze(B(nx/2,ny/2,:,1)),squeeze(B(nx/2,ny/2,:,2))],1,2)) + squeeze(ExB(nx/2,ny/2,:,2))./va);
% xline(nz/2-0.5*16)
% title('Sum','fontsize',14);legend({'Sumx';'Sumy'},'fontsize',14,'location','northwest');
% xlim([0 nz/2]);xlabel('Z');ylabel('v/va'); grid on

%Particle Initialization
%     subplot(1,2,1)
%     scatter(vthx,vthz);
%     xlabel('Vx');
%     ylabel('Vz');
%     hold on
% scatterFig = figure;


% [X,Y,Z] = meshgrid(grid_y,grid_x,grid_z);
% quiver3(X,Y,Z,squeeze(B(:,:,:,1)),squeeze(B(:,:,:,2)),squeeze(B(:,:,:,3)))
% xlabel('X'); ylabel('Y'); zlabel('Z')

