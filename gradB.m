
function [grad_B] = gradB(time,Summ_probes,probes,bx,by,bz,rx,ry,rz,R_inv)
    
    %dbx/dx
    brx = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(bx(time,probes).*rx(time,probes));
    %         bry=sum(bx(time,probes).*ry(time,probes));
    %         brz=sum(bx(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBxx = br*R_inv(:,1,time)/4^2;
    
    %dby/dx
    brx = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(by(time,probes).*rx(time,probes));
    %         bry=sum(by(time,probes).*ry(time,probes));
    %         brz=sum(by(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dByx = br*R_inv(:,1,time)/4^2;
    
    %dbz/dx
    brx = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(bz(time,probes).*rx(time,probes));
    %         bry=sum(bz(time,probes).*ry(time,probes));
    %         brz=sum(bz(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBzx = br*R_inv(:,1,time)/4^2;
    
    
    %dbx/dy
    brx = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    brx=sum(bx(time,probes).*rx(time,probes));
    bry=sum(bx(time,probes).*ry(time,probes));
    brz=sum(bx(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBxy = br*R_inv(:,2,time)/4^2;
    
    %dby/dy
    brx = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    
    %         brx=sum(by(time,probes).*rx(time,probes));
    %         bry=sum(by(time,probes).*ry(time,probes));
    %         brz=sum(by(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dByy = br*R_inv(:,2,time)/4^2;
    
    %dbz/dy
    brx = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(bz(time,probes).*rx(time,probes));
    %         bry=sum(bz(time,probes).*ry(time,probes));
    %         brz=sum(bz(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBzy = br*R_inv(:,2,time)/4^2;
    
    
    %dbx/dz
    brx = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(bx(time,probes).*rx(time,probes));
    %         bry=sum(bx(time,probes).*ry(time,probes));
    %         brz=sum(bx(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBxz = br*R_inv(:,3,time)/4^2;
    
    %dby/dz
    brx = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(by(time,probes).*rx(time,probes));
    %         bry=sum(by(time,probes).*ry(time,probes));
    %         brz=sum(by(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dByz = br*R_inv(:,3,time)/4^2;
    
    %dbz/dz
    brx = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    %         brx=sum(bz(time,probes).*rx(time,probes));
    %         bry=sum(bz(time,probes).*ry(time,probes));
    %         brz=sum(bz(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBzz = br*R_inv(:,3,time)/4^2;
    
    %form the gradient matrix
    grad_B = [dBxx, dByx, dBzx;...
        dBxy, dByy, dBzy;...
        dBxz, dByz, dBzz];
    
    %Constraint for delB = 0 ;
    %         lambda = trace(grad_B)/(trace(R_inv(:,:,time)));
    % %
    %         grad_B = grad_B - lambda*R_inv(:,:,time);
end
