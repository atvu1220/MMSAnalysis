function [r] = rotate_GPE2GSE(v,rotated_r)
%     v=[-300,50,-20]
%     
%     r = [ 66000, 500, 100]
    
    %Chu 2017 rotate from GPE back to GSE using input solar wind vector and position of HFA
    %Negated Rotation angles from GSE 2 GPE...hope it works
    v = reshape(v,1,3);
    rotated_r = reshape(rotated_r,3,1);
    theta = acos( -v(1) / sqrt( v(1)^2 + v(2)^2 ) );
    if -v(2) < 0
        theta = -abs(theta);
    else
        theta = abs(theta);
    end
    theta_rotationMatrix = [ cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    
    
    u = [v(2) / sqrt(v(1)^2 + v(2)^2), -v(1) / sqrt(v(1)^2 + v(2)^2),0 ];
    phi = acos( (v(1)^2 + v(2)^2) / ( sqrt(v(1)^2 + v(2)^2 ) * sqrt(v(1)^2 + v(2)^2 + v(3)^2 ) ) );
    
    if -v(3) < 0
        phi = abs(phi);
    else
        phi = -abs(phi);
    end
    
    R = [cos(phi) + u(1)^2*(1-cos(phi)) u(1)*u(2)*(1-cos(phi)) u(2)*sin(phi)  ;...
        u(1)*u(2)*(1-cos(phi)) cos(phi) + u(2)^2*(1-cos(phi)) -u(1)*sin(phi) ;...
        -u(2)*sin(phi) u(1)*sin(phi) cos(phi)] ;
    
    
    rotated_r = inv(R)*rotated_r;
    r = inv(theta_rotationMatrix)*rotated_r;
end

