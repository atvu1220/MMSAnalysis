function [output] = angle(a,b)
    %angle takes two vectors as an input and calcualtes the dot product and
    %finds the angle in degrees between the two vectors
    
    %output = (acosd(dot(a,b)/(norm(a)*norm(b))));
    
    
    output = atan2(norm(cross(a,b)), dot(a,b))*180/pi;
end

