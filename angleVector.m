function [RowsofAngles] = angleVector(mainVector,RowsofVectors)
    %angle takes two vectors as an input and calcualtes the dot product and
    %finds the angle in degrees between the two vectors
    
    %output = (acosd(dot(a,b)/(norm(a)*norm(b))));
    RowsofAngles = zeros(length(RowsofVectors),1);
    
    
    %Number of rows in main vector
    mainVectorSize = size(mainVector);
    if mainVectorSize(1) > 1
        
        for i=1:length(RowsofVectors)
            a = mainVector(i,:);
            b = RowsofVectors(i,:);
            RowsofAngles(i) = atan2(norm(cross(a,b)), dot(a,b))*180/pi;
        end
        
        
    else
        
        
        a=mainVector;
        for i=1:length(RowsofVectors)
            b = RowsofVectors(i,:);
            RowsofAngles(i) = atan2(norm(cross(a,b)), dot(a,b))*180/pi;
        end
        
        
    end
    
    
end

