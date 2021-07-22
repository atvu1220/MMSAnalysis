function [vector] = Unfluctuate(vector)
    %Removes those pesky spikes

%     
    average = mode(sign(vector));
    for i=1:length(vector)
        if sign(vector(i)) ~= sign(average)
            vector(i) = (-1)*vector(i);
        end
    end

 %average = mean(vector);
% for i=2:length(vector)
%     if sign(vector(i)) ~= sign(vector(i-1)) && abs(abs(vector(i)) - abs(vector(i-1))) < 0.5
%         vector(i) = -1*vector(i);
%         disp('yes')
%     end
% end

    
end