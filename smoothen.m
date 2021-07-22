function [input] = smoothen(input,n)
    
    
    
    
    
    if n > 0
        for i=1:n
            mm = mean(input);
            ii = find(abs(input - mm) > 1.2*std(input));
            input_tmp = input;
            input_tmp(ii) = mm;
            input_tmp = smooth(input_tmp);
            input(ii) = input_tmp(ii);
            %plot(sig,irf_lstyle(i));
        end
    end
    
    input = smooth(input);
    
    
    
    
end

