function [min_angle,MVA_vector_resultant] = ambig_angle(MVA_vector,other_vector)
    %angle takes two vectors as an input and calcualtes the dot product and
    %finds the angle in degrees between the two vectors
    %the first vector is from MVAB and ambiguous, thus will flip the vector's
    %sign for all three components and find the lowest angle between the two vectors
    
    %     MVA_ambig(1)=angle(MVA_vector.*[-1,1,1],other_vector);
    %     MVA_ambig(2)=angle(MVA_vector.*[1,-1,1],other_vector);
    %     MVA_ambig(3)=angle(MVA_vector.*[1,1,-1],other_vector);
    %     MVA_ambig(4)=angle(MVA_vector,other_vector);
    
    
    MVA_ambig(1)=angle(MVA_vector.*[-1,-1,-1],other_vector);
    MVA_ambig(2)=angle(MVA_vector,other_vector);
    
    [min_angle,resultant_index] = min(MVA_ambig); %Store the minimum angle and the index for which MVA normal, original or reflected.
    
    
    %pick out the MVA_vector that produces the smallest angle
    if resultant_index == 1
        
        MVA_vector_resultant = MVA_vector.*[-1,-1,-1];
        
    elseif resultant_index == 2
        
        MVA_vector_resultant = MVA_vector;
        
    end
    
    
  
end

