function [B_eigenvectors,lambda,RotationMatrix] = mvab(B,min_flag)

   
    mean_B = mean(B(:,1:3));

    m1 = mean(B(:,[1 2 3 1 1 2]).*B(:,[1 2 3 2 3 3]));
    m2 = mean_B([1 2 3 1 1 2]).*mean_B([1 2 3 2 3 3]);

    m3 = m1-m2;
    M = [[m3(1) m3(4) m3(5)]' [m3(4) m3(2) m3(6)]' [m3(5) m3(6) m3(3)]'];
    
    %Now we must calculate eigenvalues and vectors
    [V,L] = eig(M);
    [lambda,I] = sort(diag(L),'descend');
    RotationMatrix = V(:, I);
    RotationMatrix = RotationMatrix.*[1,1,1];
    B_eigenvectors = B*(RotationMatrix); %maybe need to transpose 4/10/2020

    
end

