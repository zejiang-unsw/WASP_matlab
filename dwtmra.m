function [X_DWT_MRA] = dwtmra(X, wname, lev)

% Created by Ze Jiang on 15/09/2021
% DWT MRA wavelet transform

    N = length(X); 
    
    %1-D wavelet decomposition
    [C,L] = wavedec(X,lev,wname); 

    X_DWT_MRA = nan(N,lev+1); 
    for i=1:lev
      X_DWT_MRA(:,i) = wrcoef('d',C,L,wname,i);
    end
    X_DWT_MRA(:,lev+1) = wrcoef('a',C,L,wname,lev);
      
%     disp(['Additive:' num2str(sum(abs(sum(X_DWT_MRA,2)-X)))])
%     disp(['Variance:' num2str(sum(var(X_DWT_MRA))-var(X))])
end 