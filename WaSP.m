
function [X_WaSP, C] = WaSP(Y, X, method, wname, lev)

% Created by Ze Jiang on 13/09/2021
% Y: response = m x 1
% X: predictor= (m+l) x n;  m: no. of obs, n: no. of vars, l: no. of forecasts (optional)
% method: discrete wavelet transform, including dwtmar, modwt, modwtmar, and at
% wname: wavelet filter
% lev: decomposition levels

% OUTPUT:
% X_WaSP: variance transformed X = (m+l) x n
% C: covariance vector

% USAGE:
% [X_WaSP, C] = WaSP(Y, X, method, wname, lev)

% REFERENCE:
% Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3). https://doi.org/10.1029/2019WR026962
% Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020). A wavelet-based tool to modulate variance in predictors: An application to predicting drought anomalies. Environmental Modelling & Software, 135, 104907. https://doi.org/10.1016/j.envsoft.2020.104907
% Jiang, Z., Sharma, A., & Johnson, F. (2021). Variable transformations in the spectral domain - Implications for hydrologic forecasting. Journal of Hydrology, 603, 126816. https://doi.org/10.1016/j.jhydrol.2021.126816


% get the dimension
[num_obs, num_var] = size(X) ;

% output matrix
C = nan(lev+1, num_var); 
X_WaSP = nan(num_obs,num_var) ;

% variance transformation for each variable
for i_var = 1 : num_var
    % center
    X_tmp = X(:,i_var)-mean(X(:,i_var)) ; 

    % wavelet transform
    if isequal(method,'modwt')
        X_WT = (modwt(X_tmp, wname, lev))'; 
    elseif isequal(method,'dwtmra')
        X_WT = dwtmra(X_tmp, wname, lev); 
    elseif isequal(method,'modwtmra')
        X_WT_tmp = modwt(X_tmp, wname, lev); 
        X_WT = (modwtmra(X_WT_tmp, wname))'; 
    elseif isequal(method,'at')
        X_WT = AT(X_tmp, wname, lev); 
    else
        disp('This type of wavelet transform does not apply!')
    end
    
%     disp(['Additive:' num2str(sum(abs(sum(X_WT,2)-X_tmp)))])
%     disp(['Variance:' num2str(sum(var(X_WT))-var(X_tmp))])
    
    %corr1 = corrcoef([Y X_WT(1:length(Y),:)]) ; %correct
    %corr1 = cov([Y X_WT(1:length(Y),:)]) ; %wrong
    %C(:,i_var) = corr1(1,2:lev+2); 
    
    % covariance - Eq. 10 in WRR2020 paper
    corr = 1/(length(Y)-1)*Y'*normalize(X_WT(1:length(Y),:)); 
    %disp(normalize(corr1(1,2:lev+2),'norm') - normalize(corr,'norm'))
    C(:,i_var) = corr; 
    
    % normalization to get unit norm
    C_norm = normalize(C(:,i_var),'norm') ;
%     disp(norm(C_norm))
%     disp(var(normalize(X_WT)))
    
    % variance transformation -  Eq. 9 in WRR2020 paper
    X_WaSP(:,i_var) = normalize(X_WT)*(std(X(:,i_var)).*C_norm(:)) ; 
    
    % add mean back
    X_WaSP(:,i_var) = X_WaSP(:,i_var) + mean(X(:,i_var)); 
    
end

