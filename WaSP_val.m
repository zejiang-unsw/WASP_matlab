
function [X_WaSP] = WaSP_val(X, C, method, wname)

% Created by Ze Jiang on 18/09/2021 (ze.jiang@unsw.edu.au)
% Note: This is used for applying derived C in the validation. 

% dimension of predictor
[num_obs, num_var] = size(X) ;
lev = size(C,1)-1;

% output matrix
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

    %standardization 
    X_WT_c=X_WT-repmat(mean(X_WT),num_obs,1);
    X_WT_norm=X_WT_c./repmat(std(X_WT,0,1),num_obs,1);
    
    % normalization to get unit norm
    C_norm = C(:,i_var)./sqrt(sum(C(:,i_var).^2)) ;
    
    % variance transformation -  Eq. 9 in WRR2020 paper
    X_WaSP(:,i_var) = X_WT_norm*(std(X(:,i_var)).*C_norm) ; 
    
    % add mean back
    X_WaSP(:,i_var) = X_WaSP(:,i_var) + mean(X(:,i_var)); 
	
	% maintain the original trend of the variable
	[rho, pval] = corr(X_WaSP(:,i_var),X(:,i_var)); 
    if rho < 0 && pval < 0.05
        X_WaSP(:,i_var) = -X_WaSP(:,i_var); 
    end
end
end