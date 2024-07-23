
function [X_WaSP] = WaSP_val(X, C, method, wname, flag_sign)

% Created by Ze Jiang on 18/09/2021 (ze.jiang@unsw.edu.au)
% Note: This is used for applying derived C in the validation. 

if ~exist('flag_sign','var'), flag_sign=0; end

% dimension of predictor
[num_obs, num_var] = size(X) ;
lev = size(C,1)-1;

% output matrix
X_WaSP = nan(num_obs,num_var) ;

% variance transformation for each variable
for i_var = 1 : num_var
    % center
    X_tmp = X(:,i_var)-mean(X(:,i_var)) ; 

    % increase size to fit the dimension of C
    n_rep = round(2^(lev+1)/size(X,1));
    X_tmp = [X_tmp; repmat(X_tmp, [n_rep 1])]; 

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

    X_WT = X_WT(1:num_obs,:);

    %standardization 
    X_WT_c=X_WT-repmat(mean(X_WT),num_obs,1);
    X_WT_norm=X_WT_c./repmat(std(X_WT,0,1),num_obs,1);
    
    % normalization to get unit norm
    C_norm = C(:,i_var)./sqrt(sum(C(:,i_var).^2));
    
    % variance transformation -  Eq. 9 in WRR2020 paper
    X_WaSP(:,i_var) = X_WT_norm*(std(X(:,i_var)).*C_norm) ; 
    
	% maintain the original trend of the variable
	if flag_sign
		[rho, pval] = corr(X_WaSP(:,i_var),X(:,i_var)); 
		if rho < 0 && pval < 0.05 % can change to other confidence level
			C_norm = -C_norm ; 
			X_WaSP(:,i_var) = X_WT_norm*(std(X(:,i_var)).*C_norm) ; 
		end
	end
	
    % add mean back
    X_WaSP(:,i_var) = X_WaSP(:,i_var) + mean(X(:,i_var)); 
	
end
end