function out = stepwise_VT(Ycal, Xcal, Xval, opts)
%STEPWISE_VT Stepwise high-order WASP variance transformation in MATLAB.
%
%   out = STEPWISE_VT(Ycal, Xcal, Xval, opts)

%
%   PURPOSE
%   -------
%   The normal MATLAB WaSP.m transforms all predictors independently. 
%       1. transform each remaining candidate predictor using WaSP;
%       2. calculate partial informational correlation (PIC);
%       3. select the predictor with the largest PIC;
%       4. condition the next step on already selected transformed predictors;
%       5. stop when R2 decreases, PIC is below the significance threshold,
%          no candidates remain, or nvarmax is reached;
%       6. apply the selected calibration transformations to Xval using
%          WaSP_val.m.
%
%   REQUIRED FILES IN THE SAME MATLAB PATH
%   --------------------------------------
%   Existing WASP_matlab functions:
%       WaSP.m
%       WaSP_val.m
%       dwtmra.m       required for method='dwtmra'
%       AT.m           required for method='at'
%
%   MATLAB toolbox:
%       Wavelet Toolbox.
%

%
%   INPUTS
%   ------
%   Ycal : nCal x 1 response vector.
%   Xcal : nCal x p predictor matrix for calibration.
%   Xval : nVal x p predictor matrix for validation. Use [] if not needed.
%   opts : structure with fields below. Missing fields are filled by defaults.
%
%       opts.alpha    = 0.10;       % R default significance level
%       opts.nvarmax  = 4;          % max selected predictors
%       opts.mode     = 'MRA';      % 'MRA', 'MODWT', or 'AT'
%       opts.method   = 'dwt';      % 'dwt', 'dwtmra', 'modwt', 'modwtmra', 'at'
%       opts.wf       = 'db4';      % MATLAB style 'db4' or R style 'd8'
%       opts.J        = 4;          % decomposition level
%       opts.Kres     = 5;         % K for KNN residualisation
%       opts.boundary = 'periodic'; % stored for compatibility
%       opts.flagSign = false;      % passed to WaSP/WaSP_val as flag_sign
%       opts.verbose  = true;
%
%   COMMON USE INSIDE YOUR CURRENT MAIN SCRIPT
%   ------------------------------------------
%       optsDWT.alpha    = 0.10;
%       optsDWT.nvarmax  = nvarmaxStep;
%       optsDWT.mode     = 'MRA';
%       optsDWT.method   = 'dwt';
%       optsDWT.wf       = wname;       % e.g. 'db4'
%       optsDWT.J        = lev_default; % e.g. 4
%       optsDWT.Kres     = K;
%       optsDWT.verbose  = true;
%
%       fitDWT = stepwise_VT(Yc, Xc, Xv, optsDWT);
%
%   OUTPUT
%   ------
%   out.cpy      : selected original predictor column indices.
%   out.cpyPIC   : PIC values of selected predictors.
%   out.wt       : partial weights for selected transformed predictors.
%   out.r2       : cumulative R2 after accepted selections.
%   out.dp       : selected raw calibration predictors after sequential
%                  conditioning.
%   out.dp_n     : transformed selected calibration predictors.
%   out.dp_val   : selected raw validation predictors after sequential
%                  conditioning.
%   out.dp_n_val : transformed selected validation predictors.
%   out.Zcal     : same as out.dp_n, for compatibility with user scripts.
%   out.Zval     : same as out.dp_n_val, for compatibility with user scripts.
%   out.S        : selected WaSP covariance vectors, one column per selection.
%   out.stage    : detailed information for each selected step.
%




    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    Ycal = double(Ycal(:));
    Xcal = double(Xcal);

    if nargin < 3 || isempty(Xval)
        Xval = zeros(0, size(Xcal,2));
    else
        Xval = double(Xval);
    end

    if size(Xcal,1) ~= numel(Ycal)
        error('stepwise_VT:InputSize', 'Ycal length must equal size(Xcal,1).');
    end
    if size(Xcal,2) ~= size(Xval,2)
        error('stepwise_VT:InputSize', 'Xcal and Xval must have the same number of columns.');
    end
    if any(~isfinite(Ycal)) || any(~isfinite(Xcal(:))) || any(~isfinite(Xval(:)))
        error('stepwise_VT:NonFiniteInput', 'Inputs contain NaN or Inf. Remove missing rows before calling stepwise_VT.');
    end

    opts = local_defaults(opts, numel(Ycal));
    local_check_required_files(opts);

    nCal = numel(Ycal);
    nPred = size(Xcal,2);
    nVal = size(Xval,1);

    selected = [];
    selectedPIC = [];
    r2 = [];
    Zcal = zeros(nCal,0);
    S = zeros(opts.J+1,0);
    stage = struct([]);
    remaining = 1:nPred;

    if opts.verbose
        fprintf('\nstepwise_VT started | n=%d, p=%d, method=%s, wavelet=%s, J=%d\n', ...
            nCal, nPred, opts.methodWaSP, opts.wnameMatlab, opts.J);
    end

    step = 0;
    while ~isempty(remaining) && numel(selected) < opts.nvarmax
        step = step + 1;

        % Residualise response and remaining candidate predictors against
        % already-selected transformed predictors. This is the partial part
        % of the PIC calculation.
        if isempty(Zcal)
            yIn = Ycal;
            xCand = Xcal(:, remaining);
        else
            yIn = Ycal - local_knn_loocv(Ycal, Zcal, opts.Kres);
            xCand = zeros(nCal, numel(remaining));
            for jj = 1:numel(remaining)
                x0 = Xcal(:, remaining(jj));
                xCand(:,jj) = x0 - local_knn_loocv(x0, Zcal, opts.Kres);
            end
        end

        % Use the existing MATLAB WaSP.m for the actual wavelet variance
        % transformation of all remaining candidate predictors.
        [xCandVT, Ccand] = WaSP(yIn, xCand, opts.methodWaSP, opts.wnameMatlab, opts.J, opts.flagSign);

        % PIC calculation, following the R idea:
        % PMI -> PIC = sqrt(1 - exp(-2*PMI)).
        picVals = zeros(1, numel(remaining));
        for jj = 1:numel(remaining)
            pmi = local_pmi(yIn, xCandVT(:,jj));
            if ~isfinite(pmi) || pmi < 0
                pmi = 0;
            end
            picVals(jj) = sqrt(max(0, 1 - exp(-2*pmi)));
        end

        [picMax, bestLocal] = max(picVals);
        bestCol = remaining(bestLocal);
        bestVT = xCandVT(:, bestLocal);
        bestC = Ccand(:, bestLocal);

        Ztrial = [Zcal, bestVT]; %#ok<AGROW>

        yhat_cv = local_knn_loocv(Ycal, Ztrial, opts.Kres);
        denom = sum((Ycal - mean(Ycal)).^2);
        if denom <= eps
            r2New = NaN;
        else
            r2New = 1 - sum((Ycal - yhat_cv).^2) / denom;
        end

        df = max(nCal - size(Ztrial,2), 1);
        tcrit = local_tinv(1 - opts.alpha, df);
        picThreshold = sqrt(tcrit^2 / (tcrit^2 + df));

        reject = false;
        reason = '';

        % Same stopping rule as R stepwise.VT: first variable is accepted;
        % from step 2 onward, reject if R2 decreases or PIC is too small.
        if step > 1
            if isfinite(r2New) && isfinite(r2(end)) && r2New < r2(end)
                reject = true;
                reason = 'R2 decreased';
            end
            if picMax < picThreshold
                reject = true;
                if isempty(reason)
                    reason = 'PIC below threshold';
                else
                    reason = [reason, ' and PIC below threshold']; %#ok<AGROW>
                end
            end
        elseif opts.firstMustPassPIC && picMax < picThreshold
            reject = true;
            reason = 'first PIC below threshold';
        end

        if reject
            if opts.verbose
                fprintf('Step %d rejected | predictor=%d | PIC=%.5g | threshold=%.5g | R2=%.5g | %s\n', ...
                    step, bestCol, picMax, picThreshold, r2New, reason);
            end
            break
        end

        selected = [selected, bestCol]; %#ok<AGROW>
        selectedPIC = [selectedPIC; picMax]; %#ok<AGROW>
        r2 = [r2; r2New]; %#ok<AGROW>
        Zcal = Ztrial;
        S = [S, bestC]; %#ok<AGROW>

        st = struct();
        st.step = step;
        st.originalColumn = bestCol;
        st.cpyPIC = picMax;
        st.C = bestC(:);
        st.picThreshold = picThreshold;
        st.r2 = r2New;
        st.method = opts.methodWaSP;
        st.wavelet = opts.wnameMatlab;
        st.J = opts.J;
        stage = [stage; st]; %#ok<AGROW>

        if opts.verbose
            fprintf('Step %d accepted | predictor=%d | PIC=%.5g | threshold=%.5g | R2=%.5g\n', ...
                step, bestCol, picMax, picThreshold, r2New);
        end

        remaining(bestLocal) = [];
    end

    % Calibration raw selected predictors after sequential conditioning.
    dpCal = local_condition_selected_raw(Xcal, selected, Zcal, opts.Kres);

    % Validation transformation. The selected covariance vectors are applied
    % through WaSP_val.m. Later selected validation predictors are residualised
    % against already transformed validation predictors, matching the logic of
    % R stepwise.VT.val.
    if isempty(selected)
        Zval = zeros(nVal,0);
        dpVal = zeros(nVal,0);
        wt = [];
    else
        [Zval, dpVal] = local_validation_transform(Xval, selected, S, opts);
        wt = local_pw_calc(Ycal, Zcal, selectedPIC, opts.Kres);
    end

    out = struct();
    out.cpy      = selected(:)';
    out.cpyPIC   = selectedPIC(:);
    out.wt       = wt(:);
    out.lstwet   = local_lsq_weights(Ycal, Zcal);
    out.x        = Ycal;
    out.py       = Xcal;
    out.r2       = r2(:);
    out.dp       = dpCal;
    out.dp_n     = Zcal;
    out.dp_n_cal = Zcal;
    out.dp_val   = dpVal;
    out.dp_n_val = Zval;
    out.Zcal     = Zcal;
    out.Zval     = Zval;
    out.S        = S;
    out.wavelet  = opts.wnameMatlab;
    out.method   = opts.methodWaSP;
    out.pad      = opts.pad;
    out.boundary = opts.boundary;
    out.stage    = stage;
    out.opts     = opts;

    if opts.verbose
        fprintf('stepwise_VT finished | selected predictors: %s\n\n', mat2str(out.cpy));
    end
end

%% ========================================================================
% Validation and conditioning helpers
% ========================================================================

function [Zval, dpVal] = local_validation_transform(Xval, selected, S, opts)
    nVal = size(Xval,1);
    nSel = numel(selected);
    Zval = zeros(nVal, nSel);
    dpVal = zeros(nVal, nSel);

    if nVal == 0 || nSel == 0
        return
    end

    for ii = 1:nSel
        xCurrent = Xval(:, selected(ii));

        if ii > 1
            xCurrent = xCurrent - local_knn_predict_same(xCurrent, Zval(:,1:ii-1), opts.Kres);
        end

        dpVal(:,ii) = xCurrent;
        Zval(:,ii) = WaSP_val(xCurrent, S(:,ii), opts.methodWaSP, opts.wnameMatlab, opts.flagSign);
    end
end

function dpCal = local_condition_selected_raw(Xcal, selected, Zcal, K)
    n = size(Xcal,1);
    nSel = numel(selected);
    dpCal = zeros(n, nSel);

    if nSel == 0
        return
    end

    dpCal(:,1) = Xcal(:, selected(1));
    for ii = 2:nSel
        xCurrent = Xcal(:, selected(ii));
        xCurrent = xCurrent - local_knn_loocv(xCurrent, Zcal(:,1:ii-1), K);
        dpCal(:,ii) = xCurrent;
    end
end

%% ========================================================================
% PIC / PMI helpers
% ========================================================================

function pmi = local_pmi(X, Y)
    X = double(X(:));
    Y = double(Y(:));
    keep = isfinite(X) & isfinite(Y);
    X = X(keep);
    Y = Y(keep);

    n = numel(X);
    if n < 5 || std(X) <= eps || std(Y) <= eps
        pmi = 0;
        return
    end

    pdfX = local_kde_1d(X);
    pdfY = local_kde_1d(Y);
    pdfXY = local_kde_2d([X Y]);

    ratio = pdfXY ./ max(pdfX .* pdfY, realmin);
    vals = log(max(ratio, realmin));
    vals = vals(isfinite(vals));

    if isempty(vals)
        pmi = 0;
    else
        pmi = mean(vals);
    end
end

function dens = local_kde_1d(z)
    z = double(z(:));
    n = numel(z);

    bw = 1.5 * local_bw_nrd0(z);
    if ~isfinite(bw) || bw <= eps
        bw = max(std(z), eps);
    end

    D2 = (z - z').^2;
    K = exp(-D2 ./ (2*bw^2));
    dens = sum(K, 2) ./ (sqrt(2*pi) * bw * n);
    dens = max(dens, realmin);
end

function dens = local_kde_2d(Z)
    Z = double(Z);
    n = size(Z,1);
    d = size(Z,2);

    if n < 5
        dens = realmin * ones(n,1);
        return
    end

    C = cov(Z);
    if any(~isfinite(C(:)))
        C = eye(d);
    end

    reg = 1e-10 * trace(C) / max(d,1);
    if ~isfinite(reg) || reg <= 0
        reg = 1e-10;
    end
    C = C + reg * eye(d);

    detC = det(C);
    if ~isfinite(detC) || detC <= eps
        C = C + 1e-6 * eye(d);
        detC = det(C);
    end

    sigma = 1.5 * (4/(d+2))^(1/(d+4)) * n^(-1/(d+4));
    constant = (sqrt(2*pi) * sigma)^d * sqrt(max(detC, realmin)) * n;

    dens = zeros(n,1);
    for ii = 1:n
        D = bsxfun(@minus, Z, Z(ii,:));
        q = sum((D / C) .* D, 2);
        dens(ii) = sum(exp(-q ./ (2*sigma^2))) ./ constant;
    end
    dens = max(dens, realmin);
end

function bw = local_bw_nrd0(z)
    z = sort(double(z(:)));
    n = numel(z);
    if n < 2
        bw = 1;
        return
    end

    s = std(z,0,1);
    q25 = local_quantile_sorted(z, 0.25);
    q75 = local_quantile_sorted(z, 0.75);
    iqrVal = q75 - q25;

    lo = min(s, iqrVal/1.34);
    if ~isfinite(lo) || lo <= eps
        lo = s;
    end
    if ~isfinite(lo) || lo <= eps
        lo = 1;
    end

    bw = 0.9 * lo * n^(-1/5);
end

function q = local_quantile_sorted(z, p)
    n = numel(z);
    if n == 1
        q = z;
        return
    end
    pos = 1 + (n-1)*p;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = z(lo);
    else
        q = z(lo) + (pos-lo) * (z(hi)-z(lo));
    end
end

%% ========================================================================
% KNN helpers
% ========================================================================

function yhat = local_knn_loocv(y, Z, K)
    y = double(y);
    Z = double(Z);
    n = size(Z,1);

    if isempty(Z) || size(Z,2) == 0
        yhat = repmat(mean(y,1), size(y,1), 1);
        return
    end
    if size(y,1) ~= n
        error('stepwise_VT:KNNSize', 'KNN input length mismatch.');
    end

    K = max(1, min(round(K), max(n-1,1)));
    yhat = nan(size(y));

    for ii = 1:n
        d = sum(bsxfun(@minus, Z, Z(ii,:)).^2, 2);
        d(ii) = inf;
        [~, ord] = sort(d, 'ascend');
        idx = ord(1:K);
        yhat(ii,:) = mean(y(idx,:), 1);
    end
end

function yhat = local_knn_predict_same(y, Z, K)
    % Non-LOOCV version for validation-period sequential conditioning.
    % Since y and Z are from the same validation block and there is no separate
    % training set, this follows the same practical residualisation structure
    % used by the single-block R validation code.
    y = double(y);
    Z = double(Z);
    n = size(Z,1);

    if isempty(Z) || size(Z,2) == 0
        yhat = repmat(mean(y,1), size(y,1), 1);
        return
    end

    K = max(1, min(round(K), max(n-1,1)));
    yhat = nan(size(y));

    for ii = 1:n
        d = sum(bsxfun(@minus, Z, Z(ii,:)).^2, 2);
        d(ii) = inf;
        [~, ord] = sort(d, 'ascend');
        idx = ord(1:K);
        yhat(ii,:) = mean(y(idx,:), 1);
    end
end

%% ========================================================================
% Weight helpers
% ========================================================================

function wt = local_pw_calc(x, Z, cpyPIC, K)
    Z = double(Z);
    x = double(x(:));
    cpyPIC = double(cpyPIC(:));

    nCol = size(Z,2);
    wt = zeros(nCol,1);

    if nCol == 0
        wt = [];
    elseif nCol == 1
        wt(1) = local_scale_std_ratio(x, Z(:,1), [], K) * cpyPIC(1);
    else
        for ii = 1:nCol
            others = setdiff(1:nCol, ii);
            wt(ii) = local_scale_std_ratio(x, Z(:,ii), Z(:,others), K) * cpyPIC(ii);
        end
    end
end

function ratio = local_scale_std_ratio(x, zin, zout, K)
    x = double(x(:));
    zin = double(zin(:));

    if nargin < 3 || isempty(zout)
        ratio = 1;
        return
    end

    zout = double(zout);

    xhat = local_knn_loocv(x, zout, K);
    vx = var(x);
    if vx <= eps || ~isfinite(vx)
        r1 = 1;
    else
        r1 = sqrt(max(var(x - xhat),0) / vx);
    end

    zhat = local_knn_loocv(zin, zout, K);
    vz = var(zin);
    if vz <= eps || ~isfinite(vz)
        r2 = 1;
    else
        r2 = sqrt(max(var(zin - zhat),0) / vz);
    end

    ratio = 0.5 * (r1 + r2);
end

function w = local_lsq_weights(y, Z)
    if isempty(Z) || size(Z,2) == 0
        w = [];
        return
    end
    try
        b = [ones(size(Z,1),1), Z] \ y(:);
        w = abs(b(2:end));
    catch
        w = ones(size(Z,2),1);
    end
end

%% ========================================================================
% Option and compatibility helpers
% ========================================================================

function opts = local_defaults(opts, n)
    opts = local_set_default(opts, 'alpha', 0.10);
    opts = local_set_default(opts, 'nvarmax', 4);
    opts = local_set_default(opts, 'mode', 'MRA');
    opts = local_set_default(opts, 'method', 'dwt');
    opts = local_set_default(opts, 'wf', 'db4');
    opts = local_set_default(opts, 'wname', []);
    opts = local_set_default(opts, 'J', []);
    opts = local_set_default(opts, 'lev', []);
    opts = local_set_default(opts, 'Kres', []);
    opts = local_set_default(opts, 'K', []);
    opts = local_set_default(opts, 'boundary', 'periodic');
    opts = local_set_default(opts, 'pad', 'zero');
    opts = local_set_default(opts, 'flagSign', false);
    opts = local_set_default(opts, 'verbose', true);
    opts = local_set_default(opts, 'firstMustPassPIC', false);

    if isempty(opts.wname)
        opts.wnameMatlab = local_wavelet_name_to_matlab(opts.wf);
    else
        opts.wnameMatlab = local_wavelet_name_to_matlab(opts.wname);
    end

    if isempty(opts.J)
        if ~isempty(opts.lev)
            opts.J = opts.lev;
        else
            opts.J = max(1, floor(log2(n)) - 1);
        end
    end
    opts.J = max(1, round(opts.J));

    if isempty(opts.Kres)
        if ~isempty(opts.K)
            opts.Kres = opts.K;
        else
            opts.Kres = max(1, round(sqrt(n/2)));
        end
    end

    modeLower = lower(char(opts.mode));
    methodLower = lower(char(opts.method));

    if strcmp(modeLower, 'mra')
        if strcmp(methodLower, 'dwt') || strcmp(methodLower, 'dwtmra')
            opts.methodWaSP = 'dwtmra';
        elseif strcmp(methodLower, 'modwt') || strcmp(methodLower, 'modwtmra')
            opts.methodWaSP = 'modwtmra';
        else
            opts.methodWaSP = methodLower;
        end
    elseif strcmp(modeLower, 'modwt')
        opts.methodWaSP = 'modwt';
    elseif strcmp(modeLower, 'at') || strcmp(modeLower, 'a trous') || strcmp(modeLower, 'atrous')
        opts.methodWaSP = 'at';
    else
        opts.methodWaSP = methodLower;
    end
end

function opts = local_set_default(opts, fieldName, value)
    if ~isfield(opts, fieldName) || isempty(opts.(fieldName))
        opts.(fieldName) = value;
    end
end

function wname = local_wavelet_name_to_matlab(wf)
    wf = strtrim(char(wf));

    if strcmpi(wf, 'haar')
        wname = 'haar';
        return
    end

    % R/waveslim style: d8 means Daubechies length 8 = MATLAB db4.
    tok = regexp(lower(wf), '^d(\d+)$', 'tokens', 'once');
    if ~isempty(tok)
        L = str2double(tok{1});
        if isfinite(L) && mod(L,2) == 0
            wname = sprintf('db%d', L/2);
            return
        end
    end

    % MATLAB style, e.g., db4, db8, sym4, coif2.
    wname = wf;
end

function local_check_required_files(opts)
    if exist('WaSP', 'file') ~= 2
        error('stepwise_VT:MissingWaSP', 'Cannot find WaSP.m. Put stepwise_VT.m in the WASP_matlab folder or add that folder to the MATLAB path.');
    end
    if exist('WaSP_val', 'file') ~= 2
        error('stepwise_VT:MissingWaSPVal', 'Cannot find WaSP_val.m. Put stepwise_VT.m in the WASP_matlab folder or add that folder to the MATLAB path.');
    end
    if strcmp(opts.methodWaSP, 'dwtmra') && exist('dwtmra', 'file') ~= 2
        error('stepwise_VT:MissingDWTMRA', 'method=dwtmra requires dwtmra.m from the WASP_matlab repository.');
    end
    if strcmp(opts.methodWaSP, 'at') && exist('AT', 'file') ~= 2
        error('stepwise_VT:MissingAT', 'method=at requires AT.m from the WASP_matlab repository.');
    end
end

function tq = local_tinv(p, df)
    if exist('tinv', 'file') == 2
        tq = tinv(p, df);
        return
    end

    % Fallback normal approximation for installations without the Statistics
    % and Machine Learning Toolbox. This is accurate enough for large samples,
    % but tinv is preferred when available.
    tq = sqrt(2) * erfinv(2*p - 1);
    if df > 2 && isfinite(tq)
        % Simple first-order expansion from normal to t quantile.
        tq = tq + (tq^3 + tq) / (4*df);
    end
end
