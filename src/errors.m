% errors()
% This subroutine does error propagation on the computation of carbonate system variables 
% from errors (or uncertainties) on six input 
%  - pair of carbonate system variables 
%  - nutrients (silicate and phosphate concentrations)
%  - temperature and salinity
% plus errors on dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
%
% It computes numerical derivatives then applies error propagation using the method of moments.
% The method of moments is a very general technique for estimating the second moment of a variable z
% (variance or standard deviation) based on a first-order approximation to z.
%
%**************************************************************************
%
%  **** SYNTAX:
%  [std_err,headers]=errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,..  .
%                  SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
%                  ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,r,...
%                  pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS)
% 
%  **** SYNTAX EXAMPLES:
%  [Result]          = errors(2400,2200,1,2,35,0,25,4200,0,15,1,2,2,0.01,0.01,0,0,0,0,1,4,1)
%  [Result,Headers]  = errors(2400,   8,1,3,35,0,25,4200,0,15,1,2,0.001,0,0,0,0,0,0,1,4,1)
%  [Result,Headers]  = errors(500,   8,5,3,35,0,25,4200,0,15,1,2,0.001,0,0,0,0,'',0,1,4,1)
%  [A]               = errors(2400,2000:10:2400,1,2,35,0,25,4200,0,15,2,2,0,0,0,0,'',0,1,1,4,1)
%  [A]               = errors(2400,2200,1,2,0:1:35,0,25,4200,0,15,1,2,2,0,0,0,0,'',0,1,4,1)
%  epK = [0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02];
%  [A]               = errors(2400,2200,1,2,35,0,25,0:100:4200,0,15,1,2,2,0,0,0,0,epK,0,1,4,1)
%  
%**************************************************************************
%
% INPUT:
%
%   - ePAR1, ePAR2   :  standard error (or uncertainty) on PAR1 and PAR2 of input pair of carbonate system variables
%   - eS, eT         :  standard error (or uncertainty) on Salinity and Temperature
%   - ePO4, eSI      :  standard error (or uncertainty) on Phosphate and Silicate total concentrations
%   - epK            :  standard error (or uncertainty) on all seven dissociation constants (a vector)
%   - r              :  correlation coefficient between PAR1 AND PAR2 (typicaly 0)
%   - others         :  same as input of subroutine  CO2sys() : scalar or vectors
%
% All parameters may be scalars or vectors except epK.
%   * epK must be vector of seven values : errors of pK0, pK1, pK2, pKb, pKw, pKspa and pKspc
%     these errors are assumed to be equal for all input data.
%
%     if epK is empty (= ''), default values for epK are taken
%     Default standard errors are :
%        pK0   :  0.004 
%        pK1   :  0.015
%        pK2   :  0.03
%        pKb   :  0.01    boric acid
%        pKw   :  0.01    water dissociation
%        pKspa :  0.02    solubility product of Aragonite 
%        pKspc :  0.02    solubility product of Calcite
%        TB    :  0.01    total boron
%
% In constrast, ePAR1, ePAR2, eS, eT, ePO4 and eSI, 
%   - if vectors, are errors associated with each data point
%   - if scalars, are one error value associated to all data points
% The same for parameter "r".
%
% If 'r' is not 0, having a value between -1.0 and 1.0, it indicates the correlation 
% between uncertainties of the input pair of carbonate system variables.
% By default, 'r' is zero. However, for some pairs the user may want to specify a
% different value. For example, measurements of pCO2 and pH are often anti-correlated.
% The same goes for two other pairs: 'CO2 and CO3' and 'pCO2 and
% CO3'. But even for these cases, care is needed before using non-zero values of 'r'.
% 
% When the user wishes to propagate errors for an individual
% measurement, 'r' should ALWAYS be zero if each member of the input pair is
% measured independently. In this case, we are interested in the
% correlation between the uncertainties in those measurements, not in
% the correlation between the measurments themselves. Uncertainties from
% those measurements are probably not correlated if they come from
% different instruments. Conversely, if users are interested in the
% error in the mean of a distribution of measurements (i.e., if they are
% propagating standard errors instead of standard deviations), one
% should then also account for the correlation between the measurements of
% the two variables of the input pair.
% 
% For input pairs where one member is pH, this 'errors' routine automatically
% inverses the sign of 'r'.
% The reason for that is that the associated derivatives are computed in terms of 
% the hydrogen ion concentration H+, not pH. Therefore for each of these 6
% flags, if the user wants to compute 'r' that should be done (1) using
% the H+ concentration instead of pH, and (2) the sign of that computed 'r'
% should be inversed when passing it as an argument to this routine.
% For these cases, to express perfect anticorrelation with pH, the user should 
% use 'r=1.0'. For all other flags (those without pH as a member), the sign
% of 'r' should not be inversed.
% 
%**************************************************************************
%
% OUTPUT: * an array containing standard error (or uncertainty) for the following variables
%           (one row per sample):
%         *  a cell-array containing crudely formatted headers
%
%    POS  PARAMETER        UNIT
%
%    01 - TAlk                 (umol/kgSW)
%    02 - TCO2                 (umol/kgSW)
%    03 - [H+] input           (umol/kgSW)
%    04 - pCO2 input           (uatm)
%    05 - fCO2 input           (uatm)
%    06 - HCO3 input           (umol/kgSW)
%    07 - CO3 input            (umol/kgSW)
%    08 - CO2 input            (umol/kgSW)
%    09 - OmegaCa input        ()
%    10 - OmegaAr input        ()
%    11 - xCO2 input           (ppm)
%    12 - [H+] output          ()
%    13 - pCO2 output          (uatm)
%    14 - fCO2 output          (uatm)
%    15 - HCO3 output          (umol/kgSW)
%    16 - CO3 output           (umol/kgSW)
%    17 - CO2 output           (umol/kgSW)
%    18 - OmegaCa output       ()
%    19 - OmegaAr output       ()
%    20 - xCO2 output          (ppm)
%
% Remark : if all input pairs are of the same type, standard error of input pairs are omitted
%

function [total_error, headers] = errors (PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
                  ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,r,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);

    global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi;
    
    % Input conditioning
    % ------------------

    % Determine lengths of input vectors
    veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
                length(PAR2TYPE) length(SAL) length(TEMPIN)...
                length(TEMPOUT) length(PRESIN) length(PRESOUT)...
                length(SI) length(PO4) length(ePAR1) length(ePAR2)...
                length(eSAL) length(eTEMP)...
                length(eSI) length(ePO4) length(r) length(pHSCALEIN)...
                length(K1K2CONSTANTS) length(KSO4CONSTANTS)];

    if length(unique(veclengths))>2
            disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end

    % Make column vectors of all input vectors
    PAR1         =PAR1         (:);
    PAR2         =PAR2         (:);
    PAR1TYPE     =PAR1TYPE     (:);
    PAR2TYPE     =PAR2TYPE     (:);
    SAL          =SAL          (:);
    TEMPIN       =TEMPIN       (:);
    TEMPOUT      =TEMPOUT      (:);
    PRESIN       =PRESIN       (:);
    PRESOUT      =PRESOUT      (:);
    SI           =SI           (:);
    PO4          =PO4          (:);
    ePAR1        =ePAR1        (:);
    ePAR2        =ePAR2        (:);
    eSAL         =eSAL         (:);
    eTEMP        =eTEMP        (:);
    eSI          =eSI          (:);
    ePO4         =ePO4         (:);
    r            =r            (:);
    pHSCALEIN    =pHSCALEIN    (:);
    K1K2CONSTANTS=K1K2CONSTANTS(:);
    KSO4CONSTANTS=KSO4CONSTANTS(:);

    % Find the longest column vector:
    ntps = max(veclengths);

    % Populate column vectors
    PAR1(1:ntps,1)          = PAR1(:)          ;
    PAR2(1:ntps,1)          = PAR2(:)          ;
    PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
    PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
    SAL(1:ntps,1)           = SAL(:)           ;
    TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
    TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
    PRESIN(1:ntps,1)        = PRESIN(:)        ;
    PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    ePAR1(1:ntps,1)         = ePAR1(:)         ;
    ePAR2(1:ntps,1)         = ePAR2(:)         ;
    eSAL(1:ntps,1)          = eSAL(:)          ;
    eTEMP(1:ntps,1)         = eTEMP(:)         ;
    eSI(1:ntps,1)           = eSI(:)           ;
    ePO4(1:ntps,1)          = ePO4(:)          ;
    r(1:ntps,1)             = r(:)             ;
    pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
    K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
    KSO4CONSTANTS(1:ntps,1) = KSO4CONSTANTS(:) ;

    % Default value for epK
    if (isempty(epK))
        epK = [0.004, 0.015, 0.03, 0.01, 0.01, 0.02, 0.02, 0.01];
    else
        % Check validity of epK
        if (length(epK) == 1 && epK == 0)
            % this means that the caller does not want to account for errors on dissoc. constants
            epK = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
        elseif (length(epK) != 8)
            error ('invalid parameter epK: ', epK)
        end
    end
    % names of dissociation constants
    Knames = {'K0','K1','K2','Kb','Kw','Kspa', 'Kspc', 'bor'};

    % Convert error on pH to error on [H+] concentration
    % in case where first input variable is pH
    isH = (PAR1TYPE == 3);
    
    pH = PAR1(isH);
    epH = ePAR1(isH);       % Error on pH
    H  = 10.^(-pH);         % H+ concentration
    r(isH) = -r(isH);       % Inverse sign of 'r' if PAR1 is pH

    % dpH = d(-log10[H])
    %     = d(- ln[H] / ln[10] )
    %     = -(1/ln[10]) * d (ln[H])
    %     = -(1/ln[10]) * (dH / H)
    % Thus dH = - ln[1O] * [H] dpH
    eH =  log(10) * (H .* epH);     % Remove the minus sign because all errors (sigmas) are positive by definition
    ePAR1(isH) = eH;

    % Same conversion for second variable
    isH = (PAR2TYPE == 3);
    pH = PAR2(isH);
    epH = ePAR2(isH);       % Error on pH
    H  = 10.^(-pH);         % H+ concentration
    r(isH) = -r(isH);       % Inverse sign of 'r' if PAR2 is pH

    eH =   log(10) * (H .* epH);
    ePAR2(isH) = eH;

    % initialise total square error
    sq_err = zeros(ntps,1);
        
    % Contribution of PAR1 to squared standard error
    if (any (ePAR1 != 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv1, headers] = derivnum ('PAR1',PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        err = deriv1 .* ePAR1;
        sq_err = resize(sq_err, size(err));
        sq_err = sq_err + err .* err;
    end

    % Contribution of PAR2 to squared standard error
    if (any (ePAR2 != 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv2, headers] = derivnum ('PAR2',PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        err = deriv2 .* ePAR2;
        sq_err = resize(sq_err, size(err));
        sq_err = sq_err + err .* err;
    end

    % Contribution of covariance of PAR1 and PAR2 to squared standard error
    if (any (r != 0.0) && any (ePAR1 != 0.0) && any (ePAR2 != 0.0))
        % Compute covariance from correlation coeff. & std deviations
        covariance = r .* ePAR1 .* ePAR2;
        % Contribution to squared error
        err2 = 2 * deriv1 .* deriv2 .* covariance;
        sq_err = sq_err + err2;
    end
    
    % Contribution of Silion (total dissolved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where SI = 0 
    %          because computation of sensitivity to SI fails in that case
    %
    SI_valid = (SI != 0) & (eSI != 0);
    if (any (SI_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, headers] = derivnum ('sil',PAR1(SI_valid),PAR2(SI_valid),PAR1TYPE(SI_valid),PAR2TYPE(SI_valid),...
                   SAL(SI_valid),TEMPIN(SI_valid),TEMPOUT(SI_valid),PRESIN(SI_valid),PRESOUT(SI_valid),...
                   SI(SI_valid),PO4(SI_valid),pHSCALEIN(SI_valid),K1K2CONSTANTS(SI_valid),KSO4CONSTANTS(SI_valid));
        err = deriv .* eSI(SI_valid);
        new_size = [ntps size(err)(2)];
        sq_err = resize(sq_err, new_size);
        sq_err(SI_valid,:) = sq_err(SI_valid,:) + err .* err;
    end

    % Contribution of Phosphorus (total dissoloved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where PO4 = 0 
    %          because computation of sensitivity to PO4 fails in that case
    %
    PO4_valid = (PO4 != 0) & (ePO4 != 0);
    if (any (PO4_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, headers] = derivnum ('phos',PAR1(PO4_valid),PAR2(PO4_valid),PAR1TYPE(PO4_valid),PAR2TYPE(PO4_valid),...
                   SAL(PO4_valid),TEMPIN(PO4_valid),TEMPOUT(PO4_valid),PRESIN(PO4_valid),PRESOUT(PO4_valid),...
                   SI(PO4_valid),PO4(PO4_valid),pHSCALEIN(PO4_valid),K1K2CONSTANTS(PO4_valid),KSO4CONSTANTS(PO4_valid));
        err = deriv .* ePO4(PO4_valid);
        new_size = [ntps size(err)(2)];
        sq_err = resize(sq_err, new_size);
        sq_err(PO4_valid,:) = sq_err(PO4_valid,:) + err .* err;
    end

    % Contribution of T (temperature) to squared standard error
    if (any (eTEMP != 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv, headers] = derivnum ('T',PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        err = deriv .* eTEMP;
        sq_err = resize(sq_err, size(err));
        sq_err = sq_err + err .* err;
    end

    % Contribution of S (salinity) to squared standard error
    if (any (eSAL != 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv, headers] = derivnum ('S',PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        err = deriv .* eSAL;
        sq_err = resize(sq_err, size(err));
        sq_err = sq_err + err .* err;
    end

    % Calculate dissociation constants
    data = CO2SYS (PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
        
    % Calculate [Ca++]
    % '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    % '       this is .010285.*Sali./35
    Ca = 0.02128./40.087.*(SAL./1.80655);% ' in mol/kg-SW

    % Contribution of all pKi to squared standard error
    for i = 1:length(epK)

        % if error on Ki is given
        if (epK(i) != 0.0)
            % Select Ki
            switch i
                case 1
                  Ki = data(:,53);   % K0
                case 2
                  Ki = data(:,54);   % K1
                case 3
                  Ki = data(:,55);   % K2
                case 4
                  Ki = data(:,59);   % KB
                case 5
                  Ki = data(:,58);   % KW
                case 6
                  % Recompute KAr from OmegaAr and ions [Ca++] and [CO3--] concentrations
                  OmegaAr = data(:,16);
                  CO3 = data(:,7) * 1e-6;
                  Ki = CO3.*Ca./OmegaAr;
                case 7
                  % Recompute KCa from OmegaCa and ions [Ca++] and [CO3--] concentrations
                  OmegaCa = data(:,15);
                  CO3 = data(:,7) * 1e-6;
                  Ki = CO3.*Ca./OmegaCa;
                case 8
                  Ki = data(:,79) * 1e-6 / log(10);  % TB   
                  % TB (divide by log(10) to multiply by same just below
                  % (not a pK value, unlike the others) 
            end

            % compute error on Ki from that on pKi
            eKi = - epK(i) * Ki * log(10);
            % Compute sensitivities (partial derivatives)
            [deriv, headers] = derivnum (cell2mat(Knames(1,i)),PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
            err = deriv .* eKi;
            %disp('deriv = '), disp(deriv);
            sq_err = resize(sq_err, size(err));
            sq_err = sq_err + err .* err;
        end
    end

    % Compute and return resulting total error (or uncertainty)
    total_error = sqrt (sq_err);
end
