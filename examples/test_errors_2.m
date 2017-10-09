% Test of num_deriv() with scalar input
function [result, headers] = test_errors_2 ()
    addpath ("~/Documents/CO2sys")
    
   % 2nd test : 6 records at increasing depth  with NO NUTRIENTS
   % ----------------------------------------

    ePAR1 = 2;    % Alk
    ePAR2 = 2;    % DIC
    eSAL = 0;
    eTEMP = 0;
    eSI = 4;
    ePO4 = 0.1;
 
    PAR1 = 2295;    % Alk
    PAR2 = 2154;    % DIC
    PAR1TYPE = 1;
    PAR2TYPE = 2;
    SAL = 35;
    TEMPIN = 2;
    TEMPOUT = 2; 
    PRESIN = [0:1000:5000];
    PRESOUT = PRESIN;
    SI = 0;
    PO4 = 0;
    pHSCALEIN = 1;   % total scale
    K1K2CONSTANTS = 10; % Lueker 2000
    KSO4CONSTANTS = 3;  % KSO4 of Dickson & TB of Lee 2010

    % With no errors on Ks
    epK=0;
    [result, headers] = errors (PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
                     ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);

    % With default errors on Ks
    epK = '';
    [result, headers] = errors (PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
                     ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,epK,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANTS);
end

