% PistarSetDSGE
%
% This file estimates spec JPT with TVIT T.
%   
% See preamble of SetDSGE in package VC_BayesianEstimation for information on 
% how to customize this file.
%
% .........................................................................
%
% Created: March 17, 2008 by Vasco Curdia, Ging Cee Ng
% Updated: February 11, 2015 by Vasco Curdia
% 
% Copyright (C) 2008-2015 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Preamble
clear all
tic
ttic = toc();

%% Settings
FileName.Output = 'Pistar';


%% Parallel options
UseParallel = 0;
nMaxWorkers = 4;
if UseParallel,matlabpool('open',nMaxWorkers),end


%% Allow for running only some of the blocks of actions
% Possible Actions:
%   Actions = {'All'};
%   Actions = {'Setup','MaxPost','MCMC'};
Actions = {'All'};


%% ------------------------------------------------------------------------

%% Load framework if already set
if ~any(ismember({'All','Setup'},Actions))
    save NewSettings ttic ListPathDependencies Actions UseParallel nMaxWorkers
    load(FileName.Output)
    load NewSettings
    delete NewSettings.mat
else
    
%% Setup

%% Data
FileName.Data = 'JPTData_1954q3_2009q3_3obs';
nPreSample = 132;
DateLabels.Start = '1987q3';
DateLabels.End = '2009q3'; 
TimeIdx = TimeIdxCreate(DateLabels.Start,DateLabels.End);
DateLabels.XTick = find(ismember(TimeIdx,{'1990q1','1995q1','2000q1','2005q1'}));
DateLabels.XTickLabels = {'1990','1995','2000','2005'};

%% Estimated parameters
Params = {...
    'JPTalpha', 'N', 0.3, 0.05, '\alpha';
    'JPTiotap', 'B', 0.5, 0.15, '\iota_p';
    'JPTiotaw', 'B', 0.5, 0.15, '\iota_w';
    'JPT100gamma', 'N', 0.5, 0.025, '\gamma^*';
    'JPTh', 'B', 0.5, 0.1, 'h';
    'JPTlambdap', 'N', 0.15, 0.05, '\lambda_p';
%    'JPTlambdaw', 'N', 0.15, 0.05, '\lambda_w';
    'JPTlogL', 'N', 0, 0.5, '\log L';
    'JPT100pi', 'N', 0.5, 0.1, '\pi^*';
    'JPT100DiscRate', 'G', 0.25, 0.1, 'DiscRate';
    'JPTnu', 'G', 2, 0.75, '\nu';
    'JPTxip', 'B', 0.66, 0.1, '\xi_p';
    'JPTxiw', 'B', 0.66, 0.1, '\xi_w';
    'JPTchi', 'G', 5, 1, '\chi';
    'JPTSdd', 'G', 4, 1, 'S^{''}';
    'JPTphipi', 'N', 1.7, 0.3,'\phi_{\pi}';
    'JPTphiY', 'N', 0.125, 0.05,'\phi_Y';
    'JPTphidY', 'N', 0.125, 0.05, '\phi_{dY}';
    'JPTrhoR', 'B', 0.6, 0.2, '\rho_R';
    'JPTrhoz','B', 0.6, 0.2, '\rho_z';
    'JPTrhog','B', 0.6, 0.2, '\rho_g';
    'JPTrhomu','B', 0.6, 0.2, '\rho_{\mu}';
    'JPTrhop','B', 0.6, 0.2, '\rho_p';
    'JPTrhow','B', 0.6, 0.2, '\rho_w';
    'JPTrhob','B', 0.6, 0.2, '\rho_b';
%    'JPTrhomp', 'B', 0.4, 0.2, '\rho_{mp}';
    'JPTthetap','B', 0.5, 0.2, '\theta_p';
    'JPTthetaw','B', 0.5, 0.2, '\theta_w';
    %Shocks (100*StDev)
    'JPTsigmamp', 'IG1', 0.1, 1, '\sigma_{mp}';
    'JPTsigmaz','IG1', 0.5, 1, '\sigma_z';
    'JPTsigmag','IG1', 0.5, 1, '\sigma_g';
    'JPTsigmamu','IG1', 0.5, 1, '\sigma_{\mu}';
    'JPTsigmap','IG1', 0.1, 1, '\sigma_p';
    'JPTsigmaw','IG1', 0.1, 1, '\sigma_w';
    'JPTsigmab','IG1', 0.1, 1, '\sigma_b';
    'JPTsigmapistar','IG1', 0.2, 1, '\sigma_{\pi^*}';
    'JPTrhopistar', 'B', 0.95, 0.04, '\rho_{\pi^*}';
    };
    
%% Observation variables
ObsVar = {'dlogXo';'pio';'logRo'};

%% State space variables
StateVar = {...
    % Endogenous variables
    'y';'k';'L';'rho';'w';'pi';'s';'lambda';'c';'R';'u';'phi';'i';'kbar';'gw';'x';
    % Exogenous variables
    'z';'g';'mu';'lambdapstar';'lambdawstar';'bstar';'etamp';'pistar';
    % Natural variables (rn is the real interest rate)
    'yn';'kn';'Ln';'rhon';'wn';'sn';'lambdan';'cn';'rn';'un';'phin';'in';'kbarn';
    'gwn';'xn'
    % Artificial Variables
    'Ep';'Ew';          % ARMAlambda[p/w] vars
    'dx';'dc';'di';'dw';
    'piL';'cL';'iL';'wL';'zL';
    'cnL';'inL';
    'gapn';
    };

%% Shocks
ShockVar = {'emp';'ez';'eg';'emu';'ep';'ew';'eb';'epistar'};

%% Create symbolic variables
GenSymVars

%% Auxiliary Definitions
JPTdelta = 0.025;
JPT1og = 1 - 0.22;
JPTlambdaw=0.25;

%Steady State
JPTbeta = (1+JPT100DiscRate/100)^(-1);
JPTgamma = JPT100gamma/100;
JPTeg = exp(JPTgamma);
JPT100logR = JPT100pi+100*(JPTeg/JPTbeta-1);

JPTrho = JPTeg/JPTbeta-(1-JPTdelta);
JPTw = (1/(1+JPTlambdap)*(JPTalpha^JPTalpha)*((1-JPTalpha)^(1-JPTalpha))*1/JPTrho^JPTalpha)^(1/(1-JPTalpha));
JPTkoL = JPTw/JPTrho*JPTalpha/(1-JPTalpha);
JPTFoL = JPTkoL^JPTalpha-JPTrho*JPTkoL-JPTw;
JPTyoL = JPTkoL^JPTalpha-JPTFoL;
JPTioL = (1-(1-JPTdelta)*exp(-JPTgamma))*JPTeg*JPTkoL;
JPTcoL = JPTyoL*JPT1og-JPTioL;
JPTkoy = JPTkoL/JPTyoL;
JPTFoy = JPTFoL/JPTyoL;
JPTioy = JPTioL/JPTyoL;
JPTcoy = JPTcoL/JPTyoL;

%Equation Parameters
JPTpipiF = JPTbeta/(1+JPTiotap*JPTbeta);
JPTpipiL= JPTiotap/(1+JPTiotap*JPTbeta);
JPTkappa = (1-JPTxip*JPTbeta)*(1-JPTxip)/(JPTxip*(1+JPTiotap*JPTbeta));

JPTlambdaD = (JPTeg-JPTh*JPTbeta)*(JPTeg-JPTh);
JPTheg = JPTh*JPTeg;
JPTlambdacF = JPTbeta*JPTheg/JPTlambdaD;
JPTlambdac = (JPTeg^2+JPTbeta*JPTh^2)/JPTlambdaD;
JPTlambdacL = JPTheg/JPTlambdaD;
JPTlambdaz = (JPTbeta*JPTheg*JPTrhoz-JPTheg)/JPTlambdaD;
JPTlambdab = (JPTeg-JPTh*JPTbeta*JPTrhob)/(JPTeg-JPTh*JPTbeta);

JPTphi1 = (1-JPTdelta)*JPTbeta*JPTeg^(-1);
JPTphi2 = 1-JPTphi1;

JPTe2gSdd = (JPTeg^2)*JPTSdd;

JPTkbar1 = (1-JPTdelta)*JPTeg^(-1);
JPTkbar2 = 1- JPTkbar1;

JPTwwL = 1/(1+JPTbeta);
JPTwwF = JPTbeta/(1+JPTbeta);
JPTkappaw = (1-JPTxiw*JPTbeta)*(1-JPTxiw)/(JPTxiw*(1+JPTbeta)*(1+JPTnu*(1+1/JPTlambdaw)));
JPTwpiL = JPTiotaw/(1+JPTbeta);
JPTwpi = (1+JPTbeta*JPTiotaw)/(1+JPTbeta);
JPTwpiF = JPTbeta/(1+JPTbeta);
JPTwzL = JPTiotaw/(1+JPTbeta);
JPTwz = (1+JPTbeta*JPTiotaw-JPTrhoz*JPTbeta)/(1+JPTbeta);

%Normalized Shocks
lambdap_t = lambdapstar_t/JPTkappa;
lambdaw_t = lambdawstar_t/JPTkappaw;
JPTbnorm = (1-JPTrhob)*(JPTeg-JPTbeta*JPTh*JPTrhob)*(JPTeg-JPTh)/(JPTeg*JPTh+JPTeg^2+JPTbeta*JPTh^2);
b_t = bstar_t/JPTbnorm;

%% Observation Equations
%All observables are quarterly rates, expressed in percentage. Hours are
%log deviations from the mean (1948-end), also expressed in percentage. 
ObsEq = [...
    JPT100gamma*one + dx_t + z_t - dlogXo_t;
%     JPT100gamma*one + dc_t + z_t - dlogCo_t;
%     JPT100gamma*one + di_t + z_t - dlogIo_t;
%     JPTlogL*one + L_t - logLo_t;
%     JPT100gamma*one + dw_t + z_t - dlogWPo_t;
    JPT100pi*one + pi_t - pio_t;
    JPT100logR*one + R_t - logRo_t;
    ];

%% State Equations
% Efficient economy equations ("n")
StateEq = [...
    % Production Function (y,yn)
    (1+JPTFoy)*(JPTalpha*k_t+(1-JPTalpha)*L_t)-y_t;
    (1+JPTFoy)*(JPTalpha*kn_t+(1-JPTalpha)*Ln_t)-yn_t;
    % Cost Minimization (L,Ln)
    w_t+L_t-k_t-rho_t;
    wn_t+Ln_t-kn_t-rhon_t;
    % Marginal Cost (s,sn)
    JPTalpha*rho_t+(1-JPTalpha)*w_t-s_t;
    JPTalpha*rhon_t+(1-JPTalpha)*wn_t-sn_t;
    % Phillips curve (pi,rn)
    JPTpipiF*pi_tF+JPTpipiL*piL_t+JPTkappa*s_t+lambdapstar_t-pi_t;
    JPTkappa*sn_t+0*lambdapstar_t;  %GC: kill lambdapstar for efficient economy
    % Consumption FOC (c,cn)
    JPTlambdacF*c_tF-JPTlambdac*c_t+JPTlambdacL*cL_t+...
        JPTlambdaz*z_t+JPTlambdab/JPTbnorm*bstar_t-lambda_t;
    JPTlambdacF*cn_tF-JPTlambdac*cn_t+JPTlambdacL*cnL_t+...
        JPTlambdaz*z_t+JPTlambdab*bstar_t/JPTbnorm-lambdan_t;
    % Euler Equation (lambda, lambdan)
    R_t+lambda_tF-JPTrhoz*z_t-pi_tF-lambda_t;    
    rn_t+lambdan_tF-JPTrhoz*z_t-lambdan_t;    
    % Capital Utilization FOC (rho,rhon)
    rho_t/JPTchi-u_t;
    rhon_t/JPTchi-un_t;
    % Capital FOC (phi,phin)
    JPTphi1*(phi_tF-JPTrhoz*z_t)+JPTphi2*(lambda_tF-JPTrhoz*z_t+rho_tF)-phi_t;
    JPTphi1*(phin_tF-JPTrhoz*z_t)+JPTphi2*(lambdan_tF-JPTrhoz*z_t+rhon_tF)-phin_t;   
    % Investment FOC (i,in)
    phi_t+mu_t-JPTe2gSdd*(i_t-iL_t+z_t)+JPTbeta*JPTe2gSdd*(i_tF-i_t+JPTrhoz*z_t)-lambda_t;
    phin_t+mu_t-JPTe2gSdd*(in_t-inL_t+z_t)+JPTbeta*JPTe2gSdd*(in_tF-in_t+JPTrhoz*z_t)-lambdan_t;   
    % Capital Input (k,kn)
    u_t+kbar_tL-z_t-k_t;
    un_t+kbarn_tL-z_t-kn_t;
    % Capital Accumulation (kbar,kbarn)
    JPTkbar1*(kbar_tL-z_t)+JPTkbar2*(mu_t+i_t)-kbar_t;    
    JPTkbar1*(kbarn_tL-z_t)+JPTkbar2*(mu_t+in_t)-kbarn_t;
    % Wage Phillips Curve (w,wn)
    JPTwwL*wL_t+JPTwwF*w_tF-JPTkappaw*gw_t+...
        JPTwpiL*piL_t-JPTwpi*pi_t+JPTwpiF*pi_tF+...
        JPTwzL*zL_t-JPTwz*z_t+lambdawstar_t-w_t;
    -JPTkappaw*gwn_t+0*lambdawstar_t;  %GC: kill lambdawstar for efficient economy
    % Wage Gap (gw,gwn)
    w_t-(JPTnu*L_t+b_t-lambda_t)-gw_t;
    wn_t-(JPTnu*Ln_t+b_t-lambdan_t)-gwn_t;
    % Monetary Policy Rule (R)
    JPTrhoR*R_tL+(1-JPTrhoR)*(pistar_t+JPTphipi*(pi_t-pistar_t)+JPTphiY/4*gapn_t)+etamp_t-R_t;
    % Definition of GDP (x,xn)
    y_t-JPTrho*JPTkoy*u_t-x_t;
    yn_t-JPTrho*JPTkoy*un_t-xn_t;
    % Market Clearing
    JPT1og*g_t+JPTcoy*c_t+JPTioy*i_t+JPTrho*JPTkoy*u_t-JPT1og*y_t;
    JPT1og*g_t+JPTcoy*cn_t+JPTioy*in_t+JPTrho*JPTkoy*un_t-JPT1og*yn_t;
    % Shocks
    JPTrhoz*z_tL+JPTsigmaz*ez_t-z_t;
    JPTrhog*g_tL+JPTsigmag*eg_t-g_t;
    JPTrhow*lambdawstar_tL+Ew_t-JPTthetaw*Ew_tL-lambdawstar_t;
    JPTsigmaw*ew_t-Ew_t;
    JPTrhomu*mu_tL+JPTsigmamu*emu_t-mu_t;
    JPTrhop*lambdapstar_tL+Ep_t-JPTthetap*Ep_tL-lambdapstar_t;
    JPTsigmap*ep_t-Ep_t;
    JPTrhob*bstar_tL+JPTsigmab*eb_t-bstar_t;
    JPTsigmamp*emp_t-etamp_t;
    JPTrhopistar*pistar_tL+JPTsigmapistar*epistar_t-pistar_t;
    % Lagged Variables
    x_t-x_tL-dx_t;
    c_t-c_tL-dc_t;
    i_t-i_tL-di_t;
    w_t-w_tL-dw_t;
    pi_tL-piL_t;
    c_tL-cL_t;
    i_tL-iL_t;
    w_tL-wL_t;
    z_tL-zL_t;
    cn_tL-cnL_t;
    in_tL-inL_t;
    x_t-xn_t-gapn_t;
    ];

%% Data
DataAnalysis

%% Generate Matrices for spec
MakeMats

%% Priors
PriorAnalysis

%% Generate posterior function
GenPost

%% end Setup Action if
TimeElapsed.Setup = toc();
fprintf('\n%s\n\n',vctoc([],TimeElapsed.Setup))
save(FileName.Output)
end

%% ------------------------------------------------------------------------

%% Use Starting Vector
guess=[0.2 0.21	0.15 0.53 0.85 0.25 0.12 0.5 0.75 0.1 4 0.85 0.75 5 2.5 2 0.1 ...
       0.25 0.8 0.25	0.98 0.7 0.95 0.98 0.7 0.15 0.75 0.95 0.2 0.9 0.35 5 ...
       0.15 0.2 0.05 0.2349 0.9896];
guess=[guess(1:6),guess(8:25),guess(27:end)];

%% MaxPost
if any(ismember({'All','MaxPost'},Actions))
    nMax = 20;
    MinParams.H0 = diag([Params(:).priorse].^2);
    MinParams.crit = 1e-8;
    MinParams.nit = 1000;
    MinParams.Ritmax = 30;
    MinParams.Ritmin = 10;
    MinParams.RH0 = 1;
    MinParams.Rscaledown = 0.5;
    MinParams.Guess.x0 = guess';
    DrawAll=1;
    MaxPost
    save(FileName.Output)
    save([FileName.Output,'MaxPost'])
end

%% MCMC
if any(ismember({'All','MCMC'},Actions))
    nChains = 4;
    nDrawsSearch = 5000;
    dscale = [0.2,0.05,0.01];
    BurnIn = 0.25;
    nThinning = 2;
    % nDraws = 100000;
    for nUpdate=0:1
        fprintf('\n*****************')
        fprintf('\n* MCMC Update %.0f *',nUpdate)
        fprintf('\n*****************')
        if nUpdate==0
            nDraws = 100000;
        else
            nDraws = 200000;
        end
        MCMCOptions.ScaleJumpFactor = 2.4;
        MCMCSearchScaleFactor
        save(sprintf('%sMCMCUpdate%.0f_SSF',FileName.Output,nUpdate))
        MCMC
        save(FileName.Output)
        save(sprintf('%sMCMCUpdate%.0f',FileName.Output,nUpdate))
        delete(sprintf('%sMCMCUpdate%.0f_SSF.mat',FileName.Output,nUpdate))
        MCMCAnalysis
        save(FileName.Output)
        save(sprintf('%sMCMCUpdate%.0f',FileName.Output,nUpdate))
    end
end

%% ------------------------------------------------------------------------

%% Close matlabpool
if UseParallel,matlabpool close,end

%% elapsed time
fprintf('\n%s\n\n',vctoc(ttic))

%% Save environment
save(FileName.Output)

%%------------------------------------------------------------------------
