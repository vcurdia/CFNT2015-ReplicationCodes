% GrowthPistarSetDSGE
%
% This file estimates spec TVIT T with Growth.
%   
% See preamble of SetDSGE in package VC_BayesianEstimation for information on 
% how to customize this file.
%
% .........................................................................
%
% Created: March 17, 2008 by Vasco Curdia, Ging Cee Ng
% Updated: February 10, 2015 by Vasco Curdia
% 
% Copyright (C) 2008-2015 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Preamble
clear all
tic
ttic = toc();

%% Settings
FileName.Output = 'GrowthPistar';

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
FileName.Data = 'data_greenspan_bernanke_20091204';
nPreSample = 0;
DateLabels.Start = '1987q3';
DateLabels.End = '2009q3'; 
TimeIdx = TimeIdxCreate(DateLabels.Start,DateLabels.End);
DateLabels.XTick = find(ismember(TimeIdx,{'1990q1','1995q1','2000q1','2005q1'}));
DateLabels.XTickLabels = {'1990','1995','2000','2005'};

%% Estimated parameters
Params = {...
    'omega', 'G', 1, 0.2,'\omega';
    'xi', 'G', 0.1, 0.05,'\xi';
    'eta', 'B', 0.6, 0.2,'\eta';
    'zeta', 'B', 0.6, 0.2,'\zeta';
    'rho', 'B', 0.7, 0.15,'\rho';
    'phipi', 'N', 1.5, 0.25,'\phi_\pi';
    'phidy', 'N', 0.5, 0.2,'\phi_{\Delta y}';
    'pistar', 'N', 2, 1,'\pi^*';
    'ra', 'N', 2, 1,'r^a';
    'gammaa', 'N', 3, .35,'\gamma^a';
    'rhodelta', 'B', 0.5, 0.2,'\rho_\delta';
    'rhogamma', 'B', 0.5, 0.2,'\rho_\gamma';
    'rhou', 'B', 0.5, 0.2,'\rho_u';
    'rhopistar', 'B', 0.95, 0.04,'\rho_{\pi^*}';
    'sigmadelta', 'IG1', 0.5, 2,'\sigma_\delta';
    'sigmagamma', 'IG1', 0.5, 2,'\sigma_\gamma';
    'sigmau', 'IG1', 0.5, 2,'\sigma_u';
    'sigmai', 'IG1', 0.5, 2,'\sigma_i';
    'sigmapistar', 'IG1', 0.2, 1,'\sigma_{\pi^*}';
    };
zeta=[];

%% Observation variables
ObsVar = {'gYa';'pia';'ira'};

%% State space variables
StateVar = {...
    % Regular variables
    'xtil';'YA';'pitil';'pi';'ir';'r';
    'xe';'re';'YAe';
    'delta';'gamma';'u';'pistar';
    % a couple of artificial variables
    'YAL';'YAeL';
    };

%% Shocks
ShockVar = {'edelta';'egamma';'eu';'ei';'epistar'}; 

%% create symbolic variables
GenSymVars

%% Auxiliary definitions
beta = 0.99;
gamma = gammaa/400;
r = ra/400;
phigammatil = exp(gamma)/(exp(gamma)-beta*eta);
etagammatil = exp(gamma)/(exp(gamma)-eta);
phigamma = phigammatil*etagammatil;
etagamma = eta/exp(gamma);

%% Observational Equations
ObsEq = [...
    gammaa*one+400*(YA_t-YAL_t+gamma_t) - gYa_t;
    pistar*one+400*pi_t - pia_t;
    (ra+pistar)*one+400*ir_t - ira_t];

%% State Equations
StateEq = [...
    % IS Block
    xtil_tF-phigamma^(-1)*(ir_t-pi_tF-re_t)-xtil_t;
    (xe_t-etagamma*(YAL_t-YAeL_t))-beta*etagamma*(xe_tF-etagamma*xe_t)-xtil_t;
    ir_t-pi_tF - r_t;
    % Efficient Rates
    YA_t-YAe_t-xe_t;
    YAe_tF-omega^(-1)*(gamma_tF-re_t+delta_tF)-YAe_t;
    -phigamma*(YAe_t-etagamma*(YAeL_t-gamma_t)-...
        beta*etagamma*(YAe_tF+gamma_tF-etagamma*YAe_t))+...
        beta*etagamma/(1-beta*etagamma)*delta_tF-omega*YAe_t;
    % PC Block
    beta*pitil_tF+xi*(omega*xe_t+phigamma*xtil_t)+u_t-pitil_t;
    pi_t-zeta*pi_tL-pitil_t;
    % Policy Rule
    rho*ir_tL+(1-rho)*(pistar_t+phipi*(pi_t-pistar_t)+phidy*(YA_t-YA_tL+gamma_t))+sigmai/400*ei_t - ir_t;
    % Shocks
    rhodelta*delta_tL+sigmadelta/400*edelta_t - delta_t;
    rhogamma*gamma_tL+sigmagamma/400*egamma_t - gamma_t;
    rhou*u_tL+sigmau/400*eu_t - u_t;
    rhopistar*pistar_tL+sigmapistar/400*epistar_t - pistar_t;
    % Auxiliary equations
    YAL_t - YA_tL; % defines YA_tL
    YAeL_t-YAe_tL;
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
    MaxPost
    save(FileName.Output)
    save([FileName.Output,'MaxPost'])
end

%% MCMC
if any(ismember({'All','MCMC'},Actions))
    nChains = 4;
    nDrawsSearch = 1000;
    dscale = [0.2,0.05,0.01];
    BurnIn = 0.25;
    nThinning = 1;
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
