function Mats=ReMats(x,isObsMats,isEvalGensys,fid,verbose,varargin)

% Options:
%   1. [Mats]=ReMats(x)
%   2. [Mats]=ReMats(x,isObsMats)
%   3. [Mats]=ReMats(x,isObsMats,isEvalGensys)
%   4. [Mats]=ReMats(x,isObsMats,isEvalGensys,fid)
%   5. [Mats]=ReMats(...,varargin)
%   Default: Option 1 where Mats is a structure with fields 'ObsEq','StateEq', and 'REE'
%
% Inputs:
%   x
%   Vector of parameters
%
%   isObsMats (Optional)
%   Indicates whether to output Observable matrices and/or Gensys matrices.
%   0: No observable matrices.
%   1: Function outputs observable matrices ONLY--no gensys matrices.
%   2: Function outputs observable matrices and gensys matrices, according to 'isEvalGensys' options.
%   Default: 2
%
%   isEvalGensys (Optional)
%   Flags whether to evalute gensys and determines what function will output from gensys.
%   0: Do not evaluate gensys. Function outputs gensys input matrices only.
%   1: Evalute gensys. Function outputs gensys output matrices only.
%   2: Evalute gensys. Function outputs both gensys input and output matrices.
%   Default: 1
%
%   varargin (optional)
%   Allows for arguments for verbose and flie id to pass to gensys if needed.
% Outputs:
%   Mats is a structure with fields dependent on the values of the input flags:
%      ObsEq.HBar,ObsEq.H
%          Transition matrices from observable equation
%
%      REE.GBar,REE.G1,REE.G2,REE.eu
%          Gensys output matrices, including the error code vector eu.
%
%      StateEq.Gamma0,StateEq.Gamma1,StateEq.GammaBar,StateEq.Gamma2,StateEq.Gamma3
%          Gensys inputs
%
% Gensys requires the following form:
%   Gamma0*StateVar_t = GammaBar + Gamma1*StateVar_t_1 + Gamma2*ShockVar_t + Gamma3*eta_t
% and yields system in form:
%   StateVar_t = GBar + G1*StateVar_t_1 + G2*ShockVar_t
%
% The observation equations are given by:
%   ObsVar_t = HBar + H*StateVar_t
%
% Created: 2015/2/3 with MakeMats.m 

%% Set default input values
if nargin<2, isObsMats=2; end
if nargin<3, isEvalGensys=1; end
if nargin<4, fid = 1; end

if nargin<5, verbose = 0; end

%% Map parameters
omega = x(1);
xi = x(2);
eta = x(3);
zeta = x(4);
rho = x(5);
phipi = x(6);
phix = x(7);
pistar = x(8);
ra = x(9);
gammaa = x(10);
rhodelta = x(11);
rhogamma = x(12);
rhou = x(13);
sigmadelta = x(14);
sigmagamma = x(15);
sigmau = x(16);
sigmai = x(17);

%% Create Gensys Inputs
GammaBar = [...
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    ];

Gamma0 = [...
    -1, 0, 0, -((eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400)))/exp(gammaa/200), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, (99*eta)/(100*exp(gammaa/400)), 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1, 1/omega, 1/omega, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -(99*eta*exp(gammaa/200))/(100*exp(gammaa/400)*(eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400))), (99*eta)/(100*exp(gammaa/400)*((99*eta)/(100*exp(gammaa/400)) - 1)), -(99*eta*exp(gammaa/200))/(100*exp(gammaa/400)*(eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400))), 0, 0, 0, 0;
    0, 0, -99/100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    ];

Gamma1 = [...
    -1, 0, 0, 0, -((eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400)))/exp(gammaa/200), 0, 0, ((eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400)))/exp(gammaa/200), 0, 0, 0, 0, 0, 0, 0;
    -1, 0, 0, 0, 0, 0, (99*eta^2)/(100*exp(gammaa/200)) + 1, 0, 0, 0, 0, 0, -eta/exp(gammaa/400), eta/exp(gammaa/400), 0;
    0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1/omega, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, - omega - (exp(gammaa/200)*((99*eta^2)/(100*exp(gammaa/200)) + 1))/((eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400))), 0, -(eta*exp(gammaa/200))/(exp(gammaa/400)*(eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400))), 0, 0, (eta*exp(gammaa/200))/(exp(gammaa/400)*(eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400))), 0;
    (xi*exp(gammaa/200))/((eta - exp(gammaa/400))*((99*eta)/100 - exp(gammaa/400))), 0, -1, 0, 0, 0, omega*xi, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -phipi*(rho - 1), -1, 0, -(phix*(rho - 1))/4, 1 - rho, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    ];

Gamma4 = [...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -zeta, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, rho, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, rhodelta, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rhogamma, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, rhou, 0, 0, 0;
    0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    ];

Gamma2 = [...
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, sigmai/400;
    sigmadelta/400, 0, 0, 0;
    0, sigmagamma/400, 0, 0;
    0, 0, sigmau/400, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, 0, 0;
    ];

Gamma3 = eye(15);

cv = (all(Gamma0(1:15,:)==0,2)~=0);
Gamma0(cv,:) = -Gamma1(cv,:);
Gamma1(cv,:) = Gamma4(cv,:);
Gamma3(:,cv) = [];
if ~all(all(Gamma4(~cv,:)==0,2))
  error('Incorrect system reduction')
end

%% Generate matrices for observation equations
if isObsMats>0
   HBar = [...
    gammaa;
    pistar;
    pistar + ra;
    ];
   H = [...
    0, 400, 0, 0, 0, 0, 0, 0, 0, 0, 400, 0, -400, 0, 0;
    0, 0, 0, 400, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 400, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    ];
end

%% Run gensys
if isEvalGensys>0
  [G1,GBar,G2,fmat,fwt,ywt,gev,eu]=vcgensys(Gamma0,Gamma1,GammaBar,...
  Gamma2,Gamma3,'CS',fid,verbose,varargin{:});
end
if isObsMats>0 && isEvalGensys>0
if all(GBar(:)==0)
  StateVarBar = zeros(15,1);
else
  StateVarBar = (eye(15)-G1)\GBar;
end
ObsVarBar = HBar + H*StateVarBar;

s00 = zeros(15,1);

[sig00,sig00rc]=lyapcsd(G1,G2*G2');
sig00 = real(sig00);sig00 = (sig00+sig00')/2;
if sig00rc~=0
  if verbose
    fprintf(fid,'Warning: Could not find unconditional variance!\n');
  end

end

end

%% Assign Output
if isObsMats>0
  Mats.ObsEq.HBar = HBar;
  Mats.ObsEq.H = H;
end
if isObsMats>0 && isEvalGensys>0
  Mats.KF.StateVarBar = StateVarBar;
  Mats.KF.ObsVarBar = ObsVarBar;
  Mats.KF.s00 = s00;
  Mats.KF.sig00 = sig00;
  Mats.KF.sig00rc = sig00rc;
end
if isObsMats~=1 && isEvalGensys>0
  Mats.REE.GBar=GBar;
  Mats.REE.G1=G1;
  Mats.REE.G2=G2;
  Mats.REE.eu=eu;
end
if isObsMats~=1 && isEvalGensys~=1
   Mats.StateEq.Gamma0=Gamma0;
   Mats.StateEq.Gamma1=Gamma1;
   Mats.StateEq.GammaBar=GammaBar;
   Mats.StateEq.Gamma2=Gamma2;
   Mats.StateEq.Gamma3=Gamma3;
end
