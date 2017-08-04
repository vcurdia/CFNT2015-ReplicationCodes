% Make_Figure_3
%
% Makes Figure 3 of the paper.
%
% .........................................................................
%
% Created: February 23, 2009 by Vasco Curdia
% Updated: February 17, 2015 by Vasco Curdia
%
% Copyright (C) 2009-2015 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Preamble
clear all
tic


% Generate individual series for Baseline model
PlotVar('Re','Baseline')
% PlotVar('Re','JPT')

load('JPT/Re')
load('JPT/ReMCMCDrawsUpdate2Redux','xd')

load('Recessions_1987q3_2009q3')

AltSpec = {'Baseline/Re'};

Vars2Show = {'re'};
Data2Show = {'',''};
% PlotLabels = {...
%   'JPT: W&T','JPT: W','JPT: T','Baseline: W&T';
%   };
PlotLabels = {...
  'JPT','Baseline';
  };
PlotScale = [1,1];
ShowLegend = [1,1];

nDrawsStates = 1000; %1000

isDrawStates = 1;
SecondaryAxis = 0;
ShowFig = 0;
isLegend = 1;
KeepEPS = 0;
OpenPDF = 1;
Bands2Show=[50,70,90];
useDateLabels = (exist('DateLabels','var') && isfield(DateLabels,'XTick')...
    && isfield(DateLabels,'XTickLabels'));
useDateLabels = (exist('DateLabels','var') && isfield(DateLabels,'XTick')...
    && isfield(DateLabels,'XTickLabels'));
  
ShadeColor = [0.6,0.7,0.8]; %[0.72,0.77,0.82]*0.95;
ShadeFactors = [0.2,0.5];%[0.1,0.65];
LineColor = vcColorScheme;
LineColor(2,:) = [0.25,0.25,0.25];
% LineColor(2,:) = [0,0.45,0.45];

LineStyle = {'-','--','-.','-',':',':',':'};
% LineStyle = {'-','-',':',':',':'};
LineMarker = {'none','none','none','none','o','s','^'};
LineMarkerSize = [1,1,1,1,3,3,3];
LineWidth = [3,2,3,1.5];
% LineWidth = [3,1.5];

%% DEMEAN FFR data
FFRdem = Data(:,3) - mean(Data(:,3));
FFRdem = FFRdem/100;


%% settings for parallel computing
UseParallel = 0;
nMaxWorkers = 4;

%% ------------------------------------------------------------------------

%% load the mcmc draws
xd = xd(:,randi(size(xd,2),1,nDrawsStates));
nDrawsUsed = size(xd,2);

%% draw states
StateVard = zeros(nStateVar,T,nDrawsUsed);
MatFileName = FileName.Mats;
sd = sum(100*clock);
cd JPT
for j=1:nDrawsUsed
    StateVard(:,:,j) = MakePlotsStatesFcn(xd(:,j),MatFileName,Data,nStateVar,...
                       nObsVar,nShockVar,T,sd+j,isDrawStates);
end
cd ..
    
%% Draw plots
nVars2Show = length(Vars2Show);
T = length(TimeIdx);
tid = 1:T;
nAltSpec = length(AltSpec);
jF=1;
Vars2ShowF = Vars2Show{jF};
idxVar = find(ismember(StateVar,strrep(Vars2ShowF,'re','rn')));
if ShowFig
    figure
else
    figure('Visible','off')
end
PlotData = squeeze(PlotScale(jF)*StateVard(idxVar,nPreSample+1:end,:))';
AltData = zeros(nAltSpec,T);
for jAlt=1:nAltSpec
    AltDataj = load([AltSpec{jAlt},'_Plot_',Vars2ShowF]);
    AltData(jAlt,:) = AltDataj.PlotData;
end
if strcmp(Data2Show{jF},'FFR')
    AltData(nAltSpec+1,:) = 100*FFRdem';
end
h = vcPlotDistBands(PlotData,...
                    'AltData',AltData,'Bands2Show',Bands2Show,...
                    'LineWidth',LineWidth,'LineColor',...
                    LineColor,'LineStyle',LineStyle,...
                    'LineMarker',LineMarker,'LineMarkerSize',LineMarkerSize,...
                    'ShadeColor',ShadeColor,'ShadeFactors',ShadeFactors,...
                    'ShowLegend',ShowLegend(jF),'LegendLocation','NE',...
                    'LegendString',PlotLabels(jF,:));
if useDateLabels
    set(gca,...
        'XTick',DateLabels.XTick,'XTickLabel',DateLabels.XTickLabels,...
        'FontSize',8)
end
if ShowLegend(jF) && h.Options.ShowLight
    ny = 1+size(AltData,1);
    set(h.LegendObj(ny+1:end),'Visible','off')
    for jy=1:ny
        set(h.LegendLines(:,jy),'XData',...
                          get(h.LegendObj(ny+2*(jy-1)+1),'XData'))
    end
end
vcRecessionShades(RecessionDummy)
fnSave = 'Figure_3';
vcPrintPDF(fnSave,KeepEPS,OpenPDF)
PlotData = prctile(PlotData,50,1);
save(fnSave,'PlotData')
fprintf('Saved %s to %s\n',Vars2Show{jF},fnSave);

%% close figures
if ~ShowFig, close all, end

%% Show time taken
vctoc

%% ------------------------------------------------------------------------
