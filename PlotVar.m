function PlotVar(Spec,Model)

% PlotVar
%
% Plots the evolution of variables with the recession shades.
% Plots selected Data variables
%
% .........................................................................
%
% Created: February 23, 2009 by Vasco Curdia
% Updated: February 17, 2015 by Vasco Curdia
%
% Copyright (C) 2009-2015 by Vasco Curdia

%% ------------------------------------------------------------------------

%% Preamble

if ~exist('Spec','var'), Spec = 'Baseline'; end
if ~exist('Model','var'), Model = 'Baseline'; end

fprintf('generate series for  %s/%s\n',Model,Spec)
load([Model,'/',Spec])

if strcmp(Model,'Baseline')
    if ismember(Spec,{'NoGapRe'})
        load([Model,'/',Spec,'MCMCDrawsUpdate2Redux'],'xd')
    else
        load([Model,'/',Spec,'MCMCDrawsUpdate1Redux'],'xd')
    end
end        
if strcmp(Model,'JPT')
    if ismember(Spec,{'Baseline'})
        load([Model,'/',Spec,'MCMCDrawsUpdate3Redux'],'xd')
    else
        load([Model,'/',Spec,'MCMCDrawsUpdate2Redux'],'xd')
    end
end        

load('Recessions_1987q3_2009q3')

Vars2Show = {'re','xe'};
Data2Show = {'FFR',''};
PlotLabels = {...
  'r^e','FFR (demeaned)';
  'x^e','';
  };
if strcmp(Model,'Baseline')
    PlotScale = [400,100];
else
    PlotScale = [1,1];
end
ShowLegend = [1,1];

nDrawsStates = 1000; %1000

isDrawStates = 1;
SecondaryAxis = 0;
ShowFig = 0;
isLegend = 1;
KeepEPS = 0;
OpenPDF = 0;
Bands2Show=[50,70,90];
useDateLabels = (exist('DateLabels','var') && isfield(DateLabels,'XTick')...
    && isfield(DateLabels,'XTickLabels'));
LineWidth = [2,2];

%% DEMEAN FFR data
if strcmp(Model,'Baseline')
    FFRdem = Data(:,3) - mean(Data(:,3));
    FFRdem = FFRdem;
else
    FFRdem = Data(nPreSample+1:end,3) - mean(Data(nPreSample+1:end,3));
    FFRdem = 4*FFRdem;
end


%% settings for parallel computing
UseParallel = 0;
nMaxWorkers = 4;

%% ------------------------------------------------------------------------

%% load the mcmc draws
xd = xd(:,randi(size(xd,2),1,nDrawsStates));
nDrawsUsed = size(xd,2);

%% draw states
StateVard = zeros(nStateVar,T,nDrawsStates);
MatFileName = FileName.Mats;
sd = sum(100*clock);
cd(Model)
for j=1:nDrawsUsed
    StateVard(:,:,j) = MakePlotsStatesFcn(xd(:,j),MatFileName,Data,nStateVar,...
                                          nObsVar,nShockVar,T,sd+j,isDrawStates);
end
cd ..

%% Draw plots
nVars2Show = length(Vars2Show);
tid = 1:T;
for jF=1:nVars2Show
    Vars2ShowF = Vars2Show{jF};
    if strcmp(Model,'JPT')
        Vars2ShowF = strrep(Vars2ShowF,'re','rn');
        Vars2ShowF = strrep(Vars2ShowF,'xe','xn');
    end
    idxVar = find(ismember(StateVar,Vars2ShowF));
    if ShowFig
        figure
    else
        figure('Visible','off')
    end
    PlotData = squeeze(PlotScale(jF)*StateVard(idxVar,nPreSample+1:end,:))';
    if strcmp(Data2Show{jF},'FFR')
        AltData = FFRdem';
    else
        AltData = [];
    end
    h = vcPlotDistBands(PlotData,...
                        'AltData',AltData,'Bands2Show',Bands2Show,...
                        'LineWidth',LineWidth,...
                        'ShowLegend',ShowLegend(jF),'LegendLocation','NE',...
                        'LegendString',PlotLabels(jF,:));
    if useDateLabels
        set(gca,...
            'XTick',DateLabels.XTick,'XTickLabel',DateLabels.XTickLabels,...
            'FontSize',8)
    end
    if ShowLegend(jF) && h.Options.ShowLight
        ny = (1+(~isempty(AltData)));
        set(h.LegendObj(ny+1:end),'Visible','off')
        for jy=1:ny
            set(h.LegendLines(:,jy),'XData',...
                              get(h.LegendObj(ny+2*(jy-1)+1),'XData'))
        end
    end
    vcRecessionShades(RecessionDummy)
    vcPrintPDF(sprintf('%s/%s_Plot_%s',Model,FileName.Output,Vars2Show{jF}),...
               KeepEPS,OpenPDF)
    PlotData = prctile(PlotData,50,1);
    fnSave = sprintf('%s/%s_Plot_%s',Model,FileName.Output,Vars2Show{jF});
    save(fnSave,'PlotData')
    fprintf('Saved %s to %s\n',Vars2Show{jF},fnSave);
end

%% close figures
if ~ShowFig, close all, end

%% ------------------------------------------------------------------------
