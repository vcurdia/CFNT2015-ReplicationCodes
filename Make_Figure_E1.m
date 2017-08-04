% Make_Figure_E1
%
% Makes Figure E1 of the appendix.
%
% .........................................................................
%
% Created: November 21, 2013 by Vasco Curdia
% Updated: February 17, 2015 by Vasco Curdia
%
% Copyright (C) 2013-2015 by Vasco Curdia

%% -----------------------------------------------------------------------------

%% Preamble
clear all
tic

%% Settings
Vars2Show = {'re'};
Specs2Show = {...
    'Baseline/NoGapRe';'Baseline/Baseline';
    'Baseline/HP';'Baseline/Growth';
    'Baseline/NoGapRePistar';
    'Baseline/EPiNoGapRePistar';
    'Baseline/Pi4qNoGapRePistar';
    'JPT/NoGapRe';
    };
SpecLabels = {...
    'W';'T';
    'T (HP)';'T (\Delta y)';
    'W (TVIT)';
    'W (Forward-Looking, TVIT)';
    'W (4Q Inflation, TVIT)';
    'W (JPT)';
    };


T = 89;
Fig.XTick = 11:20:71;
Fig.XTickLabel = {'1990','1995','2000','2005'};

isLegend = 1;
KeepEPS = 0;
OpenPDF = 1;
ShowFFR = 0;


Fig.Visible = 'off';
Fig.Plot.LineColor = vcColorScheme;
% Fig.Plot.LineStyle = {'-','--','-.',':',':',':'};
Fig.Plot.LineStyle = {'-','--','-.',':','-','--','-.','-'};
Fig.Plot.LineMarker = {'none','none','none','o','none','none','none','none'};
Fig.Plot.LineWidth = [1.5,2,2,2,1.25,1.5,1.5,1];
Fig.Plot.LineMarkerSize = [1,1,1,1.5,1,1,1,1];
Fig.Plot.FontSize = 8;
Fig.Plot.ShowLegend = 1;
Fig.Plot.LegendString = SpecLabels;
Fig.LegendLocation = 'SW';

%% -----------------------------------------------------------------------------

%% Generate series for plot
nSpecs = length(Specs2Show);
for jS=1:nSpecs
    Sj = Specs2Show{jS};
    idx = find(ismember(Sj,'/'));
    PlotVar(Sj(idx+1:end),Sj(1:idx-1))
end

%% Plot
nVars2Show = length(Vars2Show);
tid = 1:T;
for jVar=1:nVars2Show
    Varj = Vars2Show{jVar};
    Figj = Fig;
    PlotData = zeros(nSpecs,T);
    for jS=1:nSpecs
        s = load([Specs2Show{jS},'_Plot_',Varj],'PlotData');
        PlotData(jS,:) = s.PlotData;
    end
    if strcmp(Varj,'re') && ShowFFR
        load('Baseline\Baseline','Data')
        FFRdem = Data(:,3) - mean(Data(:,3));
        FFRdem = FFRdem;
        PlotData = [FFRdem';PlotData];
        for jS=1:nSpecs
            Figj.Plot.LegendString{jS} = ['r^e (',Fig.Plot.LegendString{jS},')'];
        end
        Figj.Plot.LegendString = {'FFR (demeaned)',Figj.Plot.LegendString{:}};
    end
    vcFigure(PlotData,Figj)
    load('Recessions_1987q3_2009q3')
    vcRecessionShades(RecessionDummy)
    vcPrintPDF('Figure_E1',KeepEPS,OpenPDF)
end

%% -----------------------------------------------------------------------------

%% close figures
if strcmp(Fig.Visible,'off')
    close all
end

%% Show time taken
vctoc

%% -----------------------------------------------------------------------------
