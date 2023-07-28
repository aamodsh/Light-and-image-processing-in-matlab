% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultAxesFontWeight','bold')

% Change default text fonts.
set(0,'DefaultTextFontname', 'Helvectica')
set(0,'DefaultTextFontSize', 16)
set(gca,'DefaultTextFontWeight','bold')

% Set default figures to 'docked' or 'normal'
set(0,'DefaultFigureWindowStyle','docked')


%Save without borders
set(gca,'position',[0 0 1 1],'units','normalized')
