%% Figure Configuration

hight = 18; width = 32;
top = 0.5; bottom = 1; left = 1.5; right = 0.5;

set(0,'defaultFigureUnits','centimeters');
set(0,'defaultFigurePosition',[5 3 width hight]);
set(0,'defaultAxesLineWidth',0.5);

set(0,'defaultAxesGridLineStyle',':');
set(0,'defaultAxesYGrid','on');
set(0,'defaultAxesXGrid','on');

set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultAxesFontSize',16);

set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultTextFontSize',16);

set(0,'defaultAxesUnits','normalized');
set(0,'defaultAxesPosition',[left/width bottom/hight (width-left-right)/width  (hight-bottom-top)/hight]);
set(0,'defaultLineLineWidth',1);
set(0,'defaultAxesColorOrder',[0 0 0]);
set(0,'defaultAxesTickDir','out');

set(0,'defaultFigurePaperPositionMode','auto');

% set(0,'defaultLegendLocation','best');
% set(0,'defaultLegendBox','on');
% set(0,'defaultLegendOrientation','horizontal');

color_2014b = [ 0         0.4470    0.7410; % blue
                0.8500    0.3250    0.0980; % orange
                0.9290    0.6940    0.1250; % yellow
                0.4940    0.1840    0.5560; % purple
                0.4660    0.6740    0.1880; % green
                0.3010    0.7450    0.9330; % cyan
                0.6350    0.0780    0.1840; % red
                ];
            
color_2014b_blue   = color_2014b(1,:);
color_2014b_orange = color_2014b(2,:);
color_2014b_yellow = color_2014b(3,:);  
color_2014b_purple = color_2014b(4,:);
color_2014b_green  = color_2014b(5,:);  
color_2014b_cyan   = color_2014b(6,:);
color_2014b_red    = color_2014b(7,:);  

