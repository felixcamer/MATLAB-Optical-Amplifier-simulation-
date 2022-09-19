clc;
clear all;




G = [9.23 15.2 11.43 12.1];
X = categorical({'RadIX','Telecom','nLight1','nLight2'});
X = reordercats(X,{'RadIX','Telecom','nLight1','nLight2'});

hB=bar(X,G,'BarWidth',0.5);          % generates one bar handle with five bars
hB.FaceColor='flat';              % ready to set color by bar w/ CData
hB.CData(1,:)=[0.8500 0.3250 0.0980];          % first and second
hB.CData(2,:)=[0 0.4470 0.7410];          % third red
hB.CData(3,:)=[0.4940 0.1840 0.5560];         % second green
hB.CData(4,:)=[0.4660 0.6740 0.1880];      % fifth white, leaves fourth default

% xtips2 = hB(1).XEndPoints;
% ytips2 = hB(1).YEndPoints;
% labels2 = {'Gs =1dB','Gs =5dB','Gs =13.14dB','Gs =1dB'};
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    % 'VerticalAlignment','bottom')

% title('Longueur optimale pour Ps = 40mw, Pp=100mw');
% ylabel('Lf(m)');


ylabel('FN(dB)');
title('FN(dB) en fonction des longueurs optimales');