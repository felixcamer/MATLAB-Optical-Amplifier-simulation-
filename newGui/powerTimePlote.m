
%Load data here
tmp = strcat('.\res_ai\EDFres9');
load(tmp);
fs=10; % fontsize
cla (handles.axes4,'reset');
axes(handles.axes4);


Color='b';
Color2='r';
[AZ,H1,H2]=plotyy(time(1,1:jj)*1e9,PsL(1,1:jj)*1e3,time(1,1:jj)*1e9,PpL(1,1:jj)*1e3);
set(H1,'Color',Color,'LineWidth',1.5);
set(H2,'Color',Color2,'LineWidth',1.5);
set(AZ(1),'Fontsize',fs,'FontName','times','YColor',Color);
set(AZ(2),'Fontsize',fs,'FontName','times','YColor',Color2);
set(AZ(1),'box','off')
xlabel(AZ(1),'time (ns)','Fontsize',fs);
ylabel(AZ(1),'input signal power','Fontsize',fs);
set(AZ(2),'box','off')
ylabel(AZ(2) , 'Pump signal power','Fontsize',fs);
grid on;
