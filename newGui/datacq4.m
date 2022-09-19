function [X,Y1,Y2] = datacq4(root)

% acquisitions données IXfiberTelecom (part 2), bande C

[num,txt,raw] = xlsread( fullfile(root, 'iXBlueTelpart2.xlsx') );

for ii=1:size(txt,1)
    B=cell2mat(txt(ii,1));
    X(1,ii)=(1e9)*str2double(B);% longueur d'onde (nm)

    C=cell2mat(txt(ii,2));
    Y1(1,ii)=str2double(C);% section efficace absorption (m^2)
    
    D=cell2mat(txt(ii,3));
    Y2(1,ii)=str2double(D); % section efficace émission (m^2)
end




end

