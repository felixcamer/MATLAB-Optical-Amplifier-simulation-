function [X,Y] = datacq3(root)

% acquisitions données IXfiberTelecom (part 1); pump wavelength

[num,txt,raw] = xlsread( fullfile(root, 'iXBlueTelpart1.xlsx') );

for ii=1:size(txt,1)
    B=cell2mat(txt(ii,1));
    X(1,ii)=(1e9)*str2double(B);  % longueur d'onde (nm)

    C=cell2mat(txt(ii,2));
    Y(1,ii)=str2double(C);  % section efficace absorption (m^2)
end



end

