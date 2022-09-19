function [X,Y] = datacq1(root)

% acquisitions données IXfiberRad (part 1) pump wavelength

chemin=strcat(root,'.txt')

[col1,col2]= textread(chemin,'%f %f','delimiter',',','headerlines',1);
X=col1';        % longueur d'onde (nm)
Y=col2';        % section efficace absorption (m^2)

end

