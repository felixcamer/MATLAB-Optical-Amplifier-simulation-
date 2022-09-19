function [X,Y,Z] = datacq2(root)

% acquisitions données IXfiberRad (part 2) C band

chemin=strcat(root,'.txt')

[col1, col2,col3]= textread(chemin,'%f %f %f','delimiter',',','headerlines',1);
X=col1';        % longueur d'onde (nm)
Y=col2';        % section efficace absorption (m^2)
Z=col3';        % section efficace émission (m^2)

end

