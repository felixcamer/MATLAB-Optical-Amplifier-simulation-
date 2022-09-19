
function [X,Y1,Y2] = my_data4(root)

A = load(root);
X = (1e9)*A(:,1)';
Y1 = A(:,2)';
Y2 = A(:,3)';

end
