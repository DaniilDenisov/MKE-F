function [points, displacement] = nafems_challenge5_shapes(model, modes, count)
%NAFEMS_CHALLENGE5_SHAPES Common physical grid, exact beam interpolation.
% Output rows: x,y at 201 points per member, shared joint retained per member.
if nargin < 3, count = 201; end
n = numel(model.elementData)/2;
points = zeros(2*count,2); displacement = zeros(4*count,size(modes,2));
for member = 1:2
    for k = 1:count
        t = (k-1)/(count-1); index = min(floor(t*n)+1,n);
        xi = t*n-(index-1);
        e = model.elementData((member-1)*n+index); L=e.length;
        local = e.transformation*modes(e.dofs,:);
        axial = (1-xi)*local(1,:)+xi*local(4,:);
        transverse = (1-3*xi^2+2*xi^3)*local(2,:) + L*(xi-2*xi^2+xi^3)*local(3,:) + ...
            (3*xi^2-2*xi^3)*local(5,:) + L*(-xi^2+xi^3)*local(6,:);
        row = (member-1)*count+k;
        points(row,:) = (1-xi)*e.nodeCoordinates(1,:)+xi*e.nodeCoordinates(2,:);
        displacement(2*row-1:2*row,:) = e.transformation(1:2,1:2).'*[axial;transverse];
    end
end
end
