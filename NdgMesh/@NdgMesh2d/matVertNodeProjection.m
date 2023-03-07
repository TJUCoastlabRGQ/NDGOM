function matVertNodeProjection( obj )
VertEleNum = ElementToVertexRelation( obj.EToV, obj.Nv );
% obj.VertToNode = zeros(max(VertEleNum)+1, obj.Nv);
% obj.NodeWeight = zeros(max(VertEleNum), obj.Nv);
VToK = zeros(max(VertEleNum),obj.Nv);
Nvc = zeros(obj.Nv,1);
for k = 1:obj.K
    v = obj.EToV(:,k);
    ind = Nvc(v) + 1 + (v-1)*max(VertEleNum);
    VToK(ind) = k;
    Nvc(v) = Nvc(v) + 1;
end
obj.VertToNode = zeros(max(Nvc)+1, obj.Nv);
obj.NodeWeight = zeros(max(Nvc), obj.Nv);
obj.VertToNode(1,:) = Nvc;
for i = 1:obj.Nv
    x = obj.vx(i); y = obj.vy(i);
    for j = 1:Nvc(i)
        ele = VToK((i-1)*max(Nvc)+j);
        for p = 1:obj.cell.Np
            if abs(obj.x(p,ele)-x)<=1.0e-14 && abs(obj.y(p,ele)-y)<=1.0e-14
                obj.VertToNode(j+1,i) = p + (ele-1)*obj.cell.Np;
                break;
            end
        end
    end
end
for k = 1:obj.K
    Area = PartialAreaCalculation( obj.vx(obj.EToV(:,k)), obj.vy(obj.EToV(:,k)), obj.cell.Nface );
    for vert = 1:obj.cell.Nface
        for j = 1:max(Nvc)
            if VToK(j, obj.EToV(vert,k)) == k
                obj.NodeWeight(j, obj.EToV(vert, k)) = Area(vert);
                break;
            end
        end
    end
end
obj.NodeWeight = obj.NodeWeight ./ sum(obj.NodeWeight);
end


%> This function is used to calculate the area, the base point is included by default.
%> x: The x coordinate of line or polygon
%> y: The y coordinate of line or polygon
function A = AreaCalculation( x, y, Num)
A = 0.0;
if Num == 2
    A = A + x(1)*y(2) - x(2)*y(1);
else
    for i = 1:Num-1
        x1 = x(i); y1 = y(i);
        x2 = x(i+1); y2 = y(i+1);
        A = A + ( x1 * y2 - x2 * y1 );
    end
    x1 = x(Num); y1 = y(Num);
    x2 = x(1); y2 = y(1);
    A = A + ( x1 * y2 - x2 * y1 );
end
end

%> This function is used to calculate the coordinates of the input polygon
%> x: The x coordinate of a polygon
%> y: The y coordinate of a polygon
%> Num: The number of points of the studied polygon
function [xc, yc] = CentroidCalculation( x, y, Num )
WeightAx = 0.0; WeightAy = 0.0;
for i = 1:Num-1
    %     x1 = x(i); y1 = y(i);
    %     x2 = x(i+1); y2 = y(i+1);
    WeightAx = WeightAx + AreaCalculation(x([i,i+1]),y([i,i+1]),2)*(x(i) + x(i+1))/3.0;
    WeightAy = WeightAy + AreaCalculation(x([i,i+1]),y([i,i+1]),2)*(y(i) + y(i+1))/3.0;
end
WeightAx = WeightAx + AreaCalculation(x([Num,1]),y([Num,1]),2)*(x(Num) + x(1))/3.0;
WeightAy = WeightAy + AreaCalculation(x([Num,1]),y([Num,1]),2)*(y(Num) + y(1))/3.0;
xc = WeightAx/AreaCalculation( x, y, Num );
yc = WeightAy/AreaCalculation( x, y, Num );
end

%> This function is used to calculate the area of the polygon with the base point as one of it vertex
%> x: The x coordinate of a polygon
%> y: The y coordinate of a polygon
%> Num: The number of points of the studied polygon
function Area = PartialAreaCalculation( x, y, Num )
[xc, yc] = CentroidCalculation( x, y, Num );
Lxc = zeros(Num, 1); Lyc = zeros(Num, 1);
Area = zeros(Num,1);
for i = 1:Num-1
    Lxc(i) = (x(i) + x(i+1))/2;
    Lyc(i) = (y(i) + y(i+1))/2;
end
Lxc(Num) = (x(1) + x(Num))/2;
Lyc(Num) = (y(1) + y(Num))/2;
Area(1) = AreaCalculation([x(1), Lxc(1), xc, Lxc(Num)], [y(1), Lyc(1), yc, Lyc(Num)], 4);
for i = 2:Num
    Area(i) = AreaCalculation([x(i), Lxc(i), xc, Lxc(i-1)], [y(i), Lyc(i), yc, Lyc(i-1)], 4);
end
end

%> This function is used to calculate the number of elements that each vertex belongs to
%> EToV: Element to vertex relation
%> Nv: Number of vertex
function VertEleNum = ElementToVertexRelation( EToV, Nv )
VertEleNum = zeros( 1, Nv );
for i = 1:numel(EToV)
    VertEleNum(EToV(i)) = VertEleNum(EToV(i)) + 1;
end
end