function matVertNodeProjection( obj )
if obj.Nz == 1
    obj.VertToNode = zeros(size(obj.mesh2d.VertToNode,1),2*size(obj.mesh2d.VertToNode,2));
    %> surface
    for i = 1:obj.mesh2d.Nv
        obj.VertToNode(1,i) = obj.mesh2d.VertToNode(1,i);
        for j =2:obj.mesh2d.VertToNode(1,i)
            obj.VertToNode(j,i) = ( ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np) - 1 ) * obj.cell.Np + obj.cell.Np - obj.mesh2d.cell.Np + obj.mesh2d.VertToNode(j,i) - ( ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np) - 1 ) * obj.mesh2d.cell.Np;
        end
    end
    %> bottom
    for i = obj.mesh2d.Nv + 1 : 2 * mesh2d.Nv
        obj.VertToNode(1,i) = obj.mesh2d.VertToNode(1,i-obj.mesh2d.Nv);
        for j = 2:obj.mesh2d.VertToNode(1,i-obj.mesh2d.Nv)
            obj.VertToNode(j,i) = (ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np) - 1)*obj.cell.Np + obj.mesh2d.VertToNode(j,i-obj.mesh2d.Nv);
        end
    end
    obj.NodeWeight = horzcat(obj.mesh2d.NodeWeight, obj.mesh2d.NodeWeight);
else
    %> The first row need to be excluded
    obj.VertToNode = zeros(2*size(obj.mesh2d.VertToNode,1)-1,obj.Nv);
    %> Data number need to be doubled
    obj.NodeWeight = zeros(2*size(obj.mesh2d.NodeWeight,1),obj.Nv);
    %> The surface
    for i = 1:obj.mesh2d.Nv
        obj.VertToNode(1,i) = obj.mesh2d.VertToNode(1,i);
        for j=2:obj.mesh2d.VertToNode(1,i)+1
            obj.VertToNode(j,i) = (ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np)-1)*obj.cell.Np*obj.Nz + obj.cell.Np - obj.mesh2d.cell.Np + obj.mesh2d.VertToNode(j,i) -  ( ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np) - 1 ) * obj.mesh2d.cell.Np;
            obj.NodeWeight(j-1,i) = obj.mesh2d.NodeWeight(j-1,i);
        end
    end
    %> The middle layers
    for L=1:obj.Nz-1
        for i = 1:obj.mesh2d.Nv
            obj.VertToNode(1,L*obj.mesh2d.Nv+i) = 2*obj.mesh2d.VertToNode(1,i);
            h1 = obj.vz((L-1)*obj.mesh2d.Nv+i)-obj.vz(L*obj.mesh2d.Nv+i);
            h2 = obj.vz(L*obj.mesh2d.Nv+i)-obj.vz((L+1)*obj.mesh2d.Nv+i);
            %> The upper part
            for j = 2:obj.mesh2d.VertToNode(1,i)+1
                obj.VertToNode(j,i + L*obj.mesh2d.Nv) = (ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np)-1)*obj.cell.Np*obj.Nz + (L-1)*obj.cell.Np + obj.mesh2d.VertToNode(j,i) -  ( ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np) - 1 ) * obj.mesh2d.cell.Np;
                obj.NodeWeight(j-1,i+L*obj.mesh2d.Nv) = h1/(h1+h2)*obj.mesh2d.NodeWeight(j-1,i);
            end
            %> The lower part
            for j = obj.mesh2d.VertToNode(1,i)+2:2*obj.mesh2d.VertToNode(1,i)+1
                obj.VertToNode(j,i+L*obj.mesh2d.Nv) = (ceil(obj.mesh2d.VertToNode(j- obj.mesh2d.VertToNode(1,i),i)/obj.mesh2d.cell.Np)-1)*obj.cell.Np*obj.Nz + L*obj.cell.Np + obj.cell.Np - obj.mesh2d.cell.Np + obj.mesh2d.VertToNode(j - obj.mesh2d.VertToNode(1,i),i) - ( ceil(obj.mesh2d.VertToNode(j - obj.mesh2d.VertToNode(1,i),i)/obj.mesh2d.cell.Np) - 1 ) * obj.mesh2d.cell.Np; ;
%                 mod(obj.mesh2d.VertToNode(j - obj.mesh2d.VertToNode(1,i),i),obj.mesh2d.cell.Np) + floor(obj.mesh2d.VertToNode(j - obj.mesh2d.VertToNode(1,i),i)/obj.mesh2d.cell.Np) * obj.mesh2d.cell.Np
                obj.NodeWeight(j-1,i+L*obj.mesh2d.Nv) = h2/(h1+h2)*obj.mesh2d.NodeWeight(j-obj.mesh2d.VertToNode(1,i)-1,i);
            end
        end
    end
    %> The bottom layer
    for i = 1:obj.mesh2d.Nv
        obj.VertToNode(1,obj.Nz*obj.mesh2d.Nv+i) = obj.mesh2d.VertToNode(1,i);
        for j = 2:obj.mesh2d.VertToNode(1,i) + 1
            obj.VertToNode(j,i+obj.Nz*obj.mesh2d.Nv) = (ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np)-1)*obj.cell.Np*obj.Nz + (obj.Nz-1)*obj.cell.Np + obj.mesh2d.VertToNode(j,i) -  ( ceil(obj.mesh2d.VertToNode(j,i)/obj.mesh2d.cell.Np) - 1 ) * obj.mesh2d.cell.Np;
            obj.NodeWeight(j-1,i+obj.Nz*obj.mesh2d.Nv) = obj.mesh2d.NodeWeight(j-1,i);
        end
    end
end
end