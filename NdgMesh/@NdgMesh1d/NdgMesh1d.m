classdef NdgMesh1d < NdgMesh
    %LINE_MESH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        dim = enumMeshDim.One
    end
    
    methods(Hidden, Access = protected)
        function [ edge ] = makeConnectNdgEdge( obj, mesh1, mid0, mid1 )
            edge = NdgEdge1d( obj, mesh1, mid0, mid1 );
        end
        
        function obj = assembleJacobiFactor( obj )
            xr = obj.cell.Dr*obj.x;
            obj.J = xr; 
            obj.rx = 1./obj.J;
            
            obj.ry = zeros(size(obj.rx));
            obj.rz = zeros(size(obj.rx));
            obj.sx = zeros(size(obj.rx));
            obj.sy = zeros(size(obj.rx));
            obj.sz = zeros(size(obj.rx));
            obj.tx = zeros(size(obj.rx));
            obj.ty = zeros(size(obj.rx));
            obj.tz = zeros(size(obj.rx));
        end
        
        function [nx, ny, nz, Js] = assembleFacialJaobiFactor( obj )
            % Define outward normals
            xb = obj.x(obj.cell.Fmask, :);
            xc = obj.GetMeshAverageValue( obj.x );
            nx = zeros(2, numel(xc));
            nx(1,:) = sign( xb(1,:) - xc );
            nx(2,:) = sign( xb(2,:) - xc );
            
            Js = ones(size(nx));
            ny = zeros(obj.cell.Nface, obj.K);
            nz = zeros(obj.cell.Nface, obj.K);
        end
        
        function faceId = assembleGlobalFaceIndex( obj )
            faceId = zeros(obj.cell.Nface, obj.K);
            for f = 1:obj.cell.Nface
                faceId(f, :) = obj.EToV(obj.cell.FToV(1,f), :);
            end
        end% func
    end% methods
    
    methods
        %> \brief refine elements
        obj = refine(obj, refine_level);
        %> \brief draw mesh
        draw( obj, zvar );
        %> construction function
        function obj = NdgMesh1d( cell, Nv, vx, K, EToV, EToR, BCToV )
            vy = zeros(size(vx)); % vy is all zeros
            vz = zeros(size(vx)); % vz is all zeros
%             obj = obj@NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR, BCToV);
            obj = obj@NdgMesh(cell, Nv, vx, vy, vz, K, EToV, EToR);
        end% func
    end
    
end

