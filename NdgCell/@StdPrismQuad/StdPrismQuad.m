classdef StdPrismQuad < handle
    
    properties ( Constant )
        %> reference element type
        type = enumStdCell.PrismQuad;
        %> num of vertices
        Nv = 8
        %> volume of reference elements
        LAV = 8
        %> vertices coordinates
        vr = [-1,  1,  1, -1, -1,  1,  1, -1]'
        vs = [-1, -1,  1,  1, -1, -1,  1,  1]'
        vt = [-1, -1, -1, -1,  1,  1,  1,  1]'
        %> num of vertices on each face
        Nfv = [4, 4, 4, 4, 4, 4];
        %> index of interpolation points on each face
        FToV = [1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8; 1,2,3,4; 5,6,7,8]';
        %> num of faces
        Nface = 6
        %> element types of each face
        faceType = [enumStdCell.Quad, ...
                    enumStdCell.Quad, ...
                    enumStdCell.Quad, ...
                    enumStdCell.Quad, ...
                    enumStdCell.Quad, ...
                    enumStdCell.Quad ]
    end
    
    properties
        %> maximum order of basis function
        N
        %> maximum polynomial order in vertical direction
        Nz
    end

    properties ( SetAccess = protected )
        %> num of interpolation points (IP)
        Np, Nph, Npz
        %> coordinates of interpolation points
        r, s, t
        %> horizontal and vertical node coordinate
        r1, s1, t1        
        %> node index of facial points
        Fmask
        %> number of facial interpolation points
        Nfp
        %> total number of face points
        TNfp
    end
    properties ( SetAccess = protected)
        %> Vandermonde matrix
        V, Vh,Vint
        %> project matrx from interpolation points to central point in
        %> vertical direction, and this is needed for GOTM implementation
        VCV        
        %> project matrx from interpolation points to quadrature points
        Vq
        %> mass matrix
        M
        %> inverse of mass matrix
        invM
        %> derivative matrix with
        %> \f$ [Dr]_{ij} = \left.\frac{\partial l_j}{\partial r}\right|_{r_i} \f$
        Dr, Ds, Dt
    end
    
    properties ( SetAccess = protected )
        %> number of gauss quadrature points
        Nq
        %> coordinate of quadrature points
        rq, sq, tq
        %> integral weights for each quadrature points
        wq
    end
    
    methods
        % construction function of triangular prism reference element
        function obj = StdPrismQuad(Nh, Nz)
            obj.N = Nh; 
            obj.Nz = Nz;
            EvaluaetNodeCoor( obj, Nh, Nz );
            
            [ obj.Nq, obj.rq, obj.sq, obj.tq, obj.wq ] ...
                = obj.quad_coor_func( Nh, Nz );
            
            AssembleVandMatrix( obj );
            
            [ obj.Vq ] = obj.assembleQuadratureMatrix( );
            [ obj.M, obj.invM ] = obj.assembleMassMatrix( );
            [ obj.Dr, obj.Ds, obj.Dt ] ...
                = obj.nodal_derivative_func(obj.r, obj.s, obj.t);
% 
            obj.Nfp = zeros(obj.Nface, 1);
            for i = 1:4
                obj.Nfp(i) = ( Nh + 1 ) * ( Nz +1 );
            end
            obj.Nfp( [ 5, 6 ] ) = obj.Nph;
            
            [ obj.TNfp ] = sum(obj.Nfp);
            [ obj.Fmask ] = obj.assembleFacialNodeIndex();
        end
        
        [ fun ] = orthogonal_func(obj, N1, N2, ind, r, s, t);
        [ node_val ] = project_vert2node(obj, vert_val);
        [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J, Jz ] = assembleJacobianMatrix( obj, x, y, z );
        [ nx, ny, nz, Js ] = assembleNormalVector( obj, x, y, z );
        
        %> @brief Evaluate all the nodal basis function values at points
        %> @param[in] obj The StdCell class
        %> @param[in] r The node coordinate
        %> @param[in] s The node coordinate
        %> @param[in] t The node coordinate
        %> @param[out] func The basis function values at points
        function [ func ] = nodal_func(obj, r, s, t)
            func = zeros(numel(r), obj.Np);
            for n = 1:obj.Np
                fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, r, s );
                fv = obj.EvaluateVerticalOrthogonalFunc( n, t );
                func(:, n) = fh .* fv;
            end
            func = func/obj.V;
        end
        
        %> @brief Project the scalar field from the interpolation nodes to the Gauss quadrature nodes
        %> @param[in] obj The StdCell class
        %> @param[in] node_val The values on these interpolation nodes
        function quad_val = project_node2quad(obj, node_val)
            quad_val = obj.Vq * node_val;
        end% func        

        %> @brief Evaluate the derivative nodal function values at points
        function [ fDr, fDs, fDt ] = nodal_derivative_func( obj, r, s, t )
            Nr = numel( r );
            Vr = zeros(Nr, obj.Np);
            Vs = zeros(Nr, obj.Np);
            Vt = zeros(Nr, obj.Np);
            for n = 1:obj.Np
                [Vr(:, n), Vs(:, n), Vt(:, n)] = ...
                    obj.derivative_orthogonal_func(obj.N, obj.Nz, n, r, s, t);
            end
            fDr = Vr/obj.V;
            fDs = Vs/obj.V;
            fDt = Vt/obj.V;
        end
    end% methods

    methods ( Access=protected )
        [ Np, Nph, Npz, r, s, t ] = node_coor_func( obj, Nh, Nv );
        [ dr, ds, dt ] = derivative_orthogonal_func( obj, Nh, Nv, ind, r, s, t );
        [ Nq, rq, sq, tq, wq ] = quad_coor_func( obj, Nh, Nv );

        %> @brief Assemble the interpolation matrix of Gauss quadrature nodes
        %> The elements of the quadratuer interpolation matrix is
        %> \f$ [V_q]_{i,j} = l_j(\xi_i) \f$
        %> where \f$ \xi_i \f$ is the ith Gauss quadrature nodes.
        function [ Vq ] = assembleQuadratureMatrix( obj )
            Vq = obj.nodal_func( obj.rq, obj.sq, obj.tq );
        end

        %> @brief Assemble the Vandermonde matrix
        %> @details The Vandermonde matrix interpolate the orthgonal basis functions
        %> to the nodal basis functions, with
        %> \f$ [V]_{i,j} = P_j(r_i) \f$
        %> where \f$ r_i \f$ is the ith interpolation nodes and \f$ P_j \f$
        %> is the jth orthgonal function.
        function V = assembleVandMatrix(obj, orthogonal_func)
            V = zeros(obj.Np, obj.Np);
            for n = 1:obj.Np
                V(:, n) = orthogonal_func(obj.N, obj.Nz, n, obj.r, obj.s, obj.t);
            end% for
        end% func

        %> @brief Assemble the mass matrix
        function [ M, invM ] = assembleMassMatrix( obj )
            invV = inv(obj.V);
            M = (invV')*invV;
            invM = obj.V * obj.V';
        end% func

        function Fmask = assembleFacialNodeIndex(obj)
            maxnfp = max(obj.Nfp);
            Fmask = zeros(maxnfp, obj.Nface);
            
            % vertical faces
            f = 1;
            ind = find( abs(obj.s + 1) < 1e-10 );
            Fmask(1:obj.Nfp(f), f) = ind;
            
            f = 2;
            ind = find( abs(obj.r - 1) < 1e-10 );
            Fmask(1:obj.Nfp(f), f) = ind;
            
            f = 3;
            ind = find( abs(obj.s - 1) < 1e-10 );
            Fmask(1:obj.Nfp(f), f) = ind;
            % horizontal faces
            f = 4;
            ind = find( abs(obj.r + 1) < 1e-10 );
            Fmask(1:obj.Nfp(f), f) = ind;
            
            f = 5;
            ind = find( abs(obj.t + 1) < 1e-10 );
            Fmask(1:obj.Nfp(f), f) = ind;
            
            f = 6;
            ind = find( abs(obj.t - 1) < 1e-10 );
            Fmask(1:obj.Nfp(f), f) = ind;            
            
        end% func
        
    end
    
end

