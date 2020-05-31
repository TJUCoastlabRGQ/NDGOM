classdef AbstractDiffSolver < handle
    %ABSTRACTDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Access = protected)
        %> the diffusion coefficient
        nv
        %> the penalty parameter for IP form
        tau
        %> the Prantal number for heat and sault diffusion
        Prantl
    end
    
    methods( Access = protected, Abstract)
        %> this function is used to update the viscosity
        matUpdateViscosity(obj)
        %> this function is used to update the penalty parameter adopted in
        %> IP form
        matUpdatePenaltyParameter(obj)
        
    end
    
    methods(Access = protected)
        function [ fm, fp ] = matEvaluateSurfValue(obj, edge, Kappa )
            [ fm, fp ] = mxEvaluateSurfValue( edge.FToE, edge.FToN1, edge.FToN2, Kappa );
        end
        
        function [ frhs ] = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS, fm, fp )
            frhs = mxEvaluateStrongFromEdgeRHS( ...
                obj.mesh.cell.invM, obj.M, ...
                obj.FToE, obj.FToN1, obj.FToN2, ...
                obj.Js, obj.mesh.J, fluxM, fluxP, fluxS, fm, fp );
        end
        
    end
    
end

