classdef NdgFlatBottomNonhydrostaticTest < NdgNonhydrostaticAbstractTest
    %NDGFLATBOTTOMNONHYDROSTATICTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgFlatBottomNonhydrostaticTest(N, cellType)
            obj = obj@NdgNonhydrostaticAbstractTest(N, cellType);
        end
        
        function [Exactfm, Exactfp] = getExactNonhydroSurfaceValue( obj )
            Num = size(obj.meshUnion(1).eidM, 1);Nf = obj.meshUnion(1).cell.Nface;
            Exactfm = ones(size(obj.meshUnion(1).eidM)); Exactfp = ones(size(obj.meshUnion(1).eidM));
            Exactfp([1:Num/Nf,2*Num/Nf+1:3*Num/Nf,3*Num/Nf+1:4*Num/Nf], 1) = - Exactfp([1:Num/Nf,2*Num/Nf+1:3*Num/Nf,3*Num/Nf+1:4*Num/Nf], 1);
            Exactfp([1:Num/Nf,Num/Nf+1:2*Num/Nf,2*Num/Nf+1:3*Num/Nf], 2) = - Exactfp([1:Num/Nf,Num/Nf+1:2*Num/Nf,2*Num/Nf+1:3*Num/Nf], 2);
        end
        
        function [ExactPenalX, ExactPenalY] = getExactPenaltyTerm(obj)
            Num = size(obj.meshUnion(1).eidM, 1);Nf = obj.meshUnion(1).cell.Nface;
            ExactPenalX = zeros(size(obj.meshUnion(1).eidM)); ExactPenalY = zeros(size(obj.meshUnion(1).eidM)); 
            ExactPenalX(3*Num/Nf+1:4*Num/Nf ,1 ) = -2; ExactPenalX(Num/Nf+1:2*Num/Nf ,2 ) = 2;
            ExactPenalY(1:Num/Nf,:) = -2;ExactPenalY(2*Num/Nf+1:3*Num/Nf,:) = 2;
        end
        
        function ExactharmonicTerm = getExactHarmonicTerm(obj)
            Num = size(obj.meshUnion(1).eidM, 1);Nf = obj.meshUnion(1).cell.Nface;
            ExactharmonicTerm = zeros(size(obj.meshUnion(1).eidM));
            ExactharmonicTerm(Num/Nf+1:2*Num/Nf ,1 ) = 1; ExactharmonicTerm(3*Num/Nf+1:4*Num/Nf ,2 ) = 1;
        end
        
        function [ExactNumfluxX, ExactNumfluxY] = getExactNonconservativeFlux(obj)
             [ExactPenalX, ExactPenalY] = obj.getExactPenaltyTerm;
             ExactharmonicTerm = obj.getExactHarmonicTerm;
             [Exactfm, Exactfp] = obj.getExactNonhydroSurfaceValue;
             ExactNumfluxX = ExactharmonicTerm .* (Exactfm./2 + Exactfp./2 + ExactPenalX./2);
             ExactNumfluxY = ExactharmonicTerm .* (Exactfm./2 + Exactfp./2 + ExactPenalY./2);
        end
        
        function  Deltaflux = getExactDeltaFlux(obj)
             [ExactPenalX, ExactPenalY] = obj.getExactPenaltyTerm;
             ExactharmonicTerm = obj.getExactHarmonicTerm;
             [Exactfm, Exactfp] = obj.getExactNonhydroSurfaceValue;
             ExactNumfluxX = ExactharmonicTerm .* (Exactfm./2 + Exactfp./2 + ExactPenalX./2);
             ExactNumfluxY = ExactharmonicTerm .* (Exactfm./2 + Exactfp./2 + ExactPenalY./2); 
             Deltaflux = obj.meshUnion(1).nx .* ( Exactfm - ExactNumfluxX ) + obj.meshUnion(1).ny .* ( Exactfm - ExactNumfluxY ) ;
        end
        
        function ExactDeltafluxY = getExactDeltaFluxY(obj)
            Num = size(obj.meshUnion(1).eidM, 1);Nf = obj.meshUnion(1).cell.Nface;
            ExactDeltafluxY = zeros(size(obj.meshUnion(1).eidM));    
            ExactDeltafluxY(1:Num/Nf ,: ) = -1;
            ExactDeltafluxY(2*Num/Nf+1:3*Num/Nf ,: ) = 1;
        end
        
        function ExactDeltafluxX = getExactDeltaFluxX(obj)
            Num = size(obj.meshUnion(1).eidM, 1);Nf = obj.meshUnion(1).cell.Nface;
            ExactDeltafluxX = zeros(size(obj.meshUnion(1).eidM));    
            ExactDeltafluxX(3*Num/Nf + 1 : 4*Num/Nf ,1 ) = -1;
            ExactDeltafluxX(Num/Nf+1:2*Num/Nf ,2 ) = 1;
        end   
        
        function ExactDivergenceX = getExactDivergenceX(obj)
            Q = ones(size(obj.meshUnion(1).x));
            ExactDivergenceX = obj.meshUnion(1).rx .* ( obj.meshUnion(1).cell.Dr * Q )+...
                obj.meshUnion(1).sx .* ( obj.meshUnion(1).cell.Ds * Q );
        end
        
        function ExactDivergenceY = getExactDivergenceY(obj)
            Q = ones(size(obj.meshUnion(1).y));
            ExactDivergenceY = obj.meshUnion(1).ry .* ( obj.meshUnion(1).cell.Dr * Q )+...
                obj.meshUnion(1).sy .* ( obj.meshUnion(1).cell.Ds * Q );
        end 
        
        
    end
    
   methods(Access = protected)
       
        function fphys = setInitialField( obj )
            fphys = cell(obj.Nmesh, 1);
            fphys{1}(:,:,4) = ones(size(obj.meshUnion(1).x));
        end
   end
    
end

