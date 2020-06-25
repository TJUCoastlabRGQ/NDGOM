function [ Adv, HorizVis, VertVis ] = initSolver( obj )
if obj.option.isKey('AdvDiffHorizontalDiffusionType')
    type = obj.getOption('AdvDiffHorizontalDiffusionType');
    if type == enumHorizontalDiffusion.Constant
        HorizVis = NdgHorizDiffSolver( obj );
    else
        HorizVis = NdgNoneHorizDiffSolver( obj );
    end
else
    HorizVis = NdgNoneHorizDiffSolver( obj );
end

if obj.option.isKey('AdvDiffVerticalDiffusionType')
    type = obj.getOption('AdvDiffVerticalDiffusionType');
    if type == enumVerticalDiffusion.Constant
        VertVis = NdgVertDiffSolver( obj );
    else
        VertVis = NdgNoneVertDiffSolver( obj );
    end
else
    VertVis = NdgNoneVertDiffSolver( obj );
end

integralType = obj.getOption('integralType');
equType = obj.getOption('equationType');

if (integralType == enumDiscreteIntegral.QuadratureFree)
    if( equType == enumDiscreteEquation.Strong )
        Adv = NdgQuadFreeStrongFormAdvSolver3d( obj );
    elseif( equType == enumDiscreteEquation.Weak )
        Adv = NdgQuadFreeWeakFormAdvSolver3d( obj );
    end
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    if( equType == enumDiscreteEquation.Strong )
        Adv = NdgGaussQuadStrongFormAdvSolver3d( obj );
    elseif( equType == enumDiscreteEquation.Weak )
        Adv = NdgGaussQuadWeakFormAdvSolver3d( obj );
    end
end
end