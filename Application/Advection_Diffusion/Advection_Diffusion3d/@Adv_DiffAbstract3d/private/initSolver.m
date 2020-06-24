function [ Adv, HorizVis, VertVis ] = initSolver( obj )
if obj.option.isKey('AdvDiffHorizontalEddyViscosityType')
    type = obj.getOption('AdvDiffHorizontalEddyViscosityType');
    if type == enumHorizontalEddyViscosity.Constant
        HorizVis = NdgHorizDiffSolver( obj );
    else
        HorizVis = NdgNoneHorizDiffSolver( obj );
    end
else
    HorizVis = NdgNoneHorizDiffSolver( obj );
end

if obj.option.isKey('AdvDiffVerticalEddyViscosityType')
    type = obj.getOption('AdvDiffVerticalEddyViscosityType');
    if type == enumVerticalEddyViscosity.Constant
        VertVis = NdgVertDiffSolver( obj );
    else
        VertVis = NdgNoneVertDiffSolver( obj );
    end
else
    VertVis = NdgNoneVertDiffSolver( obj );
end

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