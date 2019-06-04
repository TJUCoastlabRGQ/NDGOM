function [NonhydroFmPoint, NonhydroFpPoint, WetDryFaceOrder] = matAssemblePointRelatedInformation(obj, ZeroFluxBoundary, AdjacentDryCellAndFace, FToE, FToN1)
[NonhydroFmPoint, NonhydroFpPoint, WetDryFaceOrder] = mxAssemblePointRelatedInformation(ZeroFluxBoundary, AdjacentDryCellAndFace, FToE, FToN1);
end