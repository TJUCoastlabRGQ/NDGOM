function flux = matEvaluateNonconservativeFlux(obj, fm, fp, penaltyterm)
flux = (fm + fp)/2 + 1/2 * penaltyterm;
end