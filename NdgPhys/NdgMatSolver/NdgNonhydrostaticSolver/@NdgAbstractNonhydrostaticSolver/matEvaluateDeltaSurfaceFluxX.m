function deltaFluxX = matEvaluateDeltaSurfaceFluxX(obj, mesh, Qx, NumFluxx )
 deltaFluxX = mesh.nx .* ( Qx - NumFluxx );
end