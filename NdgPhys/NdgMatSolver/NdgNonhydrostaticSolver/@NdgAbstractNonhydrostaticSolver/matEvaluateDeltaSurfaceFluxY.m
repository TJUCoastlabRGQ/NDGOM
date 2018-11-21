function deltaFluxY = matEvaluateDeltaSurfaceFluxY(obj, mesh, Qy, NumFluxy )
 deltaFluxY = mesh.ny .* ( Qy - NumFluxy );
end