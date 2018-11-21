function deltaflux = matEvaluateDeltaSurfaceFlux( obj, mesh, Qx, Qy, fluxx, fluxy )
deltaflux = mesh.nx .* (Qx(mesh.eidM) - fluxx) + mesh.ny .* (Qy(mesh.eidM) - fluxy);
end