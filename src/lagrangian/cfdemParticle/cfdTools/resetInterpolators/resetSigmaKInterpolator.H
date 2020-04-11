sigmaKInterpolator_.reset(interpolation<scalar>::New(propsDict_.lookupOrDefault("sigmaKInterpolationType",word("cellPointFace")),sigmaK_).ptr());
