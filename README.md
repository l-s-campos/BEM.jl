# BEM

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://l-s-campos.github.io/BEM.jl/stable)

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> BEM

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
2. Open a Julia console and do:

   ```Julia
   using Pkg

   Pkg.activate(pwd())

   Pkg.instantiate()

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.
Then run a script from the scripts folder to run the simulation
