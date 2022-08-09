# subformula_graph
A program which reads in mass spectral data (currently netCDF) and computes the most likely set of formulae corresponding to the peaks

# How to run this program
You're going to need to install docker and make first.
Then 'make build' will build the image - it'll take a while at first, but the conda dependencies are cached
so it'll be quicker next time you build. Then 'make run' will actually run it.

