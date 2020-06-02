# sci

This repository contains libraries for scientific computing.

  - `flomat`: floating point matrices.

# Installation

Use 

    raco pkg install sci
	
to install the Racket part of `sci`.

For macOS this is enough.

For Linux/Windows check details in the documentation.

The library `flomat` relies on CBLAS and LAPACK, so they need to be installed 
in a place where they can be found by Racket.


# License

We follow the lead of the main Racket license:
Sci is distributed under the MIT license and the Apache version 2.0 license, at your option. 
