## Test environments
* local R installation (ubuntu 20.04), R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)
* R-hub windows-x86_64-devel (r-devel)
* R-hub ubuntu-gcc-release (r-release)
* R-hub fedora-clang-devel (r-devel)


## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Saint-Clair Chabert-Liddell <saint-clair.chabert-liddell@agroparistech.fr>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Barbillon (12:71)
    Chabert (12:54)
    Donnet (13:5)
    Lazega (13:16)
    Liddell (12:62)

0 errors ✓ | 0 warnings ✓ | 1 note x


* This is a new release.

* Note due to the name of the authors of the referenced article in the Description field of DESCRIPTION. 
