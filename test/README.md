# Testing the examples

In order to test a new version of ComplexMixtures against these examples, 
start `julia` in the environment of this package, adjust the version of the 
ComplexMixtures package desired (for instance `dev ComplexMixtures`), and
run the package tests.

Ideally, to test the generation of all images, download the full trajectory
files and place them at `test/trajectories`. The files are:

(200Mb): https://drive.google.com/file/d/1BuXJ8AjBeduMSD2CkDJLDNxAAD2QNNg6/view?usp=sharing

(280Mb): https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing

(365Mb): https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing

(3Gb): https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing

The tests will generate all figures again, which can be compared to the old versions
to see if everything seems fine for the new version of the package. Checking the new
version of the images against previous ones is convenient in VSCode or in a github 
pull request, which show the images side by side. 









