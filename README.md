This repository contains programs associated with the following publication:

  Likelihood-based inference in isolation-by-distance models using the spatial distribution of low-frequency alleles.  
  Novembre J, Slatkin M (2009)  
  Evolution: International Journal of Organic Evolution 63(11):2914-25
  [PMID:19624728](http://www.ncbi.nlm.nih.gov/pubmed/19624728)

The software here is called RADDLE, which is
a reddish color, but also b/c it works with the acronym: Rare Allele
Distance of DispersaL Estimation.

## To install:

1. Make sure you have `gsl` installed.  
On a redhat system you might check by using "rpm -q gsl".  If it's not installed, you will have to install it.  
If you have apt set up (ubuntu, ect), the command would be "apt-get install libgsl0ldbl".  
On OS X, similar commands work with the port package installer ("sudo port install gsl")

2. Run the following commands:
```
cd src
make RaddleMC
```

For help:
```
bin/RaddleMC
```

## Example:

Here's an example usage that calculates only the MLE for each locus based on the
input file `Example_vx200_vx100_n9_L20.xy` with 10,000 importance sampling
replicates a driving value of 100 and a heat parameter of 2.  It also outputs a
two-dimensional likelihood surface from 80 to 250 with 5 grid-points in each
dimension:

```
cd example/
../bin/RaddleMC -S 10 20 -d 2 -n 3 -M 10000 -f Example_vx200_vx100_n9_L20 -v 100 -H 2 -s 80 250 5 1 -P
```

Note that the .xy needs to be dropped when the filename is passed to
the program (possibly annoying I know...).  Also, the file has 20
loci, but I've set "-n 3" in the example so that this run will only
analyze the first three loci.

The input format for 2-dimensional data begins with a simple header
line: "LocusID f nAlleles x y."  The rest of the data set is arranged
in rows.  For each locus, there is one row for each copy of the allele
observed.  The first three columns of each row are identical for each
locus.  They contain the LocusID, f (the sampling fraction at that
locus), and the number of alleles observed at the locus.  The last two
columns (if the data are two-dimensional) are the x and y positions of
that particular allele.  For 1-d data, the header line is "Locus ID f
nAlleles x" and the data rows lack the y coordinate.

The .MLEout file will have a number of lines that begin with # marks.
These lines output the parameters that were used for running the
program as well as results on the distribution of the importance
sampling weights.  The MLE results are output as table below the #
lines.  The idea is that the file can easily be read in as a table by
R, b/c R recognizes "#" lines as comments.

Each line of the results table is organized as:
Column1: LocusID or for the joint MLE the identifier "AllLoci" is used.
Column 2,3,4: The MLE for the constrained model (Sig_x^2=Sig_y^2),
followed by the lower and upper confidence interval boundaries.
Column 5,6,7: The MLE for Sig_x^2 followed by the lower and upper
confidence interval boundaries.
Column 8,9,10: The MLE for Sig_y^2 followed by the lower and upper
confidence interval boundaries.
Column 11,12: The log likelihood for the constrained model, followed
by the log of  the standard error on the estimate of the likelihood
Column 13,14: The log likelihood for the unconstrained model, followed
by the log of  the standard error on the estimate of the likelihood

The example above uses an option that is summarized in the help as:
" -s <double> <double> <int> <bool>       Min Sig2, Max Sig2,
nGridPoints, CalcMultiDimSurf"  This option is for doing grid-based
searches of the likelihood surface.  It's useful when the
hill-climbing algorithm used by default doesn't converge (in which
case an error message will be returned)  or if you want to ouput
actual likelihood surfaces for closer inspection.  The first three
parameters determine the range over which the grid is laid out and how
many points it has in each dimension.  The last parameter is a
boolean, such that if 1, then it will calculate the 2-d likelihood
surface for the unconstrained model, if 0, it will only do a 1-d curve
(the constrained model).  The output file is .MCout, which has a
header similar to MLEout, but the data table is such that each row has
a locusID, SigX gridpoint, SigY gridpoint, the sampling fraction, then
the LogL at that gridpoint, and it's corresponding standard error.

Note:
1.  If you do not use the -P option, the program will only calculate
and output the "AllLoci" MLE.  This results in faster performance if
you're not interested in MLEs from each locus.

2. If you do not use the -s option, the program will not calculate the
grid-based likelihood surface.  It will still output MLEs based on the
hill-climbing algorithm, so if you only care about the MLE and not the
likelihood surface, this results in much faster finishing times.

3.  The standard error for the "AllLoci" estimation is always -inf,
b/c I never finished implementing it (it's a little trickier to
calculate the SE of an estimate that is a product of values that are
individually with error.  And just for perspective, I don't place much
value on the standard errors calculated b/c they're only safe to
interpret as actual standard errors when the distribution of
importance sampling weights is not skewed (which is not necessarily
typical).
