`spectralHMM`, version 1.0.0

### Table of Contents
* [SUMMARY](#SUMMARY)
* [LICENSES](#LICENSES)
* [REQUIREMENTS](#REQUIREMENTS)
* [USAGE](#USAGE)
* [EXAMPLES](#EXAMPLES)
* [BUILD](#BUILD)
* [EXTERNAL LIBRARIES](#EXTERNAL-LIBRARIES)
* [CONTACT](#CONTACT)
* [VERSION HISTORY](#VERSION-HISTORY)

#### SUMMARY

`spectralHMM` is a software program for computing likelihoods of time series allele frequencies data under given population genetic parameters (strength of selection in an arbitrary diploid model, mutation rates, population size, allele age). It assumes a model of a single di-allelic locus, with one wild type allele, and one derived allele that has the given selective advantage. These likelihoods can be used to infer the respective parameters. spectralHMM is an implementation of the method described in

Matthias Steinrücken, Anand Bhaskar, Yun S. Song (2014). [A novel spectral method for inferring general selection from time series genetic data](https://www.doi.org/10.1214/14-AOAS764). Ann. Appl. Stat. 8(4), pp. 2203-2222. [Preprint: [https://arxiv.org/abs/1310.1068v2](https://arxiv.org/abs/1310.1068v2)]


#### LICENSES

The source code is released under the GNU General Public License, version 3. The full text of the license can be found in the file `LICENSE`, which should have been included in this repository.

#### REQUIREMENTS

`SpectralHMM` is implemented in java, so a recent version of the java virtual machine (and possibly the java compiler) is required. The program was developed and tested under JDK/JSE 1.6. Furthermore, the scripts to build and run the program need python (any version greater than 2.5 should work).

When compiling the sourcecode, external libraries are required. For details, see section [EXTERNAL LIBRARIES](#EXTERNAL-LIBRARIES).

#### USAGE:

Download the content of this repository ([https://github.com/popgenmethods/spectralHMM](https://github.com/popgenmethods/spectralHMM)).

If the file `spectralHMM.jar` does not exist in the `main_dir` (the top level directory), building the sourcecode is necessary. Follow steps in section [BUILD](#BUILD) the compile the sourcecode.

In the `main_dir` (the top level directory), execute the python script to run the software with the command

```
python runSpectralHMM.py <arguments>
```

To print the usage and see the available command line arguments, execute the command

```
python runSpectralHMM.py --help
```

For some example calls, see the section [EXAMPLES](#EXAMPLES).

Note: Depending on the command line arguments, the program might require a substantial amount of memory. If the program fails due to insufficient memory, you can provide the `-Xmx` flag to have the JVM use more memory. For example:

```
python runSpectralHMM.py -Xmx10g <arguments>
```

will allow the program to use 10 GB.

#### EXAMPLES:

Note that all population genetic parameters are UNSCALED and the population size is given in terms of diploid individuals. Furthermore, all times area specified in years. To gauge the right precision and cutoff values it is helpful to run the analysis (for the extremal selection coefficients desired) several times with different values, to ensure a stable result.

1. The file `examples/single_withoutInitTime` contains a single temporal dataset and instructions on the data format. The command

```
python runSpectralHMM.py --inputFile examples/single_withOutInitTime --mutToBenef 1e-6 --mutFromBenef 1e-6 --effPopSize 10000 --yearsPerGen 5 --initFrequency 0.1 --initTime -20000 --hetF 0.000625 --homF 0.00125 --precision 40 --matrixCutoff 150 --maxM 140 --maxN 130
```

computes the likelihood of the data given (--inputFile) in the file under the following parameters: The per generation mutation probability from the wild type to the selected allele (--mutToBenef) is 1e-6, and the probability of the reverse event (--mutFromBenef) is also 1e-6. The (diploid) effective population size is 10000, and one generation corresponds to 5 years (--yearsPerGen). The initial frequency (--initFrequency 0.1) of derived advantageous alleles is 0.1 at time (--initTime) -20000 years. The fitness of a heterozygous (--hetF) individual is 0.000625 more than an individual homozygous for the wild type (reference fitness 0), and the fitness of an individual homozygous (--homF) for the derived allele is 0.00125. The computations are performed with a precision (--precision) of 10^(-40), and the size of the matrix whose eigenvalues yield the coefficients for the eigenfunctions (--matrixCutoff) is set to 150. Finally, 141 terms of the infinite sum that approximates the eigenfunctions (--maxM 140) are used, and 131 coefficients of the spectral expansion of the transition density (--maxN 130) are used.

Every line of the output that starts with a '#' denotes logging information. The relevant result is not preceded by a '#' and has 3 values on one line: The heterozygous fitness, the homozygous fitness, and the likelihood.

2. The file `examples/multi_withoutInitTime` contains several temporal dataset, one on each line of the input, and instructions on the data format. The command

```
python runSpectralHMM.py --multiplex --inputFile examples/multi_withOutInitTime --mutToBenef 1e-6 --mutFromBenef 1e-6 --effPopSize 10000 --yearsPerGen 5 --mutDriftBalance --initTime -20000 --hetF '[0.000625:0.0001:0.0008]' --homF '[0.00125:0.001:0.003]' --precision 40 --matrixCutoff 150 --maxM 140 --maxN 130 --condOnLastSegregating
```

computes the likelihood of the datasets given (--inputFile) in the file. In addition the --multiplex has to be specified to indicate multiple datasets in the input file. Most of the population genetic parameters, the precision and the cutoffs are specified as in the previous example. One difference is that instead of an initial frequency, it is specified that at time (--initTime) -20000 the allele frequency is drawn from the stationary distribution of the neutral model (--mutDriftBalance). Instead of a single parameter for the selective advantages (--hetF and --homF), now ranges are specified in the format [start:step:stop]. Finally, this command calculates the likelihood of the data conditional on the event that at least one derived allele is observed at the last sampling time point (--condOnLastSegregating).

Every line of the output that starts with a '#' denotes logging information. The relevant results are not preceded by a '#'. Each result line has 4 values: A running index for the dataset, the heterozygous fitness, the homozygous fitness, and the likelihood.

3. The file `examples/multi_withInitTime` contains several temporal dataset, one on each line of the input, and instructions on the data format. Here an initial time for each dataset is also given in the input file and does not have to be specified on the command line. The command

```
python runSpectralHMM.py --multiplex --inputFile examples/multi_withInitTime --mutToBenef 1e-6 --mutFromBenef 1e-6 --effPopSize 10000 --yearsPerGen 5 --initFrequency 0.1 --selection '[0.00125:0.001:0.003]' --dominance '[0.45:0.05:0.55]' --precision 40 --matrixCutoff 150 --maxM 140 --maxN 130 --condOnLastSegregating --ignoreBinomials
```

computes the likelihood of the datasets given (--inputFile) in the file. The parameters are mostly identical to the previous example. However, in this command, the diploid fitness values are given as a selection strength (--selection) and a dominance parameter (--dominance). Also, the likelihood is computed without the binomial factors (--ignoreBinomials) that would be required in the true probability model. Ignoring these factors does not change the shape of the likelihood surface.

Every line of the output that starts with a '#' denotes logging information. The relevant results are not preceded by a '#'. Each result line has 4 values: A running index for the dataset, the selection parameter, the dominance, and the likelihood.

#### BUILD

If the file `spectralHMM.jar` does not exist, Building from downloaded sourcecode is necessary.

Before the sourcecode can be compiled, the libraries detailed in the section [EXTERNAL LIBRARIES](#EXTERNAL-LIBRARIES) have to be downloaded.

Then in the `main_dir` (the top level directory) execute the command:

```
python build.py
```

to compile the sourcecode and create the jar-file: `spectralHMM.jar`

Once the java-executable `spectralHMM.jar` is created, see [USAGE](#USAGE) for instructions on how to build the program

#### EXTERNAL LIBRARIES

The following libraries (jar-files) are required when compiling the sourcecode:
- `JSAP-2.1.jar`.
	Download from http://sourceforge.net/projects/jsap/files/jsap/2.1/ (alt: http://www.martiansoftware.com/jsap/).
- `arpack_combined_all.jar`.
	Download from at http://en.sourceforge.jp/projects/sfnet_f2j/releases/ .
- `lapack_simple.jar`.
	Download `jlapack-0.8.tgz` from http://www.netlib.org/java/f2j/ and unpack. lapack_simple.jar can be found in this archive.

These jar files (or a symbolic link to the files) have to be put into a directory `main_dir/spectralHMM_lib/` (here `main_dir` denotes the top-level directory). If it is necessary to specify custom paths to these libraries, then the build instructions and the scripts `build.py` and `runSpectralHMM.py` have to be changed accordingly.

#### CONTACT

Please contact steinrue@uchicago.edu with bugs, comments, or questions regarding the software.


#### VERSION HISTORY

1.0.0: Initial release.
