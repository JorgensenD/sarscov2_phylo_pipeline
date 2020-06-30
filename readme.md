# Step by step guide to analyzing a regional data set with the SEIJR model

version: 2.0.0

David Jorgensen

---
**Code contributed by: Olivia Boyd, Lily Geidelberg, David Jorgensen, Manon Ragonnet, Igor Siveroni and Erik Volz**

**Affiliations:**

[sarscov2phylodynamics.org](http://sarscov2phylodynamics.org/)

[MRC Centre for Global Infectious Disease Analysis at Imperial College London](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis).

---

This walk through uses the example of Madinah, Saudi Arabia and should allow a user to reproduce updated versions of the figures presented in the [blog post on sarscov2phylodynamics.org](http://sarscov2phylodynamics.org/2020/06/07/Madinah-April-20.html). Due to restrictions on sharing gisaid data the sequence and metadata used are not provided in full alongside this walkthrough. A list of the sequences included in the analysis is [provided on github](./data/seqnames.txt). 

Cloning this github repo will allow file paths to correctly point to data and code in R. Set the working directory to the top level of the cloned repo.

### 1: Software and packages
You will need to install:
1. [BEAST2](https://www.beast2.org/) and associated package PhyDyn (for coupled MCMC runs we also use the BEASTLabs and CoupledMCMC packages)
2. [MAFFT](https://mafft.cbrc.jp/alignment/software/)
3. [IQ-TREE](http://www.iqtree.org/)
4. [TN93](https://github.com/veg/tn93)
5. [R](https://www.r-project.org/)

You will also need the CRAN R packages: `ape`, `treedater`, `phangorn`, `ggplot2`,`ggtreee`, `scales`, `lubridate`
and the utility functions package `sarscov2` from the [sarscov2Rutils github](https://github.com/emvolz-phylodynamics/sarscov2Rutils).

IQ-TREE and the BEAST2 bin/lib (windows/unix) folder will need to be on the PATH to work with the convenience wrapper functions in sarscov2Rutils. 

### 2: Alignment and filtering

Our analyses are carried out on sequence and metadata downloaded from the EpiCoV database on [gisaid.org](gisaid.org). These files are available from the downloads tab once you are registered for an account.
![](./images/gisaid_dash.PNG)

With the downloaded data we carry out the following steps, you may wish to change these:
1. Remove sequences from non-humans.
2. Remove sequences missing full date (no day / month / year given).
3. Remove gaps from sequence labels (needed for some software used).
4. Remove strings of 4+ N from sequences and replace with gaps (-).(This improves alignment but introduces other issues - may be less important if trying to align a smaller sample).
5. Drop 80% of UK sequences from March 25<sup>th</sup> onwards (helps to reduce the size of the alignment as there are so many UK sequences).
6. Align sequences with MAFFT. We currently add new sequences to the existing alignment with the --add tag rather than realigning the whole set each time.
7. Drop any position in the aligned sequences with >99% gaps (empty columns).
8. remove sequences with >10% gaps/Ns.
9. Clip sequences based on a reference (29400 bases starting from first ORF).
10. Build maximum likelihood phylogenies on subsets of the alignment (~1000 sequences) with an additional 500 sequences that are known to be well sequenced and aligned. (IQtree).
11. For each maximum likelihood phylogeny calculate cophenetic distance matrix (cophenetic.phylo in ape r package)  and drop new sequences from the alignment if their mean pairwise GD is >3 standard deviations from the phylogeny mean.
12. Calculate [TN93 distances](https://github.com/veg/tn93) and keep those <0.0001.

### 3: Editing downloaded metadata
1. Edit metadata file to include only sequences included in the final alignment.
2. Identify duplicate sequences. The majority of our analyses use only the deduplicated data with the oldest sequence from each set of duplicates retained. A binary column is added to the metadata file (inNoDups) to identify this subset of the sequences for use when building local alignments.
3. Date column standardised to sampleDate to account for differences in date format occasionally present.
4. 'region, country, division and location' columns renamed to 'Continent, Country, RegionOrState and CityOrCounty' to match previous database format, not required as these can be referenced by new column names in inclusion/exclusion rules code.


**NOTE:** As we are running analyses on any regions globally with a significant number of sequences (15-20+) we are downloading the whole GISAID database in this way. You may want to develop a simpler way to produce local alignments without first aligning the full database.

### 4: Producing BEAST2 xml

  * [Example script producing alingment and BEAST2 xml](./R/eg_xml_format.R)
  * [Template BEAST2 xml file used in our analyses](https://github.com/emvolz-phylodynamics/sarscov2Rutils/tree/sarscov2Rutils/inst/extdata)

From the complete gisaid alignment and metadata individual alignments are produced using the  `sarscov2Rutils` R package. These alignments are used to generate a set of starting trees for multiple beast runs which are inserted into xml templates available [here](https://github.com/emvolz-phylodynamics/sarscov2Rutils/tree/sarscov2Rutils/inst/extdata). These templates will be downloaded alongside the R package. 
A [sample R script](./R/eg_xml_format.R) has been produced to demonstrate these functions and their usage.
Following these steps should result in multiple BEAST2 xml files with the same alignment and different starting trees. They should also include the internal population size and start time of SEIJR dynamics set in the R script. The most recent reports use the `seijr0.1.3_skeleton_coupledMCMC.xml` template file.

The XML can also be formatted in the BEAST2 xml editor beauti with the Phydyn SEIR template.

### 5: Running BEAST2 

  * [Cluster submission guide for PBS](https://github.com/JorgensenD/BEAST_CLUSTER)

It is recommended that these beast analyses are run in parallel on an hpc platform. Instructions and example submission scripts for the Imperial PBS cluster are available [in the BEAST_CLUSTER github repo](https://github.com/JorgensenD/BEAST_CLUSTER). These scripts can be modified for other hpc platforms as needed.

The xml files can also be run with the desktop gui or command line versions of BEAST2 if preferred.

**NOTE:** By default these xml files will output identically named log, traj and trees files for each run. These will overwrite other runs if not renamed. This can be done by editing the xml for each run or by changing the names of the output files before copying back from the cluster nodes as is implemented with the `cat` commands in the [example PBS script](https://github.com/JorgensenD/BEAST_CLUSTER/blob/master/qsub_anaconda_array_resub_mc3.pbs).

### 6: Analyzing output

  * [Example script using analysis functions](./R/eg_analysis.R)
  * [Additional plotting functions](./R/combine_and_summarise0.R)

The current pipeline uses the reproducible reporting R package `orderly` to produce reports for publication.
Although versions of the figures presented in these reports can be reproduced using the `sarscov2Rutils` R package, the functions which reproduce the stylised versions of these plots found online are in [`combine_and_summarise0.R`](./R/combine_and_summarise0.R) in this repo. An [example R script](./R/eg_analysis.R) is provided to demonstrate combining log files from BEAST2 and visualizing the outputs using R. Indented sections of this code require you to have your own output from a beast2 run, other sections can be carried out using the provided example data files found in the [data folder](./data/) of this repo.

When running your own analyses we recommend checking log files manually with tracer for signs of stickiness, poor convergence etc. before running through the R code.


