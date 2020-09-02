###########################################
##                                       ##
##            EG 1: BEAST2 XML           ##
##                                       ##
###########################################

# Set WD to local clone of github repo

# Install R packages
install.packages(c('devtools', 'ape', 'phangorn', 'ggplot2', 'ggtree', 'lubridate', 'limSolve'))
devtools::install_github('emvolz-phylodynamics/sarscov2Rutils', ref = 'sarscov2Rutils')

# load packages
require(sarscov2)
require(ape)
require(phangorn)
require(ggplot2)
require(ggtree)
require(lubridate)
require(limSolve)

# define some parameters 
## path to the GISAID alignment: 
fn = 'gisaid.fas' 
## A simple name for the region:
region='Madinah'
## A path to the file containing small genetic distance pairs
distfn = 'tn93dist.txt'
## When do internal SEIR dynamics initiate: 
startTime = 2020.15
## What is the populatio size of the region of interest?
popSize = 1100000
## How many sequences to include within the region?
n_region = 50
## How many to include from the intertional reservoir?
n_reservoir = 50 # (actual number to be included will be greater since it also includes close distance matches)
## How many starting trees? A BEAST run will be carried out for each 
n_startingtrees = 20 

# Load the editied gisaid metadata 
md <- read.csv( 'md.csv', stringsAs=FALSE ) 

# Set dedup = TRUE in region sampler and exog sampler to remove duplicate sequences

# Let's sample from Madinah
regiontips = region_sampler1( md, n = n_region  , inclusion_rules = list( c('RegionOrState', '^Madinah$') ), dedup = FALSE)

# Alternatively you can include multiple inclusion criteria:
# regiontips = region_sampler1( md, n = n_region
#                               , inclusion_rules = 
#                                 list( 
#                                   c('CityOrCounty', '^KingCounty$')
#                                   ,  c('CityOrCounty', '^Washington$')
#                                 ) 
# )

# Sample a set of closely related external sequences, 
exogtips = exog_sampler2( md, n=n_reservoir, smallGDpairs='tn93dist.txt', region_sample=regiontips, exclusion_rules = list( c('RegionOrState', '^Madinah$') ), dedup= FALSE )

# You could select based on a different geographic criterion, current fields in metadata are  "Continent" "Country" "RegionOrState" "CityOrCounty"

#~ Here is a subset of the metadata that corresponds to the sample. Check this carefully- is the sample including everything you want and not including things you dont want?
mdregion = md[ match( regiontips, md$seqName ) , ]


# Now make the aligment for BEAST. This adds sample time and deme (internal or external) to tip labels:
d3 = prep_tip_labels_seijr( 'gisaid.fas', outfn = 'algn3.fasta', regiontips = regiontips, exogtips = exogtips, metadata = md  )

################################################################################################
#   We strongly recommend checking the tree and alignment producd by this step and removing    #
#                          outlying sequences caused by misalignment                           #
################################################################################################

# Make the starting trees
# Note requires IQTREE, ncpu > 1 will not work on windows
#tds = make_starting_trees('algn3.fasta', treeoutfn = "startTrees.nwk", plotout = 'MLtree.png', ntres = n_startingtrees, ncpu = 6) # have a fix for plotout not yet pushed to repo
tds = make_starting_trees('algn3.fasta', treeoutfn = "startTrees.nwk", ntres = n_startingtrees, ncpu = 6)

file.copy( system.file( package='sarscov2', 'extdata/seijr0.1.3_skeleton.xml') , paste0(region, '.xml') )

readline(prompt="At this point you should check the skeleton xml and make any necessary changes to prior distributions. You can also do this in beauti. Press any key to continue when this is done. ")

# Make the xmls 
format_xml0( xmlfn=paste0(region, '.xml'), fastafn = 'algn3.fasta', treefn='startTrees.nwk', start = startTime, susc_size = popSize) 
