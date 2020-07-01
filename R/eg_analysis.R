###########################################
##                                       ##
##      EG 2: Combine and summarize      ##
##                                       ##
###########################################




# required packages
require(sarscov2)
require(ape)
require(phangorn)
require(ggplot2)
require(ggtree)
require(scales)
require(lubridate)
require(dplyr)

# source combine and summarise file downloaded from github (set filepath to file location)
source('./R/combineAndSummarize0.R')

# Input some information for the plots and calculations
location <- 'Madinah'
date_lockdown <- NULL
first_internal <- as.Date('2020-03-29')
last_internal <- as.Date('2020-04-20')
last_sequence <- as.Date('2020-05-16')


# file path from the current working directory to the save location for the outputs
# all path_to_save or ofn commands can be set to NULL to not save to file.
path <- './eg_outputs/' 
dir.create(path)

# Input path pointing to existing data
in_path <- './data/'

################################################################################################
# The steps to combine log and traj files require log and traj files from your own beast runs  #
################################################################################################


    # we use multiple short chains and a large burnin to produce more rapid results. 
    burnin <- 50 #percentage

    # measure of lenience when including log and traj files, higher retains fewer. Setting to -1 will inlude all files
    pth <- 0.01

    # extract log and traj files produced from all of the beast runs by their file extentions 
    # If some traces show signs of stickiness or poor convergence they should be removed manually before this step
    logfns = list.files( pattern = '[0-9].xml.log$' )
    trajfns = list.files( pattern = '[0-9].xml.traj$' )


    # combine log and traj files and save each as .rds file
    combined <- combine_logs_and_traj(logfns, trajfns, burnProportion = burnin/100, ntraj=200, pth=pth,
                                                       ofn = paste0(path,'logs.rds') , ofntraj = paste0(path,'traj.rds'))


################################################################################################
# If you are using the uploaded log and traj files these functions can be used to create plots #
################################################################################################



# produce a table of R0, growth rate and doubling time and save to object R_num 

R_num <-  SEIJR_reproduction_number(paste0(in_path,'logs.rds'))

# Import and format reported cases for comparison in plots.
reported <- read.csv(paste0(in_path,'reportedcases.csv'), stringsAsFactors = F)
reported$Date <- lubridate::parse_date_time(reported$Date.of.confirmation , c("Ymd", "dmY", "mdY"))
reported$Date <- ymd(as.character(reported$Date))
# Names required by functions
colnames(reported) <- c('Date_char', 'Cumulative', 'Date')

# calculate daily change from the daily cumulative cases.
reported$Confirmed <- c(reported$Cumulative[1], diff(reported$Cumulative, lag=1))



# plot of estimated and reported cumulative cases over time
SEIJR_plot_size_returned <-  SEIJR_plot_size(trajdf = paste0(in_path,'traj.rds')
                    , case_data = reported
                    , date_limits = c( first_internal, last_sequence ) 
                    , path_to_save=paste0(path, 'size.png')
                    , last_tip = last_internal
                    , log_y_axis = T)
SEIJR_plot_size_returned$pl

# cumulative infections at last tip
infections_at_last_tip <- SEIJR_plot_size_returned$pldf[which(
  as.Date(SEIJR_plot_size_returned$pldf$Date) == last_internal &
    SEIJR_plot_size_returned$pldf$reported == FALSE )        ,]

infections_at_last_tip

# Plot of estimated and reported daily new infections 
daily_inf <- SEIJR_plot_daily_inf( paste0(in_path,'traj.rds')
                      , paste0(in_path,'logs.rds')
                      , case_data = reported
                      , date_limits = c( first_internal, last_sequence ) 
                      , path_to_save=paste0(path, 'daily_inf.png')
                      , last_tip = last_internal
                      , log_y_axis = T
)
daily_inf$pl

# plot fo estimated Rt over time
rt <- SEIJR_plot_Rt(paste0(in_path,'traj.rds')
                   , paste0(in_path,'logs.rds')
                   , gamma0 = 73
                   , gamma1 = 121.667
                   , date_limits = c(as.Date("2020-02-15"), NA)
                   , path_to_save = paste0(path, 'Rt.png')
                   , last_tip = last_internal
                   , lockdown_date = date_lockdown
)
rt$plot


# plot of proportion of total cases reported at each time point
rep <- SEIJR_plot_reporting(paste0(in_path,'traj.rds')
                , case_data = reported
                , date_limits = c(as.Date("2020-02-29"), NA)
                , path_to_save = paste0(path,'reporting.png')
                , errorbar = T)
rep$pl

################################################################################################
#       Producing an mcc tree requires your own tree files from your own beast analysis        #
################################################################################################

    ## Producing and using the mcc tree 

    treefiles <- list.files( pattern = '[0-9].xml.trees$')

    ## select only the tree files for which log and traj files were retained
    included <- paste0('\\.',c(1,10,12,13,14,15,17,18,19,2,20,3,6,7,8), '\\.')

    treefiles <- treefiles[grep(paste(included, collapse="|"), treefiles)]

    # run these through BEAST2 functions - these are R wrappers for the existing commnad line functions. 
    tree_combiner_helper(burnin=burnin, fns=treefiles, ofn = paste0(path,'combined.trees'), resample = 10000)
    # On windows tree annotator needs to be run after using log cobiner as there is no -trees option
    tree_annotator_windows(inputfile = paste0(path,'combined.trees'), outputfile = paste0(path,'mcc.nex'), lowMem = T)


################################################################################################
# The mcc tree included in the walkthrough can be used for these plots. Metadata for the tree  #
#     plot will need to be downloaded from gisaid (Continent == Region in newer datasets)      #
################################################################################################

# this nexus file can be used to plot the sampling distribution over time
sarscov2:::plot_sample_distribution(path_to_nex = paste0(in_path,'mcc.nex'), path_to_save = paste0(path, 'sample_distribution.png'))


# plotting the tree requires a table of continents for each of the included sequences.
continents <- md[c('seq_location', 'Continent')]
continents <- unique(continents)
names(continents) <- c('name', 'continent')
# Currently the function uses Americas rather than north and south america
continents$continent <- gsub('NorthAmerica|SouthAmerica', 'Americas', continents$continent)
# write to file
write.table(continents, paste0(path, "country_dict.txt"), col.names = T, row.names = F)

mcc_col_tree_plot(paste0(in_path,'mcc.nex'), mostRecentSampleDate = last_sequence, paste0(path, "country_dict.txt"), internal = location , style = 3, ofn=paste0(path,'mcc2.png'))


