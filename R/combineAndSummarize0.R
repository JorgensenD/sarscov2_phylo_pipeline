#' @export 
phylo_plot_aesthetics <- list(scale_x_datetime(labels = date_format("%B %d")),
                              theme_minimal(), 
                              xlab('')#,
                              #theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=10),
                                    #axis.title=element_text(size=10))
)


#' Compute R0, growth rate and doubling time for the SEIJR.0.0 model 
#'
#' Also prints to the screen a markdown table with the results. This can be copied into reports. 
#' The tau & p_h parameters _must_ be in the log files. If that's not the case, you can add fixed values like this: X$seir.tau <- 74; X$seir.p_h <- .2
#'
#" @param X a data frame with the posterior trace, can be produced by 'combine_logs' function
#' @param gamma0 Rate of leaving incubation period ( per capita per year )
#' @param gamma1 Rate of recovery (per capita per year)
#' @return Data frame with tabulated results and CI 
#' @export 
SEIJR_reproduction_number <- function( X, gamma0 = 73, gamma1 = 121.667,precision =3 ) {
  # tau = 74, p_h = 0.20 , 
  cat( 'Double check that you have provided the correct gamma0 and gamma1 parameters\n' )
  
  
  if(is.null(X$gamma1))
    X$gamma1 <- gamma1; 
  
  if(is.null(X$gamma0))
    X$gamma0 <- gamma0;
  
  if(is.null(X$seir.tau))
    X$seir.tau <- 74; 
  
  if(is.null(X$seir.p_h))
    X$seir.p_h <- .2
  
  Rs = ((1-X$seir.p_h)*X$seir.b/X$gamma1 + X$seir.tau*X$seir.p_h*X$seir.b/X$gamma1) 
  qR = signif( quantile( Rs, c(.5, .025, .975 )), precision )
  
  # growth rates 
  beta = (1-X$seir.p_h)*X$seir.b + X$seir.p_h*X$seir.tau*X$seir.b
  r = (-(X$gamma0 + X$gamma1) + sqrt( (X$gamma0-X$gamma1)^2 + 4*X$gamma0*beta )) / 2
  qr =  signif( quantile( r/365, c( .5, .025, .975 )), precision ) # growth rate per day 
  
  # double times 
  dr = r / 365 
  dbl = log(2) / dr 
  qdbl  = signif( quantile( dbl, c( .5, .025, .975 )), precision ) # days 
  
  #cat ( 'R0\n')
  #print( qR )
  
  O = data.frame( 
    `Reproduction number` = qR
    , `Growth rate (per day)` = qr
    , `Doubling time (days)` = qdbl
  )
  print( knitr::kable(O) )
  return(list(table = O,  parameters = X, R=Rs))
}



#' Plot the cumulative infections through time from a SEIJR trajectory sample
#' 
#" Also computes CIs of various dynamic variables 
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param case_data An optional dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Cumulative'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
SEIJR_plot_size <- function(trajdf
                            , case_data = NULL
                            , date_limits = c( as.Date( '2020-02-01'), NA ) 
                            , path_to_save='size.png'
                            , last_tip
                            , logscale = FALSE
                            , ...
) {
  library( ggplot2 ) 
  library( lubridate )
  
  if (!is.data.frame(trajdf)){
    # assume this is path to a rds 
    readRDS(trajdf) -> trajdf 
  }
  

  #~ 	r <- read.csv( '../seir21.0/weifangReported.csv', header=TRUE , stringsAsFactors=FALSE)
  #~ 	r$Date <- as.Date( r$Date )
  #~ 	r$reported = TRUE
  #~ 	r$`Cumulative confirmed` = r$Cumulative.confirmed.cases
  
  dfs <- split( trajdf, trajdf$Sample )
  taxis <- dfs[[1]]$t 
  
  
  
  
  
  if ( is.na( date_limits[2]) )
    date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
  
  qs <- c( .5, .025, .975 )
  
  # infectious
  Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
  Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
  I = Il + Ih 
  t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 
  
  # cases 
  cases <- do.call( cbind, lapply( dfs, '[[', 'infections' ))
  t(apply( cases, MAR=1, FUN=function(x) quantile(x,qs))) -> casesci 
  
  #exog 
  exog <- do.call( cbind, lapply( dfs, '[[', 'exog' ))
  t(apply( exog, MAR=1, FUN=function(x) quantile(x, qs )	)) -> exogci 
  
  #E
  E <- do.call( cbind, lapply( dfs, '[[', 'E' ))
  t(apply( E, MAR=1, FUN=function(x) quantile(x, qs )	)) -> Eci 
  
  
  #~ ------------
  pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
  pldf$`Cumulative infections` = casesci[,1]
  pldf$`2.5%` = casesci[,2]
  pldf$`97.5%` = casesci[,3] 
  
  if ( !is.null( case_data )){
    case_data$reported = TRUE
    pldf <- merge( pldf, case_data , all = TRUE ) 
  }
  
  
  pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= (date_limits[2] + 1)) , ]
  
  # pldf_test<<-pldf


  
  
  
  
  
  pl = ggplot( pldf ) +
    geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .2 , fill = "#A2A2A2") + 
    # geom_vline(xintercept =  as.POSIXct.Date(last_tip), linetype = "dashed")+
    geom_path( aes(x = Date, y = `Cumulative infections` , group = !reported), lwd=0.75, col = "gray18", linetype = "solid") +
    geom_path( aes(x = Date, y = `2.5%` , group = !reported), lwd=0.6, col = "#A2A2A2") +
    geom_path( aes(x = Date, y = `97.5%` , group = !reported), lwd=0.6, col = "#A2A2A2") +
    phylo_plot_aesthetics +
    ylab ('Cumulative estimated infections'  ) +
    # geom_label(label = round(pldf[which( as.Date(pldf$Date) == last_tip & pldf$reported == FALSE )  , "Cumulative infections"], 0),
    #           aes(x = pldf[which( as.Date(pldf$Date) == last_tip & pldf$reported == FALSE )  , "Date"] + lubridate::days(4)
    #                        , y = pldf[which( as.Date(pldf$Date) == last_tip & pldf$reported == FALSE )  , "Cumulative infections"]*1.2),
    #          fill = "#0e0a00", alpha = 0 , col = "#0e0a00", size = 5)
    
    geom_segment(aes(x = as.POSIXct.Date(last_tip), 
                     y = pldf[which( as.Date(pldf$Date) == last_tip & pldf$reported == FALSE ),"Cumulative infections"]
    ) , xend = 0, yend = 0, col = "black", linetype = "dashed") +
    
    geom_segment(aes(x = as.POSIXct.Date(last_tip), 
                     y = pldf[which( as.Date(pldf$Date) == last_tip & pldf$reported == FALSE ),"Cumulative infections"]
    ) , xend = as.POSIXct.Date(last_tip), yend = 0, col = "black", linetype = "dashed")+
    geom_point(size = 3,shape = 24,
               aes(x = as.POSIXct.Date(last_tip),
                   y = pldf[which( as.Date(pldf$Date) == last_tip & pldf$reported == FALSE )  , "Cumulative infections"]),
               col="black", fill = "gray18")
  
  if ( !is.null(case_data) ) {
    if( !(class(case_data$Date)=='Date') ){
      stop('case_data Date variable must be class *Date*, not character, integer, or POSIXct. ')
    }
    pl <- pl +
      # geom_segment(aes(x = Date, y = Cumulative, xend=Date) , size = 0.5, alpha = 0.5, yend = 0, col = "#E69F00", linetype = "dashed") +
      geom_point( aes( x = Date, y = Cumulative ) , size =2, shape = 21, fill =  "#E69F00", col = "black") +
      geom_path(aes( x = Date, y = Cumulative, group = reported), col ="#E69F00")
  }
  
  if (logscale == TRUE) {
    pl <- pl + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x)), 
                             limits = c(1, NA))+
      ylab ('Cumulative estimated infections (log scale)' )
  }

  
  
  
  
  if (!is.null(path_to_save))
    ggsave(pl, file = path_to_save)
  
  list(
    pl = pl
    , taxis = taxis 
    , Il = Il
    , Ih = Ih
    , E = E 
    , I = I 
    , cases = cases 
    , exog = exog 
    , pldf =pldf 
    , case_data = case_data 
  )
}




#' Plot the daily new infections through time from a SEIJR trajectory sample
#' 
#" Also computes CIs of various dynamic variables 
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param logdf Either a dataframe or a path to rds containing a data frame with posterior logs
#' @param case_data An optional dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Confirmed'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
SEIJR_plot_daily_inf <- function(trajdf
                                 , logdf
                                 , case_data = NULL
                                 , date_limits = c( as.Date( '2020-02-01'), NA ) 
                                 , path_to_save='daily.png'
                                 , last_tip 
                                 , logscale = FALSE
                                 , ...
) {
  library( ggplot2 ) 
  library( lubridate )
  
  # browser()
  

  if (!is.data.frame(trajdf)){
    # assume this is path to a rds
    readRDS(trajdf) -> trajdf
  }
  


  X <- logdf

  
  if(is.null(X$seir.tau))
    X$seir.tau <- 74; 
  
  if(is.null(X$seir.p_h))
    X$seir.p_h <- .2	
  
  
  dfs <- split( trajdf, trajdf$Sample )
  taxis <- dfs[[1]]$t 
  
  s <- sapply( dfs, function(df) df$Sample[1] )
  X1 <- X[match( s, X$Sample ), ]
  
  if ( is.na( date_limits[2]) )
    date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
  
  qs <- c( .5, .025, .975 )
  
  # infectious
  Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
  Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
  I = Il + Ih 
  t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 
  
  # daily new inf 
  Y = lapply( 1:length(dfs) , function(k) {
    x = X1[k, ]
    Il = dfs[[k]]$Il
    Ih = dfs[[k]]$Ih
    S = dfs[[k]]$S
    E = dfs[[k]]$E
    R = dfs[[k]]$R
    tau = X1$seir.tau[k]
    p_h = X1$seir.p_h[k] 
    b = X1$seir.b[k]
    y = (1/365) * ( b * Il + b * tau * Ih   ) *  S/ ( S + E + Il + Ih + R )
  })
  Y = do.call( cbind, Y )
  Yci = t(apply( Y, MAR=1, FUN= function(x) quantile(na.omit(x), qs ))) 
  
  #~ ------------
  pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
  pldf$`New infections` = Yci[,1]
  pldf$`2.5%` = Yci[,2]
  pldf$`97.5%` = Yci[,3] 
  
  if ( !is.null( case_data )){
    case_data$reported = TRUE
    pldf <- merge( pldf, case_data , all = TRUE ) 
  }
  
  pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
 
  # pldf_test<<-pldf
  
  
  pl = ggplot( pldf ) + 
    geom_vline(xintercept =  as.POSIXct.Date(last_tip), linetype = "dashed")+
    geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .2 , fill = "#A2A2A2") + 
    geom_path( aes(x = Date, y = `New infections` , group = !reported), lwd=0.75, col = "gray18", linetype = "solid") +
    geom_path( aes(x = Date, y = `2.5%` , group = !reported), lwd=0.6, col = "#A2A2A2") +
    geom_path( aes(x = Date, y = `97.5%` , group = !reported), lwd=0.6, col = "#A2A2A2") +
    phylo_plot_aesthetics +
    ylab ('Estimated daily new infections')  
  
  if (logscale == TRUE) {
    pl <- pl + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                             labels = scales::trans_format("log10", scales::math_format(10^.x) ))+ 
      ylab ('Estimated daily new infections (log scale) ')
  }
  
  
  
  if ( !is.null(case_data) ) {
    if( !(class(case_data$Date)=='Date') ){
      stop('case_data Date variable must be class *Date*, not character, integer, or POSIXct. ')
    }
    pl <- pl +
      # geom_segment(aes(x = Date, y = daily_cases, xend=Date) , size = 0.5, alpha = 0.5, yend = 0, col = "#E69F00", linetype = "dashed") +
      geom_point( aes( x = Date, y = daily_cases ) , size =2, shape = 21, fill =  "#E69F00", col = "black") #+
    # geom_path(aes( x = Date, y = daily_cases, group = reported), col ="#E69F00")
  }
  

  
  
  
  
  if (!is.null(path_to_save))
    ggsave(pl, file = path_to_save)
  
  list(
    pl = pl
    , taxis = taxis 
    , Il = Il
    , Ih = Ih
    , Y
    , pldf =pldf 
    , case_data = case_data 
  )
}

#' Plot reproduction number through time from a SEIJR trajectory sample
#'
#" @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param logdf Either a dataframe or a path to rds containing a data frame with posterior logs
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param gamma0 rate of becoming infectious during incubation period
#' @param gamma1 rate of recovery once infectious
#' @param path_to_save Will save a png here 
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 
SEIJR_plot_Rt <- function(trajdf
                          , logdf
                          , gamma0 = 73, gamma1 = 121.667
                          , date_limits = c( as.Date( '2020-02-01'), NA ) 
                          , path_to_save='Rt.png'
                          , last_tip
                          , school_closure_date = NULL
                          , lockdown_date = NULL
                          , ...
) {
  library( ggplot2 ) 
  library( lubridate )
  
  
  if (!is.data.frame(trajdf)){
    # assume this is path to a rds 
    readRDS(trajdf) -> trajdf 
  }
  
  X <- logdf
  
 
  
  if(is.null(X$seir.tau))
    X$seir.tau <- 74; 
  
  if(is.null(X$seir.p_h))
    X$seir.p_h <- .2	
  
  dfs <- split( trajdf, trajdf$Sample )
  taxis <- dfs[[1]]$t 
  
  s <- sapply( dfs, function(df) df$Sample[1] )
  X1 <- X[match( s, X$Sample ), ]
  
  if ( is.na( date_limits[2]) )
    date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
  
  qs <- c( .5, .025, .975 )
  
  Rs = ((1-X1$seir.p_h)*X1$seir.b/gamma1 + X1$seir.tau*X1$seir.p_h*X1$seir.b/gamma1)
  
  # infectious
  Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
  Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
  S = do.call( cbind, lapply( dfs, '[[', 'S' ))
  E = do.call( cbind, lapply( dfs, '[[', 'E' ))
  R = do.call( cbind, lapply( dfs, '[[', 'R' ))
  I = Il + Ih 
  pS <- S / ( S + E + I + R )
  Rts = t(Rs * t(pS))
  t(apply( Rts, MAR=1, FUN= function(x) quantile(na.omit(x), qs ))) -> Yci  
  
  #~ browser()
  #~ ------------
  pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
  pldf$`R(t)` = Yci[,1]
  pldf$`2.5%` = Yci[,2]
  pldf$`97.5%` = Yci[,3] 
  
  
  pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
  
  # pldf_test<<-pldf
  
  
  pl = ggplot( pldf ) + 
    # geom_path( aes(x = Date, y = `R(t)` , group = !reported), lwd=1.25) + 
    # geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .25 ) +
    
    geom_vline(xintercept =  as.POSIXct.Date(last_tip), linetype = "dashed")+
    geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`, group = !reported) , alpha = .2 , fill = "#A2A2A2") + 
    geom_path( aes(x = Date, y = `R(t)` , group = !reported), lwd=0.75, col = "gray18", linetype = "solid") +
    geom_path( aes(x = Date, y = `2.5%` , group = !reported), lwd=0.6, col = "#A2A2A2") +
    geom_path( aes(x = Date, y = `97.5%` , group = !reported), lwd=0.6, col = "#A2A2A2") +
    phylo_plot_aesthetics +
    
    geom_hline(yintercept = 1, linetype = "dotted") +
    ylab ('Effective reproduction number through time R(t)' )
  
  if (!is.null(school_closure_date))
    pl <- pl + geom_vline(xintercept =  school_closure_date, linetype = "dashed", col = "orange")
  
  if (!is.null(lockdown_date))
    pl <- pl + geom_vline(xintercept =  lockdown_date, linetype = "dashed", col = "darkred")
  
  if (!is.null(path_to_save))
    ggsave(pl, file = path_to_save)
  
  list(
    plot = pl
    , taxis = taxis 
    , Il = Il
    , Ih = Ih
    , Rts
    , pldf =pldf 
  )
}

#' Plot reporting rate through time from SEIJR trajectory sample and reported cases
#' 
#'
#' @param trajdf Either a dataframe or a path to rds containing a data frame with a posterior sample of trajectories (see combine_traj)
#' @param case_data dataframe containing reported/confirmed cases to be plotted alongside estimates. *Must* contain columns 'Date' and 'Confirmed'. Ensure Date is not a factor or character (see as.Date )
#' @param date_limits  a 2-vector containing bounds for the plotting window. If the upper bound is missing, will use the maximum time in the trajectories
#' @param path_to_save Will save a png here 
#' @return a list with derived outputs from the trajectories. The first element is a ggplot object if you want to further customize the figure 
#' @export 

SEIJR_plot_reporting <- function(trajdf
                                 , case_data
                                 , date_limits = c( as.Date( '2020-02-01'), NA ) 
                                 , path_to_save='reporting.png'
                                 , ...
) {
  library( ggplot2 ) 
  library( lubridate )
  library( dplyr )
  
  if (!is.data.frame(trajdf)){
    # assume this is path to a rds 
    readRDS(trajdf) -> trajdf 
  }
  
  dfs <- split( trajdf, trajdf$Sample )
  taxis <- dfs[[1]]$t 
  
  if ( is.na( date_limits[2]) )
    date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
  
  qs <- c( .5, .025, .975 )
  
  # infectious
  Il = do.call( cbind, lapply( dfs, '[[', 'Il' ))
  Ih = do.call( cbind, lapply( dfs, '[[', 'Ih' ))
  I = Il + Ih 
  t(apply( I, MAR=1, FUN= function(x) quantile(x, qs ))) -> Ici 
  
  # cases 
  cases <- do.call( cbind, lapply( dfs, '[[', 'infections' ))
  t(apply( cases, MAR=1, FUN=function(x) quantile(x,qs))) -> casesci 
  
  #exog 
  exog <- do.call( cbind, lapply( dfs, '[[', 'exog' ))
  t(apply( exog, MAR=1, FUN=function(x) quantile(x, qs )	)) -> exogci 
  
  #E
  E <- do.call( cbind, lapply( dfs, '[[', 'E' ))
  t(apply( E, MAR=1, FUN=function(x) quantile(x, qs )	)) -> Eci 
  
  
  #~ ------------
  pldf <- data.frame( Date = ( date_decimal( taxis ) ) , reported=FALSE )
  pldf$`Cumulative infections` = casesci[,1]
  pldf$`2.5%` = casesci[,2]
  pldf$`97.5%` = casesci[,3] 
  
  if( !(class(case_data$Date)=='Date') ){
    stop('case_data Date variable must be class *Date*, not character, integer, or POSIXct. ') 
  }
  
  case_data$reported = TRUE
  
  pldf <- merge( pldf, case_data , all = TRUE ) 
  pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
  pldf$reporting <- lead(pldf$Cumulative, n=1L)/pldf$`Cumulative infections`
  pldf$`rep97.5` <- pmin(lead(pldf$Cumulative, n=1L)/pldf$`2.5%`, 1)
  pldf$`rep2.5` <- lead(pldf$Cumulative, n=1L)/pldf$`97.5%`
  
  pldf_test<<-pldf
  
  
  pl = ggplot( pldf, aes(x = Date, ymin = `rep2.5`*100, ymax = `rep97.5`*100, group = !reported) ) + 
    geom_errorbar( width = 2, 
                   col = "black", alpha = 0.6,size=0.8)+
    
    geom_point(aes(x = Date, y = `reporting`*100, 
                   group = !reported),   size =2, shape = 21, fill =  "#E69F00", col = "black") + 
    
    ylab("Percentage of cases reported (%)")+
    phylo_plot_aesthetics
  
  if (!is.null(path_to_save))
    ggsave(pl, file = path_to_save)
  
  list(
    pl = pl
    , taxis = taxis 
    , Il = Il
    , Ih = Ih
    , E = E 
    , I = I 
    , cases = cases 
    , exog = exog 
    , pldf =pldf 
    , case_data = case_data 
  )
}

#' Sampling distribution of sequences
#' 
#' @param path_to_nex Path to mcc nexus tree
#' @param path_to_save PNG saved here 
plot_sample_distribution = function(path_to_nex, path_to_save = 'sample_distribution.png' ) {
  
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(ape)
  
  #parse dates/location 
  nex <- read.nexus(path_to_nex)
  algn3 = data.frame(seq_id = nex$tip.label) %>% 
    separate(seq_id, c("hcov", "country", "region_code", "year", "EPI", "Date", "Date_2", "Region"), "[\\/\\|//]") %>% 
    mutate(Region = ifelse(Region == "_Il", "Local", "Global")) %>% 
    mutate(Date = as.Date(Date, "%Y-%m-%d"))
  
  pl = ggplot(algn3, aes(x = Date, color= Region, fill = Region)) +
    geom_histogram(aes(y=..density..), position="identity", 
                   alpha=0.5, bins = as.numeric(max(algn3$Date)-min(algn3$Date)) + 1) +
    geom_density(alpha=0.4) + 
    theme_minimal() +
    scale_color_manual(values=c( "#999999", "#E69F00"))+
    scale_fill_manual(values=c( "#999999", "#E69F00"))+
    labs(x="", y = "Sampling density")
  
  
  if (!is.null(path_to_save))
    ggsave(pl, file = path_to_save)
  
  return(pl)
}