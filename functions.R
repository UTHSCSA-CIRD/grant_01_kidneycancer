#' Return a factor containing the top-N levels, and the rest binned
#' into the specified existing or new level. That level is listed last
#' in the resulting factor.
#' @param xx A \code{vector} (required).
#' @param topn \code{numeric} for the top how many levels to keep (optional, default =4)
#' @param binto \code{character} which new or existing level to dump the other values in
cl_bintail <- function(xx,topn=4,binto='other'){
  if(!is.factor(xx)) xx <- factor(xx);
  counts <- sort(table(xx),decreasing = T);
  if(is.numeric(binto)) binto <- names(counts)[binto];
  keep <- setdiff(names(counts)[1:min(length(counts),topn)],binto);
  droplevels(
    factor(
      ifelse(
        xx %in% keep, as.character(xx), binto
        ),levels=c(keep,binto)));
}

#' Take a character vector and perform multiple search-replace 
#' operations on it.
#' @param xx A \code{vector} of type \code{character} (required)
#' @param searchrep A \code{matrix} with two columns of type \code{character} (required). The left column is the pattern and the right, the replacement.
#' @param method One of 'partial','full', or 'exact'. Controls whether to replace only the matching regexp, replace the entire value that contains a matching regexp, or replace the entire value if it's an exact match.
submulti <- function(xx,searchrep,method=c('partial','full','exact')){
  # if no method is specified by the user, this makes it take the first value
  # if a method is only partially written out, this completes it, and if the
  # method doesn't match any of the valid possibilities this gives an informativ
  # error message
  method<-match.arg(method);
  # rr is a sequence of integers spanning the rows of the searchrep matrix
  rr <- 1:nrow(searchrep);
  # oo will hold the result that this function will return
  oo <- xx;
  switch(method
         ,partial = {for(ii in rr)
           oo <- gsub(searchrep[ii,1],searchrep[ii,2],oo)}
         ,full =    {for(ii in rr)
           oo[grepl(searchrep[ii,1],oo)]<-searchrep[ii,2]}
         ,exact = {for(ii in rr)
           oo[grepl(searchrep[ii,1],oo,fixed=T)]<-searchrep[ii,2]}
           #oo <- gsub(searchrep[ii,1],searchrep[ii,2],oo,fixed = T)}
         );
  oo;
}

#' Take a data.frame and create a PL/SQL table definition for it, ready to paste into DBVis and such
#' @param xx A \code{data.frame} (required)
sql_create_table <- function(xx,tname,sqltemplate='CREATE TABLE %s (%s);'){
  mysearchrep <- cbind(c('integer','character','numeric','Date','logical','factor')
    ,c('NUMBER(32)','VARCHAR2(nn)','NUMBER(32)','DATE','NUMBER(1)','VARCHAR2(nn)'));
  if(missing(tname)) tname <- as.character(substitute(xx));
  thenames <- sapply(xx, class);
  thenames <- submulti(thenames, mysearchrep, 'exact');
  strlengthres<-apply(nsqip,2,function(xx){max(str_length(xx))});
  thenames <- apply(cbind(round(strlengthres*1.5), thenames),1,function(xx){gsub('nn',xx[1],xx[2])});
  thenames <- paste(' ',names(thenames), as.character(thenames), sep=' ', collapse= ',\n');
  thenames <- gsub('.', '_', thenames, fixed=TRUE);
  thenames <- gsub('_{1,}','_',thenames);
  thenames <- gsub('_ ',' ',thenames);
  thenames <- sprintf(sqltemplate, tname, thenames);
  cat(thenames);
  invisible(thenames);
  # TODO: collect the column names and classes as we did in console
  # TODO: use the above and submulti() to convert R classes to SQL data types
  # TODO: paste() and/or sprintf() to create the string that will be the output
}

#' Take a data.frame or character vector and a vector of grep targets and return
#' the values that match (for data.frame, column names that match). If no 
#' patterns given just returns the names
#' @param xx A \code{data.frame} or character vector (required)
#' @param patterns A character vector of regexp targets to be OR-ed
grepor <- function(xx,patterns='.') {
  if(is.list(xx)) xx <-names(xx);
  grep(paste0(patterns,collapse='|'),xx,val=T);
}

#' Delete all the junk in your environment, for testing
clearenv <- function(env=.GlobalEnv) rm(list=setdiff(ls(all=T,envir=env),'clearenv'),envir=env);

#' From ... http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  # Make the panel
  # ncol: Number of columns of plots
  # nrow: Number of rows needed, calculated from # of cols
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))}
  if (numPlots==1) print(plots[[1]]) else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} 

#' ### Specific to RAI-A
#' 

#' Calculates the RAI score
#' 
raiscore <- function(xx){
  with(xx, 5*(gender=='Male') +
    # note that there are two different multipliers for body weight loss
    (5+4)*(x_loss_bw_6_months_prior_surg=='Yes') +
    6*(currently_dialysis=='Yes'|isTRUE(serum_creatinine>3)) +
    4*(chr_30_dy_prior_surg=='Yes') +
    8*(dyspnea == 'At Rest') +
    # Points to addres: 
    # 1. We do not have impaired sensorium or other cog decline 
    # available. Do we assume NO cognitive decline?  Or do we 
    # take a weighted average of the two health status scores 
    # based on the population average expected for cognitive 
    # decline?
    # 2. How do we treat unknown functional status? Currently treated as independent.
    as.numeric(
      as.character(
        factor(functnal_heath_status
               ,levels=c('Independent','Partially Dependent','Totally Dependent','Unknown')
               # with cognitive decline the labels would have been...
               #,labels=c(-2,10,21,'0.0')
               ,labels=c(0,8,16,'0.0')
               ))) +
    # 3. We do not have impaired sensorium or other cog decline here
    # transfer status
    8*(origin_status!='Not transferred (admitted from home)') +
    ifelse(disseminated_cancer=='Yes'
           # scoring for age when there is cancer
           ,as.numeric(cut(dat1$age_at_time_surg,breaks = c(0,69,74,79,84,89,99,Inf),labels = c(20,19,18,17,15,14,13)))
           # scoring for age when there is no cancer
           ,as.numeric(cut(dat1$age_at_time_surg,breaks = c(0,69,74,79,84,89,94,99,Inf),labels = 2:9))
    ))
  }
  
raiscore.bak <- function(xx){
  with(xx,
       rep(0,nrow(xx)) +
         ifelse(gender=='Male',5,0) +
         ifelse(vent_dependent=='Yes',1,0) +
         ifelse(currently_dialysis=='Yes'|acute_renal_failure=='Yes',6,0) +
         ifelse(chr_30_dy_prior_surg=='Yes',4,0) +
         ifelse(x_loss_bw_6_months_prior_surg=='Yes',5,0) +
         ifelse(dyspnea == 'At Rest',8,0) +
         ifelse(preop_tnsf_rbc_72h_prior_surg=='Yes',8,0) +
         ifelse(functnal_heath_status=='Independent',-2,0) +
         ifelse(functnal_heath_status=='Partially Dependent',10,0) +
         ifelse(functnal_heath_status=='Totally Dependent',21,0) +
         ifelse(disseminated_cancer=='Yes'
                ,as.numeric(cut(dat1$age_at_time_surg,breaks = c(0,69,74,79,84,89,99,Inf),labels = c(20,19,18,17,15,14,13)))
                ,as.numeric(cut(dat1$age_at_time_surg,breaks = c(0,69,74,79,84,89,94,99,Inf),labels = 2:9))
                )
       );
}

#' Returns a list of column names from the data dictionary for which the column
#' named in the first argument is true. The first arg can be either a string or 
#' a name. The second must be a data.frame
v <- function(var,dictionary=dct0) {cc<-substitute(var);na.omit(dictionary[dictionary[[as.character(cc)]],'dataset_column_names'])[[1]]}
