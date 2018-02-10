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

#' Take an object name \code{obj}, check to see if it  exists in environment \code{env}
#' and if it does not, run \code{expression} \code{EXPR} and assign its result to \code{obj}
#'
#' @param obj   A \code{character} string (required) naming the variable in env
#'   to which the results of running \code{EXPR} will be assigned.
#' @param EXPR  An \code{expression} to execute and assign to the object named
#'   by \code{obj} if it doesn't already exist.
#' @param env   An \code{environment} to check for the existence of \code{obj}
#'   and where to create it if it doesn't exist.
#'
#' @return Does not return anything.
#'
#' @examples `checkrun('dat3',{group_by(dat1,idn_mrn) %>% summarise_all(first)});`
checkrun <- function(obj,EXPR,env=as.environment(-1)){
  env<-env;
  if(length(ls(env,all=T,pattern=obj))==0){
    browser();
  }
}


#' Delete all the junk in your environment, for testing
clearenv <- function(env=.GlobalEnv) rm(list=setdiff(ls(all=T,envir=env),'clearenv'),envir=env);

#' Plots in the style we've been doing (continuous y, discrete x and optionally z)
#' 
#' Instead of creating new tables for 'All', just set xx=T
#' 
#' To suppress plotting of a legend (but still use fill) set fill.name=NA
#' Likewise, to suppress printing y-axis labels make yy.name=NA
#' To merely omit printing a name, set the name in question to ''
autoboxplot <- function(pdata, xx, yy, zz, subset=T
                        , type=c('box','violin')
                        , title=sprintf('%s vs %s\n by %s',xx,yy,zz)
                        , xx.name=if(xx==TRUE) 'All' else xx, xx.breaks=if(xx==TRUE) xx else unique(pdata[[xx]])
                        , xx.labels=if(xx==TRUE) '' else xx.breaks
                        , yy.name=yy, yy.labels
                        , fill.name, fill.breaks, fill.labels
                        , counts=T
                        ,...){
  subset <- substitute(subset);
  plot_type <- switch(match.arg(type)
                      ,box=geom_boxplot(coef=100)
                      ,violin=geom_violin());
  pdata <- subset(pdata,subset=eval(subset));
  if(!missing(zz)){
    if(missing(fill.name)) fill.name <- zz;
    if(missing(fill.breaks)) fill.breaks <- unique(pdata[[zz]]);
    if(missing(fill.labels)) fill.labels <- fill.breaks;
    fill <- if(is.atomic(fill.name)&&is.na(fill.name)) scale_fill_discrete(guide=F) else {
      scale_fill_discrete(name=fill.name,breaks=fill.breaks,labels=fill.labels);
    }
  }
  if(missing(yy.labels)){
    yy.labels <- if(is.na(yy.name)) scale_y_continuous(name='',labels=NULL) else {
      scale_y_continuous(name=yy.name,labels = comma);
      };
    };
  out <- ggplot(data = pdata, aes_string(x = xx, y = yy)) + plot_type + yy.labels +
    scale_x_discrete(name=xx.name,breaks=xx.breaks,labels=xx.labels) +
    labs(title=title)
  if(!missing(zz)) out <- out + aes_string(fill=zz) + fill;
  if(counts){
    ccrds<-ggplot_build(out)$layout$panel_ranges[[1]];
    ann.label <- if(xx==TRUE && missing(zz)) nrow(pdata) else if(xx==TRUE){
      paste0(table(pdata[,zz]),collapse=' \t ') } else if(missing(zz)){
        paste0(table(pdata[,xx]),collapse=' \t ') } else {
          apply(table(pdata[,c(xx,zz)]),1,function(ii) paste0(ii[!ii%in%list(0,'')],collapse=' \t '));
        }
    #ann.label <- gsub('\\b0\\b','',ann.label);
    out <- out + annotate('text',x=ccrds$x.major_source,y=ccrds$y.range[1]
                          ,label=ann.label[ccrds$x.major_source]);
  }
  attr(out,'call') <- match.call();
  out;
}

getCall.data.frame <- getCall.gg <- function(xx) attr(xx,'call');


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

mktabsie <- function(data,subsets=list(Full=T),pw
                     ,vars
                     ,filepath='.'
                     ,filename='survSave.rdata'
                     ,serverTitle='TABSIE'
                     ,serverStatement=bquote(h4("Welcome to TABSIE"))){
  serverData <- sapply(subsets,function(ii) subset(data[,vars],eval(ii)),simplify=F);
  serverDataDic <- names(serverData);
  serverHash <- digest::digest(pw,algo='sha512',ascii=TRUE);
  save(serverStatement,serverData,serverDataDic,serverTitle,serverHash
       ,file=paste0(filepath,'/',filename));
}

#' Returns a vector of column names that contain data elements of a particular type
#' as specified by the user: "integer","POSIXct" "POSIXt", "numeric", "character", 
#' "factor" and "logical". 
vartype <- function(dat, ctype) {
  xx <- unlist(sapply(dat, class));
  idx <- which(xx %in% ctype);
  res <- names(xx)[idx];
  return(res)
}

#' ## Functions for RMarkdown and ggplot2
#' 
#' Return a commit hash (for inclusion in reports for example) after first making
#' sure all changes are committed and pushed
#' TODO: instead of auto-committing, error if uncommited changes, needs to be 
#' a deliberate process, otherwise we have tons of meaningless auto-commit
#' messages that will make future maintenance harder
gitstamp <- function(production=T) {
  if(production){
    if(length(gitdiff<-system("git update-index --refresh && git diff-index HEAD --",intern = T))!=0) stop(sprintf(
      "\ngit message: %s\n\nYou have uncommitted changes. Please do 'git commit' and then try again."
      ,gitdiff));
    system("git push && git log --pretty=format:'%h' -n 1",intern=T);
  } else system("git log --pretty=format:'%h' -n 1",intern=T);
}

#' This function can be called from `stat_summary()` as the
#' `fun.data=` argument. It will cause group counts to be 
#' over-printed on a `geom_boxplot()` (probably other similar
#' plots too) if `stat_summary()` is added to it.
n_fun <- function(xx) data.frame(y=mean(quantile(xx,c(.5,.75))),label=as.character(length(xx)));

#' take a list of subset criteria and return a list of data.frames
ssply<-function(dat,...) sapply(sys.call()[-(1:2)],function(ii) subset(dat,eval(ii)),simplify=F);

#' save a named list of tables, including names
savetablelist <- function(lst,fileprefix,filesuffix=paste0(format(Sys.Date(),'%m-%d-%Y-%H_%M_%S'),'.tsv')
                            ,filepath='./',outfile=paste0(filepath,fileprefix,filesuffix)
                          ,sep='\t',row.names=F
                          ,singletabname=as.character(substitute(lst))
                          ,replaceifexists=T,...) {
  if(is.data.frame(lst)) lst<-setNames(list(lst),singletabname);
  if(replaceifexists&&file.exists(outfile)) file.remove(outfile);
  for(ii in names(lst)){
  write(ii,file=outfile,append=T);
  write.table(lst[[ii]],outfile,sep=sep,row.names=row.names,append=T);
  }
}

#' Usage: `xx<-mapnames(xx,lookup)` where lookup is a named character vector
#' the names of the elements in the character vector are what you are renaming
#' things TO and the values are what needs to be matched, i.e. what renaming things
#' FROM. If you set namesonly=T then it just returns the names, not the original
#' object.
mapnames<-function(xx,lookup,namesonly=F,...){
  xnames <- names(xx);
  idxmatch <- na.omit(match(xnames,lookup));
  newnames <- names(lookup)[idxmatch];
  if(namesonly) return(newnames);
  names(xx)[xnames %in% lookup] <- newnames;
  xx;
}

#' Forget modifying the columns to TRUE/FALSE... some of them have 'Unknown',
#' NA, who knows what else. We can just do it dynamically when we need to. 
#' This might even let us get rid of a couple of analytic columns.

#' Example of using R methods dispatch
#' 
#' The actual usage is: `truthy(foo)` and `truthy()` itself figures
#' out which method to actually dispatch.
truthy <- function(xx,...) UseMethod('truthy');
truthy.logical <- function(xx,...) xx;
truthy.factor <- function(xx,...) truthy.default(as.character(xx),...);
truthy.numeric <- function(xx,...) xx>0;
truthy.default <- function(xx,truewords=c('TRUE','true','Yes','T','Y','yes','y')
                           ,...) xx %in% truewords;

countfrac <- function(xx,outcomes
                      # set to TOTAL in order to do _only_ the total
                      ,groupcols='rai_range',sortby=groupcols
                      # set to 'none' if don't want to sort
                      ,dir=c('desc','asc','none')
                      ,summfns=c(n=sum,frac=mean)
                      # set to '' to disable
                      ,totalrow='Total'){
  # callback, so this function can call itself even if it's renamed
  # or nameless, for purposes of rerunning to get the totals row
  if(totalrow!='') thisfn <- sys.function();
  oo <- xx;
  # is the function being invoked for the purpose of calculating a total row?
  if(doingtotalrow <- tolower(trimws(groupcols[1]))=='total'){
    oo[[groupcols]] <- totalrow;
  };
  oo <- group_by_(oo,groupcols);
  # we coerce everything to logical, whatever it might have been
  # but leave alone columns other than those specified in the 'outcomes' variable
  oo[,outcomes] <- mutate_all(oo[,outcomes],truthy);
  # Two different summary tables, the second one applies the same set of
  # functions to everything, and the first one is just a count that doesn't 
  # rely on any one specific column. So we do them separately and then cbind()
  oo <- cbind(summarise(oo,bin_n=n()),summarise_all(oo[,c(groupcols,outcomes)],summfns)[,-1]) %>%
    # creating cumul_count, as in the current code
    mutate(cumul_count=rev(cumsum(rev(bin_n))));
  # sort, if desired, just as in the current code
  if(!doingtotalrow) oo <- switch(match.arg(dir)
                                  ,none=oo
                                  ,desc=arrange_(oo,sprintf('desc(%s)',groupcols))
                                  ,asc=arrange_(oo,groupcols));
  # now we insure that the column order is the same as the order of the groupcols
  # argument-- first the grouping column name, then the 'bin_n' (used to be called
  # rai_n, but this isn't specific to RAI, could group by anithing, so renamed)
  # cumul_count...
  oo <- oo[,c(groupcols,'bin_n','cumul_count'
              # ...and then this: it pastes the suffixes as obtained from the names 
              # of the summfns argument but orders them in the same way as they were
              # given, rather than first the first summary function and then the second
              ,c(t(outer(outcomes,names(summfns),paste,sep='_'))))];
  if(!doingtotalrow&&totalrow!='') { 
    tot <- thisfn(xx,outcomes,groupcols=totalrow,summfns=summfns,totalrow='');
    names(tot)<-names(oo);
    tot[[groupcols]] <- totalrow;
    oo <- rbind(oo,tot);
  }
  attr(oo,'call') <- match.call();
  oo;
}

#' ### Specific to RAI-A
#' 

#' Return a tableone object formatted just the way we like it
#' @param xx     A \code{data.frame} (required).
#' @param vars   Vector of variable names to pass to \code{CreatTableOne} (optional)
#' @param ...    Named expressions from which to create the strata variable
#'               for \code{CreatTableOne} (only tested for one variable)
stratatable <- function(xx,vars=NULL,...){
  nmx <- names(xx);
  xx <- transform(xx,...);
  mystrata <- setdiff(names(xx),nmx);
  res <- tableone::CreateTableOne(data=xx, vars= vars, strata=mystrata) %>% print;
  # res[,'p'] %>%  gsub("<","", x=.) %>% trimws() %>%
  #   as.numeric() %>% p.adjust() %>% cbind(res,padj=.) %>%
  #   data.frame %>% dplyr::select(-test) -> res;
  res[,'p'] %>%  gsub("<","", x=.) %>% trimws() %>%
    as.numeric() %>% p.adjust() %>% cbind(res[,1:2],padj=.) -> res;
  return(res);
}


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
v <- function(var,dictionary=dct0) {cc<-substitute(var);na.omit(dictionary[dictionary[[as.character(cc)]],'colname'])}
