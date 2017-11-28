#' ---
#' title: "Kidney Cancer in Cancer Registry Libraries and Global Settings"
#' author: "Bokov"
#' date: "11/27/2017"
#' ---
#' 
#' ## Load libraries
#+ warning=FALSE, message=FALSE
rq_libs <- c('compiler'                                   # just-in-time compilation
             ,'MatchIt','DHARMa'                          # propensity scores and glm residuals
             ,'pscl'                                      # zero-inflated poisson, sigh
             ,'survival','MASS','Hmisc','zoo','coin'      # various analysis methods
             ,'readr','dplyr','stringr','magrittr'        # data manipulation & piping
             ,'ggplot2','ggfortify','grid','GGally'       # plotting
             ,'stargazer','broom', 'tableone','janitor'); # table formatting
rq_installed <- sapply(rq_libs,require,character.only=T);
rq_need <- names(rq_installed[!rq_installed]);
if(length(rq_need)>0) install.packages(rq_need,repos='https://cran.rstudio.com/',dependencies = T);
sapply(rq_need,require,character.only=T);
#' Turn JIT to max: pre-compile all closures, `for`, `while`, and `repeat` loops
enableJIT(3);
#' ## Load local config file
#' 
source('./config.R');
#source('./metadata.R'); 
source('./functions.R');



#' ## Set generic variables
#' 
#' data dictionary:
dctfile <- datadct;
tcrcodes <- tcrcodes;
#' saved session data (not used right now)
session <- 'session.rdata';
