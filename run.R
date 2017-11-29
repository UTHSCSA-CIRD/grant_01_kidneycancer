#' ---
#' title: "Kidney Cancer in Cancer Registry"
#' author: "Bokov"
#' date: "11/27/2017"
#' ---
#' 
source('global.R');

#' ## Load data if it exists 
#' 
#' (useful later, right now don't bother saving sessions)
#'if(session %in% list.files()) load(session);
#' Load your data. Notice that we're using `read_csv()` from the readr library.
#' It is a little smarter than the built-in `read.csv()`
dat0 <- read_csv(inputdata,na=c('(null)',''));
#' Read in the data dictionary
dct0 <- read_csv(dctfile,na = '');
dct0$naaccr <- dct0$rule=='naaccr';
dct1 <- read_tsv(tcrcodes,na='');
#colnames(dat0) <- tolower(colnames(dat0));
#' Create copy of original dataset
dat1 <- with(dat0,dat0[order(patient_num,age_at_visit_days),]) %>% 
  group_by(patient_num) %>% 
  # record for each patient their age at first diagnosis
  mutate(a_aad=age_at_visit_days[which(!is.na(v001_Dt_of_Dgns))[1]]);
#' Get rid of the extra quotes that somehow got into data
sapply(intersect(names(dat1),v(removequotes)),function(ii) .GlobalEnv$dat1[[ii]] <- gsub('\"','',dat1[[ii]]));
#' Remove the non-informative or redundant columns
dat1[,v(remove)] <- NULL;
#' Turn NA/code columns into T/F
.junk <- sapply(intersect(names(dat1),v(factor)),function(ii) .GlobalEnv$dat1[[ii]] <- factor(dat1[[ii]]));
dat1$v029_Mrtl_Sts <- gsub('DEM|MARITAL:','',dat1$v029_Mrtl_Sts,fixed = T);
#' For the diabetes and BMI, drop the values that occur after diagnosis
for(ii in c(v(trimlate),'start_date')) {
  dat1 <- mutate_(dat1
                  ,.dots=setNames(
                    list(parse(text=sprintf('ifelse(age_at_visit_days>a_aad,NA,%s)',ii))[[1]]),ii));
}
#' ## When multiple codes, get rid of the unknown indicators
dat1$v011_Grd <- gsub(",9|9,",'',dat1$v011_Grd);
dat1$v021_TNM_Cln_Stg_Grp <- gsub(",9|9,",'',dat1$v021_TNM_Cln_Stg_Grp);
dat1$v022_TNM_Cln_M <- gsub(",88|88,",'',dat1$v022_TNM_Cln_M);

dat2 <- summarize_all(dat1,function(xx) tail(c(NA,na.omit(xx)),1));
.junk <- sapply(intersect(names(dat2),v(tf)),function(ii) .GlobalEnv$dat2[[ii]]<- !is.na(dat2[[ii]]));
dat2$start_date <- as.Date(dat2$start_date);
dat3 <- subset(dat2
               ,start_date>as.Date('2008-01-01') &
                 start_date<as.Date('2013-12-31'));
dat3[,c('birth_date','age_at_visit_days')] <- NULL;

#' ## Look at which codes have more than one value
sapply(intersect(v(naaccr),names(dat3)),function(xx) class(dat3[[xx]])) -> foo;
sapply(grep('Cmrbd_Cmplctn',names(foo)[foo=='character'],invert = T,val=T),function(ii) any(grepl(",",dat3[[ii]]))) -> bar;
sapply(names(bar[bar]),function(jj) table(grep(',',dat3[[jj]],val=T)));

dct0a <- subset(dct0,colname %in% names(dat3));

write_tsv(dat3,na='',path='local/out/HSC20170563N_171128_kcTCR_data.csv');
write_tsv(dct0a,na='',path='local/out/HSC20170563N_171128_kcTCR_columns.csv');

#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX) %in% c('20,01','20,62,02','10','20','21','60,10')]<-'Private';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX) %in% c('31','35','64','64,02','31,35')] <-'Medicaid';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX) %in% c('60','61','62','63','60,02','60,62')] <- 'Medicare';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='01'] <- 'Uninsured';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='02'] <- 'Self-Pay';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='65'] <- 'TriCare';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='67'] <- 'VA';
#' levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='99'] <- 'Unknown';
#' 
#' levels(dat1$v008_Spnsh_Hspnc)[levels(dat1$v008_Spnsh_Hspnc) == '0'] <- 'Non Hispanic'
#' levels(dat1$v008_Spnsh_Hspnc)[levels(dat1$v008_Spnsh_Hspnc) %in% c('9','0,9')] <- 'Unknown'
#' levels(dat1$v008_Spnsh_Hspnc)[!levels(dat1$v008_Spnsh_Hspnc) %in% c('Non Hispanic','Unknown')] <- 'Hispanic'
#' 
#' possiblyblank <- c('v000_Cmrbd_Cmplctn','v002_Brthplc','v004_Mrtl_Sts_DX','v005_Rc','v006_Rc','v008_Spnsh_Hspnc','v009_Sx','v010_Cmrbd_Cmplctn','v011_Grd','v012_Cmrbd_Cmplctn','v013_Prmr_Pr_at_DX','v014_Cls_of_Cs','v015_Cmrbd_Cmplctn','v016_TNM_Cln_N','v017_SR_Smr_Stg','v018_TNM_Cln_T','v021_TNM_Cln_Stg_Grp','v022_TNM_Cln_M','v023_TNM_Edtn_Nmbr','v025_Cmrbd_Cmplctn','v026_Cmrbd_Cmplctn','v027_Cmrbd_Cmplctn','v028_Cmrbd_Cmplctn','v029_Mrtl_Sts','v030_Cmrbd_Cmplctn','v031_Cmrbd_Cmplctn','v032_Bd_Ms_Indx_num','v032_Bd_Ms_Indx_info');
#' cbind(apply(dat1[,possiblyblank],1,function(xx) all(is.na(xx))),apply(dat1[,v(tf)],1,function(xx) !any(xx))) %>% 
#'   apply(1,all) -> emptyrows;
#' dat2 <- dat1[!emptyrows,];
#' #' TODO: exclude all non-static EMR facts happening at a greater age than age at first diagnosis
#' #' TODO: then, keep the non-static EMR fact that is most recent
