#' ---
#' title: "Kidney Cancer in Cancer Registry"
#' author: "Bokov"
#' date: "11/27/2017"
#' ---
#' 
source('global.R');

#' ## Define metadata
cls_tf <- c(
'_Spnsh_Hspnc'
,'_Mlgnt_nplsm'
#,'_Mlgnt_nplsm_inactive'
,'_Vtl_Sts'
,'_Mlgnt_nplsm'
#,'_Mlgnt_nplsm_inactive'
,'_Dbts_mlts'
#,'_Dbts_mlts_inactive'
,'_Dcsd_pr_SS'
,'_Dbts_mlts'
#,'_Dbts_mlts_inactive'
,'_Hspnc_or_Ltn'
,'_Kdn_NS');

finclass_rename <- rbind(
  c("101","Carelink"),
  c("102","Governmental"),
  c("103","Contract"),
  c("10","Blue Shield"),
  c("12","Other"),
  c("1","Commercial"),
  c("2","Medicare"),
  c("3","Medicaid"),
  c("4","Self-pay"),
  c("5","Worker's Comp"),
  c("0","Financial Class Not Recorded")
  );

cls_forlinkage <- c('patient_num','start_date','birth_date','sex_cd'
                    ,'language_cd','race_cd','_Dcsd_pr_SS'
                    ,'_Cs_of_Dth','_Vtl_Sts','_Brthplc','_Sx'
                    ,'_Spnsh_Hspnc','_Mrtl_Sts_DX','_TNM_Cln_Dscrptr','_Rc'
                    ,'_Kdn_NS');
cls_canceremr <- '_Mlgnt_nplsm';
cls_comorb <- '_Cmrbd_Cmplctn';

#' locators for which the summarize function should be 'any'
cls_summ_any <- 'a_diabetes';

pnum <- 'patient_num'; time <- 'age_at_visit_days';



#' ## Load data if it exists 
#' 
#' (useful later, right now don't bother saving sessions)
#'if(session %in% list.files()) load(session);
#' Load your data. Notice that we're using `read_csv()` from the readr library.
#' It is a little smarter than the built-in `read.csv()`
dat0 <- group_by(read_csv(inputdata,na=c('(null)','')),patient_num);
#' Read in the data dictionary
dct0 <- read_csv(dctfile,na = '') %>% subset(colname %in% names(dat0));
#' add information about the standard columns all datafinisher
#' files have which datafinisher currently does not include in
#' its meta file
dct0 <- rbind(data.frame(colname=setdiff(names(dat0),dct0$colname)
                         ,colname_long=setdiff(names(dat0),dct0$colname)
                         ,rule='demog',stringsAsFactors = F),
              dct0);
#' add locator information, since datafinisher doesn't yet do it
#' (a locator is a stable, but non unique, part of a column name
#' useful for manipulating groups of columns without knowing what
#' their exact names will be in future revisions of the dataset)
#' All it is is removing the vXXX prefix where it exists
dct0$locator <- gsub('v[0-9]{3}','',dct0$colname);
dct0$class <- sapply(dat0,function(xx) class(xx)[1])
dct0$nunique <- sapply(dat0,function(xx) length(unique(na.omit(xx))));
dct0$nna <- sapply(dat0,function(xx) sum(is.na(xx)));
dct0$c_naaccr <- dct0$rule=='naaccr';
dct0$c_extra <- grepl('_inactive$|_info$',dct0$colname);
dct0$c_notanalytic <- grepl('_inactive$|_info$|^patient_num$',dct0$colname);
dct0$c_date <- dct0$class %in% c('POSIXct','Date');
dct0$c_tf <- with(dct0,(nunique==1|rule=='diag')&!c_notanalytic);
dct0$c_canceremr <- dct0$locator %in% cls_canceremr;
dct0$c_char <- dct0$class == 'character';
dct0$c_dm <- dct0$locator == '_Dbts_mlts';
dct0$c_rcc <- dct0$locator == '_Mlgnt_nplsm';
dct0$c_fin <- dct0$locator == '_Fncl_Cls';
if('naaccr' %in% dowhat){
  dct0$c_todrop <- dct0$locator %in% '_Hspnc_or_Ltn';
} else if('epic00' %in% dowhat){
  dct0$c_todrop <- with(dct0, locator %in% c('_Dcsd_pr_SS') | 
    c_dm | c_rcc | grepl('_inactive$|_info$', dct0$colname));
}
lkup <- setNames(dct0$colname,dct0$colname_long);
dct1 <- read_tsv(tcrcodes,na='');
dat1 <- dat0[order(dat0[[pnum]],dat0[[time]]),];
#' Remove those annoying quote marks
for(ii in v(c_char)) dat1[[ii]] <- gsub('"','',dat1[[ii]]);
if('naaccr' %in% dowhat){
  # drop the non-naaccr fields that are only there as a workaround placeholder
  dat1[,v(c_todrop)] <- NULL;
  # drop the rows that were only there because of the non-naaccr fields
  dat1 <- dat1[apply(dat1[,-(1:8)],1,function(xx) any(!is.na(xx))),];
  dct0 <- subset(dct0,colname %in% names(dat1));
  var_tt <- subset(dct0,locator %in% '_Dt_Lst_Cntct')$colname[1];
}

#' For naaccr data
if('naaccr' %in% dowhat){
  #' Mark cases where diabetes is listed as any of the comorbidities
  dat1$a_diabetes <- apply(dat1[,subset(dct0,locator %in% cls_comorb)$colname]
                           ,1,function(xx) any(grepl('^250|,250',xx)));
  dat2 <- summarise_all(dat1, function(xx) last(na.omit(xx)));
  dat2[,c('a_diabetes','t0','t1','age0','tt','cc')] <- summarise_(dat1
                                    ,a_diabetes='any(a_diabetes)'
                                    ,t0='min(v018_Dt_of_Dgns,na.rm=T)'
                                    ,t1="max(v022_Dt_Lst_Cntct,na.rm=T)"
                                    ,tt='max(as.numeric(t1 - t0,units="days"),na.rm=T)'
                                    ,tt='if(is.infinite(tt)) NA else tt'
                                    ,age0='max(as.numeric(t0 - v012_Dt_of_Brth,units="days"),na.rm=T)'
                                    ,age0='if(is.infinite(age0)) NA else age0'
                                    ,cc="any(gsub('1,0','1',v024_Vtl_Sts)=='0',na.rm=T)"
                                    )[,c('a_diabetes','t0','t1','age0','tt','cc')];
  dat2$v024_Vtl_Sts <- gsub('1,0','1',dat2$v024_Vtl_Sts);
  # the for-export, human-readable tables
  dat3 <- mapnames(dat2,lkup)[-(2:8)];
  dct03 <- dct0[-(2:8),];
} else if('epic00' %in% dowhat) {
  dat1$a_Diabetes <- apply(dat1[,v(c_dm)],1,function(xx) any(!is.na(xx)));
  dat1$a_Kidney_Cancer <- apply(dat1[,v(c_rcc)],1,function(xx) any(!is.na(xx)));
  dat1[,v(c_todrop)] <- NULL;
  # updated dictionary that matches the remaining columns in dat1
  # naming dat3 for consistency with naaccr case
  dct03 <- dct01 <- subset(dct0,colname %in% names(dat1));
  formals(.GlobalEnv$v)$dictionary <- quote(dct01);
  # convert indicator variables
  for(ii in v(c_tf)) dat1[[ii]] <- !is.na(dat1[[ii]]);
  # trim down the financial class variables
  dat1[[v(c_fin)]] <- gsub('[^0-9,]','',dat1[[v(c_fin)]]);
  # human-readable labels for financial class
  browser();
  dat1$a_Finclass <- submulti(dat1[[v(c_fin)]],finclass_rename);
  # for export, naming dat3 for consistency with the naaccr case
  dat3 <- setNames(dat1,gsub(' ','_',submulti(names(dat1),dct01[
    ,c('colname','colname_long')])));
}
.savefile <- gsub('.csv$','.out.tsv',inputdata);
write_tsv(dat3,path=.savefile);
.savefile <- gsub('.csv$','.datadict.tsv',datadct);
write_tsv(dct03,path=.savefile);

#colnames(dat0) <- tolower(colnames(dat0));
#' Create copy of original dataset
#dat1 <- with(dat0,dat0[order(patient_num,age_at_visit_days),]) %>% 
#  group_by(patient_num); #%>% 
#dat1[,v(c_tf)] <- lapply(dat1[,v(c_tf)],function(xx) !is.na(xx));
#dat1$a_canceremr <- apply(dat1[,v(c_canceremr)],1,any);

#' # Starting here is old stuff, not necessarily usable for the new dataset

#' ## Dataset for Linking
#' 
#' Linking = 
#' * deathdates from UHS/SSDMF
#' * annotating NAACCR codes,
#' * adding EPIC IDs
#' * adding addresses, zipcodes, cities
dat2 <- dat0[,dct0$colname[dct0$locator %in% cls_forlinkage]];
dat2 <- dat2[apply(dat2[,dct0$colname[dct0$locator %in% cls_forlinkage[-(1:6)]]],1
              ,function(xx) !all(is.na(xx))),];
dct20 <- subset(dct0,colname %in% colnames(dat2));
#' 
#' 
#' 
#' 
# record for each patient their age at first diagnosis
#mutate(a_aad=age_at_visit_days[which(!is.na(v001_Dt_of_Dgns))[1]]);
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
