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
dct1 <- read_tsv(tcrcodes,na='');
#colnames(dat0) <- tolower(colnames(dat0));
#' Create copy of original dataset
dat1 <- dat0;
#' Get rid of the extra quotes that somehow got into data
sapply(intersect(names(dat1),v(removequotes)),function(ii) .GlobalEnv$dat1[[ii]] <- gsub('\"','',dat1[[ii]]))
#' Remove the non-informative or redundant columns
dat1[,v(remove)] <- NULL;
#' Turn NA/code columns into T/F
sapply(intersect(names(dat1),v(tf)),function(ii) .GlobalEnv$dat1[[ii]]<- !is.na(dat1[[ii]]));
sapply(intersect(names(dat1),v(factor)),function(ii) .GlobalEnv$dat1[[ii]] <- factor(dat1[[ii]]));
 
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX) %in% c('20,01','20,62,02','10','20','21','60,10')]<-'Private';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX) %in% c('31','35','64','64,02','31,35')] <-'Medicaid';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX) %in% c('60','61','62','63','60,02','60,62')] <- 'Medicare';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='01'] <- 'Uninsured';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='02'] <- 'Self-Pay';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='65'] <- 'TriCare';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='67'] <- 'VA';
levels(dat1$v013_Prmr_Pr_at_DX)[levels(dat1$v013_Prmr_Pr_at_DX)=='99'] <- 'Unknown';

levels(dat1$v008_Spnsh_Hspnc)[levels(dat1$v008_Spnsh_Hspnc) == '0'] <- 'Non Hispanic'
levels(dat1$v008_Spnsh_Hspnc)[levels(dat1$v008_Spnsh_Hspnc) %in% c('9','0,9')] <- 'Unknown'
levels(dat1$v008_Spnsh_Hspnc)[!levels(dat1$v008_Spnsh_Hspnc) %in% c('Non Hispanic','Unknown')] <- 'Hispanic'

possiblyblank <- c('v000_Cmrbd_Cmplctn','v002_Brthplc','v004_Mrtl_Sts_DX','v005_Rc','v006_Rc','v008_Spnsh_Hspnc','v009_Sx','v010_Cmrbd_Cmplctn','v011_Grd','v012_Cmrbd_Cmplctn','v013_Prmr_Pr_at_DX','v014_Cls_of_Cs','v015_Cmrbd_Cmplctn','v016_TNM_Cln_N','v017_SR_Smr_Stg','v018_TNM_Cln_T','v021_TNM_Cln_Stg_Grp','v022_TNM_Cln_M','v023_TNM_Edtn_Nmbr','v025_Cmrbd_Cmplctn','v026_Cmrbd_Cmplctn','v027_Cmrbd_Cmplctn','v028_Cmrbd_Cmplctn','v029_Mrtl_Sts','v030_Cmrbd_Cmplctn','v031_Cmrbd_Cmplctn','v032_Bd_Ms_Indx_num','v032_Bd_Ms_Indx_info');
cbind(apply(dat1[,possiblyblank],1,function(xx) all(is.na(xx))),apply(dat1[,v(tf)],1,function(xx) !any(xx))) %>% 
  apply(1,all) -> emptyrows;
dat2 <- dat1[!emptyrows,];
#' TODO: exclude all non-static EMR facts happening at a greater age than age at first diagnosis
#' TODO: then, keep the non-static EMR fact that is most recent
#' ## Column names of primary relevance
c_modelvars <- c(v(c_rai),v(c_postop),'income_final','hispanic_ethnicity');

#' Backup up the modified cnopatos column
#' ...because it's easier if patos-subtracted columns are modified in place
c_canbepatos <- v(c_canbepatos);
c_patos <- v(c_patos);
dat1[,paste0('bak_',c_canbepatos)] <- dat1[,c_canbepatos];
#' Since the postop & accompanying columns are counts, we just subtract patos
#' from their postop counterparts.
#dat1[,c_canbepatos] <- mapply(function(xx,yy){ifelse(xx,0,yy)}, dat1[,carepatos], dat1[,cnopatos]);
dat1[,c_canbepatos] <- dat1[,c_canbepatos] - dat1[,c_patos];


#' Create binned versions of certain numeric vars.
#' (commented out until we can put a c_num2bin or something into dct0)
# dat1[,paste0('bin_',cnum2bin)] <- sapply(dat1[,cnum2bin],function(ii){
#   qii <- c(0,quantile(ii,c(.25,.5,.75),na.rm=T),Inf);
#   cut(ii,breaks = qii);
# })

#' ## Create response variables
#' 
#' Create a column that is sum of all complications. Lets name analytically
#' created colums with an `a_` prefix. This way you could make it binary
#' via `dat1$a_allcomp > 0` or leave it as an integer and use number of 
#' different complications as a proxy for severity of outcome.
#dat1$a_allcomp <- rowSums(dat1[,csrscomp]);
c_postop_yesno <- setdiff(v(c_postop),v(c_count));
c_postop_count <- intersect(v(c_postop),v(c_count));
dat1$a_postop <- rowSums(dat1[,c_postop_count]) +
  apply(dat1[,c_postop_yesno],1,function(xx) sum(na.omit(xx %in% c('Yes','Positive'))));

#' -Hack the values of these variables to be binary for now.-
#dat1$sepsis_sirs_sepsis_sepshk_48h <- dat1$sepsis_sirs_sepsis_sepshk_48h != 'None';
#dat1$first_unp_ret_or <- dat1$first_unp_ret_or == 'Yes';
#dat1$hisp <- dat1$hispanic_ethnicity == 'Yes';
c_cd4_yesno <- setdiff(v(c_cd4),v(c_count));
c_cd4_count <- intersect(v(c_cd4),v(c_count));


#' Do the same as above but just for the `ccd4` complications
#dat1$a_cd4 <- rowSums(dat1[,c_cd4]);
dat1$a_cd4 <- rowSums(dat1[,c_cd4_count]) + 
  apply(dat1[,c_cd4_yesno],1,function(xx) sum(na.omit(xx=='Yes')));

#' TRUE/FALSE variables for postop and cd4 in both cases
#' indicating whether or not the patient had _any_ complications
dat1$a_any_cd4 <- factor(dat1$a_cd4>0);
dat1$a_any_postop <- factor(dat1$a_postop>0);

dat1$a_transfer <- dat1$origin_status!='Not transferred (admitted from home)';
dat1$a_readm_30_dy <- ifelse(dat1$readm_30_dy==0, 'FALSE', 'TRUE'); 

#' Time from surgery to adverse outcome if any
dat1$a_t <- with(dat1
                 ,pmin(
                   dt_death
                   ,dt_first_readm
                   ,dt_first_unpl_ret_or
                   ,dt_second_unp_ret,na.rm = T) %>% 
                   difftime(proc_surg_finish,units='days') %>%
                   as.numeric());
#' Censor the variables at 30 days
dat1$a_t[dat1$a_t>30] <- 30;
dat1$a_t[is.na(dat1$a_t)] <- 30;
dat1$a_c <- dat1$a_t!=30;


#' Obtain the RAI score
dat1$a_rai <- raiscore(dat1);
dat1$rai_range <- cut(dat1$a_rai, breaks=c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55));
dat1$a_discrete_rai <- cut(dat1$a_rai
                           ,c(0,15,21,Inf)
                           ,right = F
                           ,labels = c('Non Frail','Pre Frail','Frail'));
dat1$a_rai_hisp <- with(dat1,interaction(a_discrete_rai,hispanic_ethnicity));


#' ## The Rockwood Scale
dat1$a_rockwood <- with(dat1,(
  as.numeric(bmi>=25)+
    as.numeric(origin_status!='Not transferred (admitted from home)')+
    as.numeric(diabetes_mellitus!='No')+
    as.numeric(current_smoker_within_1_year=='Yes')+
    as.numeric(dyspnea!='No')+
    as.numeric(!functnal_heath_status%in%c('Independent','Unknown'))+
    as.numeric(vent_dependent=='Yes')+
    as.numeric(history_severe_copd=='Yes')+
    as.numeric(ascites_30_dy_prior_surg=='Yes')+
    as.numeric(hypertensn_req_medicatn=='Yes')+
    as.numeric(acute_renal_failure=='Yes')+
    as.numeric(currently_dialysis=='Yes')+
    as.numeric(disseminated_cancer=='Yes')+
    as.numeric(open_wound=='Yes')+
    as.numeric(steroid_immunosupp=='Yes')+
    as.numeric(x_loss_bw_6_months_prior_surg=='Yes')+
    as.numeric(bleeding_disorder=='Yes')+
    as.numeric(chr_30_dy_prior_surg=='Yes')+
    as.numeric(isTRUE(serum_creatinine>3))
)/(
  19-
    as.numeric(functnal_heath_status=='Unknown')-
    as.numeric(is.na(serum_creatinine))
));

c_tabsievars <- c(v('c_tabsie')
                ,'a_postop','a_any_postop','a_cd4','a_any_cd4'
                ,'a_rai','a_discrete_rai','a_rockwood');
#' ## Transform Rows
#'
#' ### Sort the rows by patient ID and then by date of surgery, ascending
#' 
dat1 <- dat1[order(dat1$proc_surg_start),];
#' 
#' ### Drop patients without an income
#dat1 <- dat1[!is.na(dat1$income_final),];

#' ### Make several subsets of dat1 all at once
#' 
#' for later use to make multiple versions of the same table and multiple
#' versions of the same graph, as for item #3 of the 10/13/2017 PKS email.

 dat1subs <- ssply(dat1
                  ,full=T
		  ,all_elective=elective_surg=='Yes'
		  ,all_urgent=elective_surg=='No' & emergency_case=='No'
		  ,all_emergency=emergency_case=='Yes'
                  ,all_colon_all=cpt_code %in% v(c_all_colon,dct1)
		  ,all_colon_elective=cpt_code %in% v(c_all_colon,dct1) &
			 elective_surg=='Yes'
		  ,all_colon_urgent=cpt_code %in% v(c_all_colon,dct1) &
			 elective_surg=='No' & emergency_case=='No'
		  ,all_colon_emergency=cpt_code %in% v(c_all_colon,dct1) &
			 emergency_case=='Yes'
		  ,open_colon_all=cpt_code %in% v(c_open_colon,dct1)
		  ,open_colon_elective=cpt_code %in% v(c_open_colon,dct1) &
			 elective_surg=='Yes'
		  ,open_colon_urgent=cpt_code %in% v(c_open_colon,dct1) &
			 elective_surg=='No' & emergency_case=='No'
		  ,open_colon_emergency=cpt_code %in% v(c_open_colon,dct1) &
			 emergency_case=='Yes'
		  ,lapa_colon_all=cpt_code %in% v(c_lapa_colon,dct1)
		  ,lapa_colon_elective=cpt_code %in% v(c_lapa_colon,dct1) &
			 elective_surg=='Yes'
		  ,lapa_colon_urgent=cpt_code %in% v(c_lapa_colon,dct1) &
			 elective_surg=='No' & emergency_case=='No'
		  ,lapa_colon_emergency=cpt_code %in% v(c_lapa_colon,dct1) &
			 emergency_case=='Yes'                  
);


#' ### Create a version of the dataset that only has each patient's 1st encounter
#' 
#' (you need to have specified the name of the ID column in `metadata.R`)
dat2 <- group_by(dat1,idn_mrn) %>% summarise_all(first);

#' Filter down to only NHW and hispanic
dat3<-subset(dat2,hispanic_ethnicity!='Unknown'&(hispanic_ethnicity=='Yes'|race=='White'));
dat3$hispanic_ethnicity<-factor(dat3$hispanic_ethnicity);
#'  creating tables similar to the tables that Dan MacCarthy creates for the VASQIP data:


#' Both of our main non time-to-event responses
c_resps <- c('a_postop','a_cd4');

#' ### Create your random sample
source('random_seed.R');
pat_samp <- sample(dat3$idn_mrn,1000,rep=F);
dat4 <- subset(dat3,idn_mrn %in% pat_samp);

