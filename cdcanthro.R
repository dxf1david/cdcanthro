
set_cols_first <- function (DT, cols, intersection = TRUE) # thanks to hutils
   {
       if (intersection) {
           return(setcolorder(DT, c(intersect(cols, names(DT)),
               setdiff(names(DT), cols))))
       }
       else {
           return(setcolorder(DT, c(cols, setdiff(names(DT), cols))))
       }
   }

z_score=function(var, l, m, s){ # LMS formula with modified (m) z-scores
      ls=l*s; invl=1/l
      z = (((var/m) ^ l) -1) / (ls) # z-score forumula
      sdp2 = (m * (1 + 2*ls) ^ (invl)) - m; # modified z-score (+2)
      sdm2 = m - (m * (1 - 2*ls) ^ (invl));
      mz=fifelse(var < m, (var - m)/(0.5*sdm2),
                   (var - m)/(sdp2*0.5) )
      list(z, mz)
   }

cdcanthro <- function(data, age=age_in_months, wt=weight_kg, ht=height_cm, bmi=bmi)
{
   age_in_months <- weight <- height <- seq_ <- sex <- agey <- bz <- bp <-
      lwt2 <- mwt2 <- swt2 <- lbmi2 <- mbmi2 <- sbmi2 <- lht2 <- mht2 <- sht2 <-
      lwt1 <- mwt1 <- swt1 <- lbmi1 <- mbmi1 <- sbmi1 <- lht1 <- mht1 <- sht1 <-
      mbmi <- lbmi <- sbmi <- mref <- sref <- denom <- weight_kg <- height_cm <-
      bmiz <- l <- m <- s <- waz <- haz <- z1 <- z0 <- p95 <- bmip <-
      ebp <- ebz <- agemos <- agemos1 <- agemos2 <- NULL

   setDT(data)
   data$seq_ <- 1L:nrow(data) # needed for merging back with original data
   dorig <- copy(data)

   age <- deparse(substitute(age)); data$age <- data[[age]]
   wt <- deparse(substitute(wt)); data$wt <- data[[wt]]
   ht <- deparse(substitute(ht)); data$ht <- data[[ht]]
   bmi <- deparse(substitute(bmi)); data$bmi <- data[[bmi]]

   data <- data[between(age,24,240) & !(is.na(wt) & is.na(ht)),
                    .(seq_, sex,age,wt,ht,bmi)];
   v1 <- c('seq_','id','sex','age','wt','ht','bmi')

   # 'dref' is CDCref_d.csv,  https://www.cdc.gov/nccdphp/dnpao/growthcharts/resources/sas.htm
   dref <- ref_data[`_AGEMOS1`>23 & denom=='age']
   names(dref) <- tolower(names(dref))
   names(dref) <- gsub('^_', '', names(dref))

   d20 <- dref[agemos2==240,
            .(sex,agemos2,lwt2,mwt2,swt2,lbmi2,mbmi2,sbmi2,lht2,mht2,sht2)]
   names(d20) <- gsub('2','',names(d20));

   dref <- dref[,.(sex,agemos1,lwt1,mwt1,swt1,lbmi1,mbmi1,sbmi1,lht1,mht1,sht1)]
   names(dref) <- gsub('1','',names(dref));

   dref=rbindlist(list(dref,d20))
   adj_bmi_met <- dref[agemos==240,.(sex,mbmi,sbmi)]
   setnames(adj_bmi_met,c('sex','mref','sref'))
   adj_bmi_met <- adj_bmi_met[,':=' (mref=round(mref,5), sref=round(sref,5))]

   dref <- dref[adj_bmi_met, on='sex']
   v=c('sex','age','wl','wm','ws','bl','bm','bs','hl','hm','hs','mref','sref');
   setnames(dref,v)

   # interpolate reference data to match each agemos in input data
   v=c('sex','age','wl','wm','ws','bl','bm','bs','hl','hm','hs','mref','sref')

   if (length(setdiff(data$age,dref$age))>0) {
   uages=unique(data$age); uages

   db <- dref[sex==1]
        fapp <- function(v,...)approx(db$age,v,xout=uages)$y
        db <- sapply(db[,..v],fapp)
   dg <- dref[sex==2]
        fapp <- function(v,...)approx(dg$age,v,xout=uages)$y
        dg <- sapply(dg[,..v],fapp)
   dref <- setDT(data.frame(rbind(db,dg)))
   }

   du <- unique(data[,.(sex,age)],by=c('sex','age'))
   dref <- dref[du, on=c('sex','age')]

   setkey(data,sex,age); setkey(dref,sex,age)
   dt <- dref[data];

   dt[,c('waz', 'mwaz'):= z_score(dt$wt, dt$wl, dt$wm, dt$ws)]
   dt[,c('haz', 'mhaz'):= z_score(dt$ht, dt$hl, dt$hm, dt$hs)]
   dt[,c('bz', 'mbz'):= z_score(dt$bmi, dt$bl, dt$bm, dt$bs)]

   setDT(dt);  setnames(dt,c('bl','bm','bs'),c('l','m','s'))
   dt[,c('wl','wm','ws','hl','hm','hs'):=NULL]

   dt[,':=' (
      bp=100*pnorm(bz),
      p50= m * (1 + l*s*qnorm(0.5))^(1 / l),
      p85= m * (1 + l*s*qnorm(0.85))^(1 / l),
      p95= m * (1 + l*s*qnorm(0.95))^(1 / l),
      p97= m * (1 + l*s*qnorm(0.97))^(1 / l),
      wp=100*pnorm(waz),  hp=100*pnorm(haz),

     # other BMI metrics -- PMID 31439056
      z1=((bmi/m) - 1) / s,  # LMS formula when L=1: ((BMI/M)-1)/S
      z0 = log(bmi/m)/s # LMS transformation with L=0, note these end in '0'
     )
     ][,':=' (
      dist = z1 * m * s, # unadjusted distance from median with L=1
      adist = z1 * sref * mref, # Adjusted (to age 20y) dist from median
      perc = z1 * 100 * s, # unadjusted %distance from median
      aperc = z1 * 100*sref, # adj %distance from median
      perc0 = z0 * 100 * s, # unadj. %distance from median with L=0
      aperc0 = z0 * 100* sref,  # adj %distance from median w L=0
      bmip95=100*(bmi/p95),
      obese=1L*(bmi>=p95),
      sev_obese=1L*(bmi>=1.2*p95)
   )]

   ## now create Extended z-score for BMI >=95th P
    dt[,':=' (ebz=bz, ebp=bp, agey=age/12)]
    dt[, sigma:=fifelse(sex==1, 0.3728 + 0.5196*agey - 0.0091*agey^2,
                               0.8334 + 0.3712*agey - 0.0011*agey^2)]
    dt[bp>=95, ebp:=90 + 10*pnorm((bmi - p95) / sigma)]
    dt[bp>=95, ebz:=qnorm(ebp/100)]
    dt[bp>99 & is.infinite(ebz), ebz:=8.21] # highest possible value is 8.20945

   x <- c('agey','mref','sref','sex','wt','ht','bmi'); dt[,(x):=NULL]

   setnames(dt,c('adist',  'aperc',  'bp',   'bz',  'mbz',     'mwaz',   'mhaz',
                 'ebp',    'ebz',    'l',    'm',    's', 'wp', 'hp',
                 'perc', 'perc0'),
                c('adj_dist','adj_perc','bmip','bmiz','mod_bmiz','mod_waz','mod_haz',
                  'ext_bmip','ext_bmiz', 'bmi_l','bmi_m','bmi_s','wap', 'hap',
                  'perc_median', 'log_perc_median')
   )

   v=c('seq_', 'bmiz', 'bmip', 'waz', 'wap', 'haz', 'hap', 'p50', 'p95', 'p97', 'bmip95',
       'mod_bmiz', 'mod_waz', 'mod_haz',  'ext_bmip', 'ext_bmiz', 'adj_dist', 'dist',
       'adj_perc', 'perc_median', 'log_perc_median')
   dt <- dt[,..v]
   # can also use dt <- dt[,select(.SD,a:d,j:m)]
   setkey(dt,seq_); setkey(dorig,seq_)
   dtot <- dt[dorig]
   set_cols_first(dtot,names(dorig))
   dtot[,seq_:=NULL]

   dtot[]
}

