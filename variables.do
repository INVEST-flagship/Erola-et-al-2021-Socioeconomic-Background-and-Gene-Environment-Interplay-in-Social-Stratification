*** Socioeconomic Background and Gene–Environment Interplay in Social Stratification across the Early Life Course 

*** This do file defines all the used variables 
*** Run this do-file before all the models 
*** date:  2.1.2020
*** Hannu Lehti 



*Master file 
do "W:\Hannu L\twin correlations\twin 1\01_do\master.do"
*Import data 
use "$data\twinfinal10122019", clear

replace isei_max =. if isei_max==9999 // 1,189 no isei for twins // 1,571 missing twin and siblings total
replace isei_mean =. if isei_max==9999 
replace inc = . if (inc == 50  | inc==0) // no inc 523 // total isei inc missing = 651 after 

codebook inc isei_max
sum id if inc==. & isei_max==. // both missing ISEI and inc 519 total 
sum id if inc==. & isei_max==. & twin==1 // both missing ISEI and inc 405 for twins 

order isei_max inc eduy, alpha 
egen omiss=rowmiss(isei_max inc eduy)
replace omiss = 1 if omiss > 1

forvalues t = 1/5 { // 278 missing total 
replace pedu`t' = . if pedu`t'==9999 
replace faminc`t' =. if faminc`t'==50  | faminc`t'==0
replace pisei`t' = . if pisei`t' ==9999
}

order pedu1-pedu5 faminc1-faminc5 pisei1-pisei5, alpha
egen pmiss = rowmiss(pedu1-pedu5 faminc1-faminc5 pisei1-pisei5)
replace pmiss = 1 if pmiss>1
ta pmiss if twin==1

*omit missing values 
mark nomiss 
markout nomiss inc eduy isei_max pedu? faminc? pisei? id famid
keep if nomiss==1
drop total

bysort famid: gen total= _N 
bysort famid: egen twintotal = total(twin) if twin==1

ta total
ta twintotal

*keep twins only 
keep if twintotal == 2 
keep if twin==1

*std. education 
gen std_eduy=.
forvalues t=0/1 { 
*edu at 28
sum eduy if female==`t' & twin==1
replace std_eduy = (eduy-r(mean))/r(sd) if female==`t'
}
bys female: sum std_eduy if twin==1

*std. isei 
gen std_isei=.
forvalues t=0/1 { 
sum isei_max if female==`t' & twin==1
replace std_isei = (isei_max-r(mean))/r(sd) if female==`t' & twin==1
}
bys female: sum std_isei if twin==1

*std. log-income  
gen linc = log(inc)
gen std_inc=.
forvalues t=0/1 { 
sum linc if female==`t' & twin==1
replace std_inc = (linc-r(mean))/r(sd) if female==`t' & twin==1
}
bys female: sum std_isei if twin==1

*same sex correction 
ta same_sex if twin==1,matcell(x)
matrix list x 
scalar os = x[1,1]
scalar ss = x[2,1]
scalar sscor = (ss-os)/ss + 0.5*(os/ss)
dis sscor

*log family income and parental education into dummies
forvalues t = 1/5 { 
gen lfaminc`t' = log(faminc`t') 

foreach x in p m f {
recode `x'edu`t' (6/9=6),  copyrest
replace `x'edu`t' = . if `x'edu`t' ==9999 
ta `x'edu`t', gen(`x'edu`t')
	
	}	
}

sum id famid isei_max inc eduy famid id pedu? faminc? pisei?
