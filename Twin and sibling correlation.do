*** Socioeconomic Background and Gene–Environment Interplay in Social Stratification across the Early Life Course 
*** Twin and sibling correlations (appendix)
*** date: 24.10.2019
*** Hannu Lehti 

**** NOTE! USE STATA 16 ***** 

do "W:\Hannu L\twin correlations\01_do\master.do"
use "$data\twinfinal10122019", clear // new data with siblings


*replace isei_max = r(mean) if isei_max==9999
replace isei_max =. if isei_max==9999
replace isei_mean =. if isei_mean==9999
replace inc = . if inc == 50 | inc== 0 // no inc 438 // total isei inc missing = 1627


forvalues t = 1/5 { // 278 missing total 
replace pedu`t' = . if pedu`t'==9999 
replace fedu`t' = . if fedu`t'==9999 
replace medu`t' = . if medu`t'==9999 
replace faminc`t' =. if faminc`t'==50  | faminc`t'==0
replace pisei`t' = . if pisei`t' ==9999
}

*omit missing values 
mark nomiss 
markout nomiss inc eduy isei_max pedu? faminc? pisei? id famid
keep if nomiss==1
drop total
bysort famid: gen total= _N 
bysort famid: egen twintotal = total(twin) if twin==1

ta total
ta twintotal

drop if total ==1 | twintotal==1 // otetaan perheet missä yksi tai enemmän kaksosia ja sisaruksia
ta twin // 2240 observation deleted due to missing inc or isei 
ta total

* keep only twins
keep if twin==1
drop total
bysort famid: gen total= _N 

*save twin data 
save "$temp\onlytwins.dta", replace 

use "$data\twinfinal10122019", clear
keep if sibling ==1 
drop total
bysort famid: gen total = _N
drop if total == 1

*same sex and opposite sex siblings
forvalues t = 2/7 {
bysort famid: egen sibsex`t' = total(female) if sibling==1 & total==`t'
  
}

*two sibling families
gen sibsex = 1 if sibsex2 == 0 // same sex Male
replace sibsex = 2 if sibsex2 == 2 // same sex female
replace sibsex = 3 if sibsex2 == 1 // opposite sex siblings

*three sibling families
replace sibsex = 1 if sibsex3 == 0 // same sex Male
replace sibsex = 2 if sibsex3 == 3 // same sex female
replace sibsex = 3 if inlist(sibsex3,1,2) // opposite sex siblings

*four sibling families
replace sibsex = 1 if sibsex4 == 0 // same sex Male
replace sibsex = 2 if sibsex4 == 4 // same sex female
replace sibsex = 3 if inlist(sibsex4,1,2,3) // opposite sex siblings

*five sibling families
replace sibsex = 1 if sibsex5 == 0 // same sex Male
replace sibsex = 2 if sibsex5 == 5 // same sex female
replace sibsex = 3 if inlist(sibsex5,1,2,3,4) // opposite sex siblings

*keep siblings 
keep if total==2 & twin==.

*labels  
label def sibsex 1 "SS Male sib" 2 "SS Female sib" 3 "OS siblings"
label val sibsex sibsex

*Append twins and siblings 
append using "$temp\onlytwins.dta",


*new variable for twins and siblings 
recode twinsex (0=1 "SS Male twin") (2=2 "SS Female twin") (1=3 "OS twin"), gen(sex_ts)
replace sex_ts = 4 if sibsex==1 // same sex Male
replace sex_ts = 5 if sibsex==2 // same sex female
replace sex_ts = 6 if sibsex==3 // opposite sex siblings
label def sex_ts 4 "SS Male sibs" , modify 
label def sex_ts 5 "SS Female sibs", modify 
label def sex_ts 6 "OS sibs", modify 


*income, education and ISEI for all 
gen linc = log(inc)
sum linc isei_max eduy
drop if isei_max==9999
*own isei std 
gen std_isei=.
gen std_linc=.
gen std_eduy=. 


forvalues t=0/1 { 
*ISEI
sum isei_max if female==`t' 
replace std_isei = (isei_max-r(mean))/r(sd) if female==`t' 
*income 
sum linc if female==`t' 
replace std_linc = (linc-r(mean))/r(sd) if female==`t' 
*edu at 30
sum eduy if female==`t' 
replace std_eduy = (eduy-r(mean))/r(sd) if female==`t' 

}


***************************************************************
															  *	
*** TWIN AND SIBLING CORRELATIONS *** 						  *			
															  * 
***************************************************************
forvalues t=1/6 { 
foreach x in y l h { 
gen `x'`t'=. 
	}
}
* income 
forvalues t = 1/6 {
eststo inc`t': mixed std_linc if sex_ts==`t' || famid:, reml var 
estat icc 
estadd scalar icc = r(icc2)
estadd scalar se = r(se2)
scalar icc`t' = r(icc2)
scalar h`t' = r(ci2)[1,2]
scalar l`t' = r(ci2)[1,1]
replace y1 = scalar(icc`t') if sex_ts==`t'
replace l1 =  scalar(l`t') if sex_ts==`t'
replace h1 =  scalar(h`t') if sex_ts==`t'
}


mixed std_isei if sex_ts==3 || famid:, mle var 
estat icc 


* ISEI 
forvalues t =1/6 {
eststo isei`t': mixed std_isei if sex_ts==`t' || famid:, reml var 
estat icc 
estadd scalar icc = r(icc2)
estadd scalar se = r(se2)
*these for figures  
scalar icc`t' = r(icc2)
scalar h`t' = r(ci2)[1,2]
scalar l`t' = r(ci2)[1,1]
replace y2 = scalar(icc`t') if sex_ts==`t'
replace l2 =  scalar(l`t') if sex_ts==`t'
replace h2 =  scalar(h`t') if sex_ts==`t'
}


* EDU 
forvalues t =1/6 {
eststo edu`t': mixed std_eduy if sex_ts==`t' || famid:, reml var 
estat icc 
estadd scalar icc = r(icc2)
estadd scalar se = r(se2)
*these for figures  
scalar icc`t' = r(icc2)
scalar h`t' = r(ci2)[1,2]
scalar l`t' = r(ci2)[1,1]
replace y3 = scalar(icc`t') if sex_ts==`t'
replace l3 =  scalar(l`t') if sex_ts==`t'
replace h3 =  scalar(h`t') if sex_ts==`t'

}


*Tables 
esttab inc1 inc2 inc3 inc4 inc5 inc6 using "$tables\ICC.rtf", s(icc se) drop(_cons) b(2) se(2) nopar nonum nodep compress mtitle("SS Male twin" "SS Female twin" "OS twin" "SS Male sibs" "SS Female sibs" "OS sibs") replace title("Income")
esttab isei1 isei2 isei3 isei4 isei5 isei6 using "$tables\ICC.rtf", s(icc se) drop(_cons) b(2) se(2) nopar nonum nodep compress mtitle("SS Male twin" "SS Female twin" "OS twin" "SS Male sibs" "SS Female sibs" "OS sibs") title("ISEI") append
esttab edu1 edu2 edu3 edu4 edu5 edu6 using "$tables\ICC.rtf", s(icc se) drop(_cons) b(2) se(2) nopar nonum nodep compress mtitle("SS Male twin" "SS Female twin" "OS twin" "SS Male sibs" "SS Female sibs" "OS sibs") title("Education") append




*Figures 
	set scheme plotplain

	twoway scatter y1 sex_ts,  mc(black) ms(o) || rcap l1 h1 sex_ts, lc(black) xlabel(1(1)6, valuelabels angle(45) labsize(small)) /// 
	xtitle("") yscale(r(0(0.1)0.6))  ylabel(0(0.1) 0.6) ytitle("ICC") legend(off) title("Intra Class Correlations for Income") name(inc,replace)

	twoway scatter y2 sex_ts,  mc(black) ms(o) || rcap l2 h2 sex_ts, lc(black) graphr(c(white)) xlabel(1(1)6, valuelabels angle(45) labsize(small)) /// 
	xtitle("") yscale(r(0(0.1)0.6))  ylabel(0(0.1) 0.6) ytitle("ICC", c(white)) legend(off) title("Intra Class Correlations for ISEI") name(isei,replace)

	twoway scatter y3 sex_ts,  mc(black) ms(o) || rcap l3 h3 sex_ts, lc(black) graphr(c(white)) xlabel(1(1)6, valuelabels angle(45) labsize(small)) /// 
	xtitle("") yscale(r(0(0.1)0.6))  ylabel(0(0.1) 0.6) ytitle("ICC") legend(off) title("Intra Class Correlations for Education") name(edu,replace)

gr combine inc isei edu , 
*gr export "$tables\ICC.pdf",replace
gr export "$tables\ICC.png",replace

