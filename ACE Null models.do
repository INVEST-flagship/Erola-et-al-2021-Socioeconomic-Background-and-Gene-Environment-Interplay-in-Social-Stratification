*** Socioeconomic Background and Gene–Environment Interplay in Social Stratification across the Early Life Course 
*** Null (not controlled) ACE models for education, ISEI and income 
*** date: 22.10.2019 
*** Hannu Lehti 

do "W:\Hannu L\twin correlations\01_do\variables.do"

*** EMPTY ACE MODELS **** 

foreach t in std_eduy std_isei std_inc {
local x = sscor 
eststo `t': acelong `t' zyg famid id, vce(cluster famid) mzc(`x') mp
matrix b=r(b)
matrix se = r(se)
estadd scalar Avar = b[1,1], replace
estadd scalar Cvar = b[2,1], replace 
estadd scalar Evar = b[3,1], replace
estadd scalar Avse = se[1,1], replace
estadd scalar Cvse = se[2,1], replace 
estadd scalar Evse = se[3,1], replace
estadd scalar total = b[4,1], replace
estadd scalar total_se = se[4,1], replace
estadd scalar A = b[5,1], replace
estadd scalar C = b[6,1], replace
estadd scalar E = b[7,1], replace
estadd scalar Ase = se[5,1], replace
estadd scalar Cse = se[6,1], replace
estadd scalar Ese = se[7,1], replace
}

*Print table 
esttab std_eduy std_isei  using "$tables\ACE_edu_isei_inc_null.rtf", /// 
scalars("Avar A var" "Avse s.e." "Cvar C var" "Cvse s.e." "Evar E var" "Evse s.e." "total Total variance" "total_se s.e." /// 
"A A%" "Ase s.e." "C C%" "Cse s.e." "E E%" "Ese s.e." "N_clust N pairs") /// 
cells(Avar Avse Cvar Cvse Evar Evse total total_se A Ase C Cse E Ese N N_clust) sfmt(2) replace  /// 
mtitle("Edu at age 28" "ISEI" "Income") ///
collabels(none) nobase compress nogaps nopar nonum obslast 


*null model according to assortative mating 



*Calculate assortative mating according to parental education 
local x = sscor 
acelong std_eduy zyg famid id, vce(cluster famid) mzc(`x') mp
matrix ACE = r(b)
spearman medu5 fedu5
local pc =r(rho)
scalar dzc_e = .5 + .5 * (ACE[5,1]/100)*`pc'
dis dzc_e

*model for education 
local x = sscor 
local y = dzc_e
eststo astd_eduy: acelong std_eduy zyg famid id, vce(cluster famid) mzc(`x') dzc(`y') mp
matrix b=r(b)
matrix se = r(se)
estadd scalar Avar = b[1,1], replace
estadd scalar Cvar = b[2,1], replace 
estadd scalar Evar = b[3,1], replace
estadd scalar Avse = se[1,1], replace
estadd scalar Cvse = se[2,1], replace 
estadd scalar Evse = se[3,1], replace
estadd scalar total = b[4,1], replace
estadd scalar total_se = se[4,1], replace
estadd scalar A = b[5,1], replace
estadd scalar C = b[6,1], replace
estadd scalar E = b[7,1], replace
estadd scalar Ase = se[5,1], replace
estadd scalar Cse = se[6,1], replace
estadd scalar Ese = se[7,1], replace

*Calculate assortative mating according to parental ISEI 
local x = sscor 
acelong std_isei zyg famid id, vce(cluster famid) mzc(`x') mp
matrix ACE = r(b)
corr misei5 fisei5
local pc =r(rho)
scalar dzc_i = .5 + .5 * (ACE[5,1]/100)*`pc'
dis dzc_i


*model for ISEI 
local x = sscor 
local y = dzc_i
eststo astd_isei: acelong std_isei zyg famid id, vce(cluster famid) mzc(`x') dzc(`y') mp
matrix b=r(b)
matrix se = r(se)
estadd scalar Avar = b[1,1], replace
estadd scalar Cvar = b[2,1], replace 
estadd scalar Evar = b[3,1], replace
estadd scalar Avse = se[1,1], replace
estadd scalar Cvse = se[2,1], replace 
estadd scalar Evse = se[3,1], replace
estadd scalar total = b[4,1], replace
estadd scalar total_se = se[4,1], replace
estadd scalar A = b[5,1], replace
estadd scalar C = b[6,1], replace
estadd scalar E = b[7,1], replace
estadd scalar Ase = se[5,1], replace
estadd scalar Cse = se[6,1], replace
estadd scalar Ese = se[7,1], replace

*For income the model does not alter because no C (same as AE model)
local x = sscor 
local y = dzc_i
eststo astd_linc: aelong std_linc zyg famid id, vce(cluster famid) mzc(`x') mp
matrix b=r(b)
matrix se = r(se)
estadd scalar Avar = b[1,1], replace
*estadd scalar Cvar = b[2,1], replace 
estadd scalar Evar = b[2,1], replace
estadd scalar Avse = se[1,1], replace
*estadd scalar Cvse = se[2,1], replace 
estadd scalar Evse = se[2,1], replace
estadd scalar total = b[3,1], replace
estadd scalar total_se = se[3,1], replace
estadd scalar A = b[4,1], replace
*estadd scalar C = b[6,1], replace
estadd scalar E = b[5,1], replace
estadd scalar Ase = se[4,1], replace
*estadd scalar Cse = se[6,1], replace
estadd scalar Ese = se[5,1], replace


*print table for assortative mating models 
esttab astd_eduy astd_isei astd_linc using "$tables\ACE_edu_isei_inc_null_mating.rtf", /// 
scalars("Avar A var" "Avse s.e." "Cvar C var" "Cvse s.e." "Evar E var" "Evse s.e." "total Total variance" "total_se s.e." /// 
"A A%" "Ase s.e." "C C%" "Cse s.e." "E E%" "Ese s.e." "N_clust N pairs") /// 
cells(Avar Avse Cvar Cvse Evar Evse total total_se A Ase C Cse E Ese N N_clust) sfmt(2)  /// 
mtitle("Edu at age 28" "ISEI" "Income") ///
collabels(none) nobase compress nogaps nopar nonum obslast replace title("Assortative mating adjusted")

