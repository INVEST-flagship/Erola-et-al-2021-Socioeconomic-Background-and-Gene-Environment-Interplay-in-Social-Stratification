*** Socioeconomic Background and Gene–Environment Interplay in Social Stratification across the Early Life Course
*** Interaction ACE models by children's age 
*** date: 22.06.2020  
*** Hannu Lehti 

do "W:\Hannu L\twin correlations\01_do\variables.do"


*parental education in years 
forvalues x=1/5 {
foreach t in f m { 
gen `t'edu_years`x'=.

*gen `t'edu_years= . 
replace `t'edu_years`x' = 9 if `t'edu`x'==1 // basic edu 
replace `t'edu_years`x' = 11 if `t'edu`x'==2 // vocational 
replace `t'edu_years`x' = 12 if `t'edu`x'==3 // upper sec 
replace `t'edu_years`x' = 14 if `t'edu`x'==4 // opisto
replace `t'edu_years`x' = 15 if inlist(`t'edu`x',5,6) // polytechnics and bachelor uni
replace `t'edu_years`x' = 17 if inrange(`t'edu`x',7,9) // masters uni or higher
replace `t'edu_years`x' = 7 if `t'edu`x'==1 & `t'yb<1950 // kansakoulu
replace `t'edu_years`x' = 9 if `t'edu`x'==2 & `t'yb<1950 
replace `t'edu_years`x' = 8 if `t'edu`x'==1 & `t'yb>1949 & `t'yb<1965 // kansakoulu
replace `t'edu_years`x' = 10 if `t'edu`x'==2 & `t'yb>1949 & `t'yb<1965
	
	}
}

*highest parental edu in years
forvalues x=1/5 {
order fedu_years`x' fedu_years`x'
egen peduy`x' = rowmax(fedu_years`x' medu_years`x')
}

*pca
reshape long lfaminc pisei peduy, i(id) j(age)
pca lfaminc peduy pisei, com(1)
predict y // standardized latent factor for every year
keep lfaminc pisei peduy id age std_eduy std_isei std_inc famid zyg y
reshape wide y lfaminc pisei peduy, i(id) j(age) // reshape back long 

*ses categorical 5 periods 
xtile ses1 = y1, nq(5)
xtile ses2 = y2, nq(5)
xtile ses3 = y3, nq(5)
xtile ses4 = y4, nq(5)
xtile ses5 = y5, nq(5)

******************************** 

*create constant 
gen cons=1

*hetorogeneity at level 1 and level 2 for std_eduy 
forvalues t = 1/5 {
eq var2: var2 
eq var3: var3
eq het: y`t' cons 
eq slope: y`t' 
cons def 1 [pai2_1]var3 = [mem1_1]var2 

gllamm std_eduy, i(member pair) constr(1)  adapt eqs(var2 var3 slope) nrf(1 2) nocorr s(het) nip(7) ip(m) 
matrix list e(b)
nlcom  (E_cons: exp(2*[lns1]cons)) (E: 2*[lns1]y`t') (A_cons: ([pai2_1]var3)^2) (A: ([pai2_2]y`t')^2), post  // varianssit 
eststo e`t'
}

esttab e1 e2 e3 e4 e5
esttab e1 e2 e3 e4 e5 using "$tables\SES by age.rtf" , b(3) se(3) nostar replace title("Education") nopar compress nogaps nodep nonum /// 
mtitle("AGE: 0-5" "AGE: 6-10" "AGE: 11-15" "AGE: 16-20" "AGE: 21-25")

*hetorogeneity at level 1 and level 2 for ISEI by children's age 
forvalues t = 1/5 {
eq var2: var2 
eq var3: var3
eq het: y`t' cons 
eq slope: y`t' 
cons def 1 [pai2_1]var3 = [mem1_1]var2 

gllamm std_isei, i(member pair) constr(1)  adapt eqs(var2 var3 slope) nrf(1 2) nocorr s(het) nip(7) ip(m) 
nlcom  (E_cons: exp(2*[lns1]cons)) (E_slope: 2*[lns1]y`t') (A_cons: ([pai2_1]var3)^2) (A_slope: ([pai2_2]y`t')^2), post  // varianssit 
eststo p`t'
}

esttab p1 p2 p3 p4 p5 using "$tables\SES by age.rtf", se(3) nostar append title("ISEI") nopar compress nogaps nodep nonum /// 
mtitle("AGE: 0-5" "AGE: 6-10" "AGE: 11-15" "AGE: 16-20" "AGE: 21-25")

*hetorogeneity at level 1 and level 2 for income by children's age 
forvalues t = 1/5 {
eq var2: var2 
eq var3: var3
eq het: y`t' cons 
eq slope: y`t' 
cons def 1 [pai2_1]var3 = [mem1_1]var2 

gllamm std_inc, i(member pair) constr(1)  adapt eqs(var2 var3 slope) nrf(1 2) nocorr s(het) nip(7) ip(m) 
nlcom  (E_cons: exp(2*[lns1]cons)) (E_slope: 2*[lns1]y`t') (A_cons: ([pai2_1]var3)^2) (A_slope: ([pai2_2]y`t')^2), post  // varianssit 
eststo i`t'
}

esttab i1 i2 i3 i4 i5 using "$tables\SES by age.rtf", se(3) nostar append  title("Income") nopar compress nogaps nodep nonum ///  
mtitle("AGE: 0-5" "AGE: 6-10" "AGE: 11-15" "AGE: 16-20" "AGE: 21-25")



