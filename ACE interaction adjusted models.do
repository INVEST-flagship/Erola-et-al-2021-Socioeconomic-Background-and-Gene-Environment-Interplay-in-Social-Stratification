*** Socioeconomic Background and Gene–Environment Interplay in Social Stratification across the Early Life Course
*** Adjusted interaction ACE models for children's income and ISEI
*** date: 22.06.2020  
*** Hannu Lehti 

do "W:\Hannu L\twin correlations\01_do\variables.do"

* Variables that is needed for PCA *  

*parental education in years 
forvalues x=1/5 {
foreach t in f m { 
gen `t'edu_years`x'=. 
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

*take highest parental education in years
forvalues x=1/5 {
order fedu_years`x' fedu_years`x'
egen peduy`x' = rowmax(fedu_years`x' medu_years`x')
}

*principal componen analysis (PCA)
reshape long lfaminc pisei peduy, i(id) j(age)
pca lfaminc peduy pisei, com(1)
predict y // standard latent factor every year 
keep lfaminc pisei peduy id age std_eduy std_isei std_inc famid zyg y pedu
reshape wide y lfaminc pisei peduy , i(id) j(age) // reshape back long 

*parental edu 4 categories 
recode pedu3 (1=1 "Basic") (2/3 = 2 "Secondary") (4=3 "Post-secondary") (5/6 = 4 "Tertiary"),gen(pedu)
label var pedu "Parental education"

recode pedu3 (1 = 1 "Basic") (2/3=2 "Secondary") (4/6 = 3 "Tertiary"), gen(paredu)
label var paredu "Parental education"


* take quintiles of parental variables 
xtile ses = y3, nq(5)
xtile finc = lfaminc3, nq(5) 
xtile zpisei = pisei3, nq(5)
xtile zpisei10 = pisei3, nq(10)


*Parental SES quintiles at 5 periods 
xtile ses1 = y1, nq(5)
xtile ses2 = y2, nq(5)
xtile ses3 = y3, nq(5)
xtile ses4 = y4, nq(5)
xtile ses5 = y5, nq(5)

forvalues t = 1/5 { 
ta ses`t',gen(se`t')	
	
}

*code dummy variables for the interaction ACE models 
ta ses, gen(s)
ta finc, gen(f)
ta zpisei, gen(z)
ta pedu, gen(p)
ta paredu, gen(par)

label var p1 "Basic"
label var p2 "Sec."
label var p3 "Post-sec."
label var p4 "Tertiary"


label var par1 "Basic"
label var par2 "Secondary"
label var par3 "Tertiary"

foreach t in s f z {
	forvalues y =1/5 {
	
	label var `t'`y' "`y' Q."
	
	}
}
*family and individual level variables for ACE interaction models 
gen mid = famid if zyg==1 
replace mid = id if zyg==2 
gen pair = famid
gen member = mid

*Calculate zygosity weights for ACE interaction models 
ta zyg, gen(m)
gen var3 = m1+sqrt(0.5)*(1-m1) 
gen var2 = sqrt(sscor-0.5)*(1-m1) 


************************************************************************
*Adjusted AE interaction models for children's income and ISEI´********* 
************************************************************************

*baseline AE interaction models for ISEI and INC (without control variables) 
foreach t in  std_isei std_inc  {  
foreach x in  s  {
eq var2: var2 
eq var3: var3
eq het:   `x'1 `x'2 `x'3 `x'4 `x'5  
eq inter: `x'1 `x'2 `x'3 `x'4 `x'5 
cons def 1 [pai2_1]var3 = [mem1_1]var2 
*model 
gllamm `t', i(member pair) constr(1) adapt eqs(var2 var3 inter) nrf(1 2) s(het) nip(7) ip(m) nocorr 
eststo `t'_`x' // save estimates 
nlcom  (`x'1: exp(2*[lns1]`x'1))  (`x'2: exp(2*[lns1]`x'2)) (`x'3: exp(2*[lns1]`x'3)) (`x'4: exp(2*[lns1]`x'4)) (`x'5: exp(2*[lns1]`x'5)), post  
eststo `t'_`x'_e

estimates restore `t'_`x'
nlcom  (`x'1: (([pai2_1]var3^2) + ([pai2_2]`x'1^2))) ///
	   (`x'2: ([pai2_1]var3^2 + ([pai2_2]`x'1)^2 * ([pai2_2]`x'2)^2)) ///
	   (`x'3: ([pai2_1]var3^2 + ([pai2_2]`x'1)^2 * ([pai2_2]`x'3)^2))  ///
	   (`x'4: ([pai2_1]var3^2 + ([pai2_2]`x'1)^2 * ([pai2_2]`x'4)^2)) /// 
	   (`x'5: ([pai2_1]var3^2 + ([pai2_2]`x'1)^2 * ([pai2_2]`x'5)^2)), post  
eststo `t'_`x'_a
	}
} 


*AE interaction model, for ISEI children's education controlled for  
eq var2: var2 
eq var3: var3
eq inter: s1 s2 s3 s4 s5
eq het: s1 s2 s3 s4 s5 
cons def 1 [pai2_1]var3 = [mem1_1]var2 
*model 
gllamm std_isei std_eduy, i(member pair) constr(1) adapt eqs(var2 var3 inter) s(het) nrf(1 2) ip(m) nocorr
eststo s 
nlcom  (s1: exp(2*[lns1]s1)) (s2: exp(2*[lns1]s2)) (s3: exp(2*[lns1]s3)) (s4: exp(2*[lns1]s4)) (s5: exp(2*[lns1]s5)), post  
eststo ses_e 
est restore s
nlcom  	(t1: ([pai2_1]var3^2) + ([pai2_2]s1^2)) /// 
		(t2: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s2)^2)) /// 
	    (t3: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s3)^2)) /// 
	    (t4: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s4)^2)) /// 
		(t5: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s5)^2)), post 
eststo ses_a 

*print table  
esttab std_isei_s_a std_isei_s_e ses_a ses_e using  "$tables\isei_edu_controlled_categories.rtf",  /// 
ci(2) b(3) replace nopar nostar compress nogaps nodep nonum title("ISEI") mtitle("Baseline A" "Baseline E" "Own edu controlled A"  "Own edu controlled E") 

*AE interaction for income, children's education controlled for 
eq var2: var2 
eq var3: var3
eq inter: s1 s2 s3 s4 s5
eq het: s1 s2 s3 s4 s5 
cons def 1 [pai2_1]var3 = [mem1_1]var2 
  
gllamm std_inc std_eduy, i(member pair) constr(1) adapt eqs(var2 var3 inter) s(het) nrf(1 2) ip(m) nocorr // oma koulutus kontrolloitu 
eststo inc 
nlcom  (s1: exp(2*[lns1]s1)) (s2: exp(2*[lns1]s2)) (s3: exp(2*[lns1]s3)) (s4: exp(2*[lns1]s4)) (s5: exp(2*[lns1]s5)), post  
eststo inc_e1
est restore inc
nlcom  	(s1: ([pai2_1]var3^2) + ([pai2_2]s1^2)) /// 
		(s2: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s2)^2)) /// 
	    (s3: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s3)^2)) /// 
	    (s4: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s4)^2)) /// 
		(s5: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s5)^2)), post
eststo inc_a1 

*AE interaction model for income, children's education and ISEI controlled for
eq var2: var2 
eq var3: var3
eq inter: s1 s2 s3 s4 s5
eq het: s1 s2 s3 s4 s5 
cons def 1 [pai2_1]var3 = [mem1_1]var2 
gllamm std_inc std_eduy std_isei, i(member pair) constr(1) adapt eqs(var2 var3 inter) s(het) nrf(1 2) ip(m) nocorr
eststo inc2 
nlcom  (s1: exp(2*[lns1]s1)) (s2: exp(2*[lns1]s2)) (s3: exp(2*[lns1]s3)) (s4: exp(2*[lns1]s4)) (s5: exp(2*[lns1]s5)), post  
eststo inc_e2 
est restore inc2
nlcom  	(s1_1: ([pai2_1]var3^2) + ([pai2_2]s1^2)) /// 
		(s2_1: 100*([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s2)^2)) /// 
	    (s3_1: 100*([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s3)^2)) /// 
	    (s4_1: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s4)^2)) /// 
		(s5_1: ([pai2_1]var3^2 +  ([pai2_2]s1)^2 * ([pai2_2l]s5)^2)), post
nlcom (s1: _b[s1_1]) (s2: _b[s2_1]/100) (s3: _b[s3_1]/100) (s4: _b[s4_1]) (s5: _b[s5_1]) , post 		
eststo	inc_a2	 

*Print table 
esttab  std_inc_s_a std_inc_s_e inc_a1 inc_e1 inc_a2 inc_e2 using  "$tables\inc_edu_controlled_categories.rtf",  ci(3) b(3) replace nopar nostar compress nogaps nodep nonum title("Income") /// 
mtitle("Baseline A" "Baseline E" "Own edu controlled A"  "Own edu controlled E" "Own edu and ISEI controlled A"  "Own edu and ISEI controlled E") 

