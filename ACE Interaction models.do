*** Socioeconomic Background and Gene–Environment Interplay in Social Stratification across the Early Life Course
*** Interaction ACE models
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

******************************** 


*GLLAMM models for interaction


* ACE interaction model by parental SES and income  
foreach t in std_eduy std_isei std_inc  {  
foreach x in  s f  {
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

* ACE interaction model by parental ISEI 
foreach t in  std_eduy std_isei std_inc   {  
foreach x in  z  {
eq var2: var2 
eq var3: var3
eq het:   `x'5 `x'1 `x'2 `x'3 `x'4   
eq inter: `x'5 `x'1 `x'2 `x'3 `x'4 
cons def 1 [pai2_1]var3 = [mem1_1]var2 
*model 
gllamm `t', i(member pair) constr(1) adapt eqs(var2 var3 inter) nrf(1 2) s(het) nip(7) ip(m) nocorr 
eststo `t'_`x' // save estimates 
nlcom  (`x'1: exp(2*[lns1]`x'1))  (`x'2: exp(2*[lns1]`x'2)) (`x'3: exp(2*[lns1]`x'3)) (`x'4: exp(2*[lns1]`x'4)) (`x'5: exp(2*[lns1]`x'5)), post  
eststo `t'_`x'_e
estimates restore `t'_`x'
nlcom  (`x'1: ([pai2_1]var3^2 + ([pai2_2]`x'5)^2 * ([pai2_2]`x'1)^2)) /// 
	   (`x'2: ([pai2_1]var3^2 + ([pai2_2]`x'5)^2 * ([pai2_2]`x'2)^2)) ///
	   (`x'3: ([pai2_1]var3^2 + ([pai2_2]`x'5)^2 * ([pai2_2]`x'3)^2))  ///
	   (`x'4: ([pai2_1]var3^2 + ([pai2_2]`x'5)^2 * ([pai2_2]`x'4)^2)) /// 
	   (`x'5: (([pai2_1]var3^2) + ([pai2_2]`x'5^2))), post  
eststo `t'_`x'_a
	}
} 

*ACE interaction model by 4 categorical parental education 
foreach t in std_eduy std_isei std_inc {  
eq var2: var2 
eq var3: var3
eq het:  p1 p2 p3 p4 
eq inter: p1 p2 p3 p4 
cons def 1 [pai2_1]var3 = [mem1_1]var2 
*model 
gllamm `t', i(member pair) constr(1) adapt eqs(var2 var3 inter) nrf(1 2) ip(m) s(het) nip(7)  nocorr
eststo `t'  // save estimates 
nlcom  (p1: exp(2*[lns1]p1)) (p2: exp(2*[lns1]p2)) (p3: exp(2*[lns1]p3)) (p4: exp(2*[lns1]p4))  , post  
eststo `t'_p_e
estimates restore `t'
nlcom  (p1: (([pai2_1]var3^2) + ([pai2_2]p1^2))) (p2: ([pai2_1]var3^2 +  ([pai2_2]p1)^2 * ([pai2_2l]p2)^2)) (p3: ([pai2_1]var3^2 +  ([pai2_2]p1)^2 * ([pai2_2l]p3)^2)) (p4: ([pai2_1]var3^2 +  ([pai2_2]p1)^2 * ([pai2_2l]p4)^2)), post 
eststo `t'_p_a 

}

***********************
*Tables               * 
***********************
esttab std_eduy_s_a std_eduy_s_e std_isei_s_a std_isei_s_e std_inc_s_a std_inc_s_e  /// 
using "$tables\Interactions all.rtf" , mtitle("Edu A" "Edu E" "ISEI A" "ISEI E" "Income A" "Income E") title("Parental SES") modelwidth(8) varwidth(11) nonum nostar label se(3) b(3) nopar replace 

esttab std_eduy_p_a std_eduy_p_e std_isei_p_a std_isei_p_e std_inc_p_a std_inc_p_e  /// 
using "$tables\Interactions all.rtf" , mtitle("Edu A" "Edu E" "ISEI A" "ISEI E" "Icome A" "Income E") title("Parental Education") modelwidth(8) varwidth(11) nonum nostar label se(3) b(3) nopar append 

esttab std_eduy_z_a std_eduy_z_e std_isei_z_a std_isei_z_e std_inc_z_a std_inc_z_e  /// 
using "$tables\Interactions all.rtf" , mtitle("Edu A" "Edu E" "ISEI A" "ISEI E" "Income A" "Income E") title("Parental ISEI") modelwidth(8) varwidth(11) nonum nostar label se(3) b(3) nopar append 

esttab std_eduy_f_a std_eduy_f_e std_isei_f_a std_isei_f_e std_inc_f_a std_inc_f_e  /// 
using "$tables\Interactions all.rtf" , mtitle("Edu A" "Edu E" "ISEI A" "ISEI E" "Icome A" "Income E") title("Parental Income") modelwidth(8) varwidth(11) nonum nostar label se(3) b(3) nopar append 

*************************************************************
*Figures by all the parental variables 			            *
*************************************************************

*Parental SES 		
coefplot (std_eduy_s_a, label("A") offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_eduy_s_e, label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), /// 
title("Education",  box bexpand fc(gs13)) plotr(lc(black)) name(s1,replace) vertical  recast(connected)  cirecast(rcap) ytitle("Variance", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) xtitle("Parental SES", size(small))

coefplot (std_isei_s_a, label("A") offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_isei_s_e, label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), /// 
title("ISEI", box bexpand fc(gs13)) plotr(lc(black)) name(s2,replace) vertical  recast(connected)  cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) xtitle("Parental SES", size(small))

coefplot (std_inc_s_a, label("A") offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_inc_s_e, label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))) , ///
title("Income", box bexpand fc(gs13)) plotr(lc(black)) name(s3,replace) vertical recast(connected) cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) xtitle("Parental SES", size(small))

*Parental income 
coefplot (std_eduy_f_a, label("A") offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_eduy_f_e,  label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), ///
title("Education",  box bexpand fc(gs13)) plotr(lc(black)) name(f1,replace) vertical  recast(connected)  cirecast(rcap) ytitle("Variance", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) xtitle("Family income", size(small))

coefplot (std_isei_f_a, label("A")  offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_isei_f_e, label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), /// 
title("ISEI",  box bexpand fc(gs13)) plotr(lc(black)) name(f2,replace) vertical  recast(connected)  cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) ///
xtitle("Family income", size(small))

coefplot (std_inc_f_a, label("A")  offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_inc_f_e, label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))) , /// 
title("Income",  box bexpand fc(gs13)) plotr(lc(black)) name(f3,replace) vertical  recast(connected)  cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) ///
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) xtitle("Family income", size(small))

*Parental ISEI 
coefplot (std_eduy_z_a, label("A") offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_eduy_z_e,  label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), /// 
title("Education", box bexpand fc(gs13)) plotr(lc(black)) name(z1,replace) vertical  recast(connected)  cirecast(rcap) ytitle("Variance", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall))  yscale(r(0(0.1)1.4)) xtitle("Parental ISEI", size(small))

coefplot (std_isei_z_a, label("A")  offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_isei_z_e,  label("E")  offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), /// 
title("ISEI", box bexpand fc(gs13)) plotr(lc(black)) name(z2,replace) vertical  recast(connected)  cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) ///
ylabel(0(0.1)1.4, labsize(vsmall))  yscale(r(0(0.1)1.4)) xtitle("Parental ISEI", size(small))

coefplot (std_inc_z_a, label("A")  offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_inc_z_e,  label("E")  offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), /// 
title("Income", box bexpand fc(gs13)) plotr(lc(black)) name(z3,replace) vertical  recast(connected) cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall))  yscale(r(0(0.1)1.4)) xtitle("Parental ISEI", size(small))

*Parental edu 
coefplot (std_eduy_p_a, label("A") offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_eduy_p_e, label("E") offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), ///
title("Education", box bexpand fc(gs13)) plotr(lc(black)) name(p1,replace) vertical  recast(connected)  cirecast(rcap) ytitle("Variance", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) ///
ylabel(0(0.1)1.2, labsize(vsmall)) yscale(r(0(0.1)1.4)) xtitle("Parental Education", size(small)) xlabel(,labsize(small))

coefplot (std_isei_p_a, label("A")  offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_isei_p_e, label("E")  offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))), ///
title("ISEI", box bexpand fc(gs13)) plotr(lc(black)) name(p2,replace) vertical  recast(connected)  cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale( r(0(0.1)1.4)) xtitle("Parental Education", size(small))  xlabel(,labsize(small))

coefplot (std_inc_p_a, label("A")  offset(0) lp(solid) mc(red) ms(o) lcolor(red) ciopts(lp(solid) color(red))) (std_inc_p_e, label("E")  offset(0.03) lp(solid) mc(blue) ms(o) lcolor(blue) ciopts(lp(solid) color(blue))) , /// 
title("Income",  box bexpand fc(gs13)) plotr(lc(black) margin(top)) name(p3,replace) vertical  recast(connected)  cirecast(rcap) ytitle("", size(vsmall) margin(zero)) legend(row(1) size(small) pos(6)) /// 
ylabel(0(0.1)1.4, labsize(vsmall)) yscale(r(0(0.1)1.4)) xtitle("Parental Education", size(small)) xlabel(,labsize(small))

*Print interaction figures by parental variable 
grc1leg s1 s2 s3, row(1) name(ses,replace) // parental SES 
gr export "$tables\int_par_ses.pdf", replace
grc1leg f1 f2 f3, row(1) name(inc,replace) // parental income 
gr export "$tables\int_par_inc.pdf", replace
grc1leg z1 z2 z3, row(1) name(sei, replace) // parental ISEI 
gr export "$tables\int_par_isei.pdf", replace		
grc1leg p1 p2 p3, row(1) name(edu, replace) // parental edu  
gr export "$tables\int_par_edu.pdf", replace	

*Print all the figures in the same graph   	
grc1leg p1 p2 p3 z1 z2 z3 f1 f2 f3 s1 s2 s3, row(4) 
gr display, xsize(10) ysize(13) 
gr export "$tables\int_all.pdf", replace	

