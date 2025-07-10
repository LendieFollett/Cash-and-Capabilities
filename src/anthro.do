* PROJECT: A new kind of impact evaluation
* BY: Lendie Follett, Heath Henderson, and Smriti Tiwari
* LAST UPDATE: 6/26/2025
* LOCATION: /Users/hendersonhl/Documents/Articles/Impact-Evaluation/Cleaning/CT-OVC
* PURPOSE: To examine negative impacts on child anthropometrics

* Set up
clear all
cd "/Users/hendersonhl/Documents/Articles/Impact-Evaluation/Cleaning/CT-OVC/"
import delimited "ctovc_final"
gen interact = treated*year

* Results using OLS
reg stunting treated year interact caregiver-district6, cluster(location)
reg wasting treated year interact caregiver-district6, cluster(location)
reg underweight treated year interact caregiver-district6, cluster(location)

* Heterogeneous effects by age
sum caregiver if sample_nutrition==1, detail   // Median age is 48
foreach var of varlist stunting-underweight {
    reg `var' treated year interact caregiver-district6 ///
	if caregiver<48, cluster(location)
    reg `var' treated year interact caregiver-district6 ///
	if caregiver>=48, cluster(location)		
}

*** Outliers

* Winsorize anthropometric variables
winsor stunting, p(0.05) gen(stunting2)   
winsor wasting, p(0.05) gen(wasting2)
winsor underweight, p(0.05) gen(underweight2)

* Regressions with winsorized variables
reg stunting2 treated year interact caregiver-district6, cluster(location)
reg wasting2 treated year interact caregiver-district6, cluster(location)
reg underweight2 treated year interact caregiver-district6, cluster(location)

* Dichotomize anthropometric variables
gen stunting3 = .
replace stunting3=0 if stunting>-2 & stunting!=.
replace stunting3=1 if stunting<=-2 & stunting!=.
gen wasting3 = .
replace wasting3=0 if wasting>-2 & wasting!=.
replace wasting3=1 if wasting<=-2 & wasting!=.
gen underweight3 = .
replace underweight3=0 if underweight>-2 & underweight!=.
replace underweight3=1 if underweight<=-2 & underweight!=.

* Regressions with dichotomized variables
reg stunting3 treated year interact caregiver-district6, cluster(location)
reg wasting3 treated year interact caregiver-district6, cluster(location)
reg underweight3 treated year interact caregiver-district6, cluster(location)

*** Price effects

* Open wave 3 data
clear all
use community_w3.dta

* Generate treatment and wave variable
gen treated = .
replace treated = 1 if cw3_h3==1
replace treated = 0 if cw3_h3==2
gen year=1

* Keep, order, and save
renpfix cw3_e_
keep treated year price1-price10
order treated year price1-price10
tempfile cw3
save `cw3'

* Open wave 1 data
* Note: First get treatment information from HH survey
clear all
use hh_w1.dta
keep commidw1 hhw1_househol
gen treated = 0
replace treated=1 if hhw1_househol=="A"
collapse treated, by(commidw1)
replace treated=1 if treated>0
tempfile hhw1
save `hhw1'

* Open community dataset
clear all
use community_w1.dta

* Keep price data 
keep commidw1 cw1_price1-cw1_price10
renpfix cw1_

* Merge in treatment information
merge 1:1 commidw1 using `hhw1'
drop _merge

* Miscellaneous cleaning 
gen year=0
drop commidw1
order treated year price*

* Append wave 3 data
append using `cw3'

* Remove missing values
foreach var of varlist price1-price10 {
	replace `var' = . if `var'==9998 | `var'==9999 | `var'==998
}

* Difference-in-differences regressions
gen interact = treated*year
foreach var of varlist price1-price10 {
	reg `var' treated year interact	
}

*** Household size effects

* Open wave 3 data
clear
use hh_w3.dta
gen treated = 0
replace treated = 1 if hhw3_j9==1 | hhw3_j9==2
gen year = 1

* Generate HH size variable and save
collapse treated year (count) indivcode, by(hhcode)
rename indivcode size
tempfile hhw3
save `hhw3'

* Open wave 1 data
clear
use hh_w1.dta
gen treated = 0
replace treated = 1 if hhw1_househol=="A"
gen year = 0

* Generate HH size variable
collapse treated year (count) indivcode, by(hhcode)
rename indivcode size

* Append wave 1 dataset and run regressions
* Note: SEs are not clustered here, but it will not change conclusion.
append using `hhw3'
gen interact = treated*year
reg size treated year interact




























