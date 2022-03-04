** QQ plot of EPO meta-analysis **
** Charli Harlow **

* make sure p-value column is called “p"
rename pvalue p
sort p
gen float log_observed_p = -log10(p)

* determine number of rows (e.g. 13259219)
range expected_p 1/25090314 1
generate float log_expected_p = -log10(expected_p)

range up_n 1 25090314 25090314
range down_n 25090314 1 25090314
generate float upper_p = invibeta(up_n,down_n,0.025)
generate float lower_p = invibeta(up_n,down_n,0.975)
generate double log_upper_95_p = -log10(upper_p)
generate double log_lower_95_p = -log10(lower_p)

twoway (rarea log_upper_95_p log_lower_95_p log_expected_p, fcolor(gs13) lcolor(gs13)) (line log_expected_p log_expected_p, lcolor(gs6) lwidth(vthin)) (scatter log_observed_p log_expected_p, mcolor(blue) msize(small) msymbol(circle)) if p < 0.01,  ytitle (Observed -log10 P) xtitle(Expected -log10 P)  ylabel(, nogrid) ymtick(, nogrid) xscale(line) xlabel(, nogrid) xmtick(, nogrid) legend(off) graphregion(fcolor(white) lcolor(white) lwidth(thin) lpattern(solid) ifcolor(white) ilcolor(white) ilwidth(thin) ilpattern(solid)) plotregion(fcolor(white) lcolor(white) lwidth(thin) lpattern(solid) ifcolor(white) ilcolor(white) ilwidth(thin) ilpattern(solid))
