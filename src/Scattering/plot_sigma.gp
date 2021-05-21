reset
set term pngcairo
set termoption enhanced
set datafile separator '\t'
set output 'total_CrossSection.png'
set xlabel 'E [meV]'
set ylabel '{/Symbol S}_{tot} [A^2]' offset 1,0 rotate by 90

unset key
set grid 

plot 'sigma_tot.csv' using 1:2 with lines lt 3 title "Total cross-section"
system('eog total_CrossSection.png')