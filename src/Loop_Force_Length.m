lcerel_length=[0.5:0.15:1.40];

for i=1:7

    [AAA, BBB] = ce_fl_simple(lcerel_length(i),parms);
    f_iso_rel(i)=AAA;

end

plot(lcerel_length,f_iso_rel)