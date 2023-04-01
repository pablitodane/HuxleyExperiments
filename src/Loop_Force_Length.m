%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isometric force (you can comment and uncomment it)

lce_length = [0.009 : 0.0002 : 0.02];

lcerel=lce_length/parms.lceopt;

for i=1:length(lcerel)

    [AAA, BBB] = ce_fl_simple(lcerel(i),parms);
    f_iso_rel(i)=AAA;

end

f_real = f_iso_rel*parms.Fmax;

plot(lce_length,f_real)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force F_pe

% lpe_length = [0: 0.0005: 0.0186];
% 
% for idx=1:length(lpe_length)
% 
%     [fse_loop, fpe_loop, kse_loop, kpe_loop] = CEEC_simple(nan,lpe_length(idx),parms);
%     fpe_loop_rel(idx) = fpe_loop;
% 
% end
% 
% plot(lpe_length,fpe_loop_rel)
