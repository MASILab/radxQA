function [Shewart,message] = ShewartEval(Val_sorted)

% 1: mean+/- 2D of previous measures
if abs(Val_sorted(end) - mean(Val_sorted(end-20:end-1))) > 2*std(Val_sorted(end-20:end-1))
    Shewart = 1; message = 'Measure out of 2*SD limits';
else 
    Shewart = 0; message = 'No Action Required: Within baseline measures';
end

if Shewart == 1
    % 2: exceeds 3D
    if abs(Val_sorted(end) - mean(Val_sorted(end-20:end-1))) > 3*std(Val_sorted(end-20:end-1))
        Shewart = 2; message = 'Out of 3*SD limits: May need evaluation';
    end
    % 3: two consecutive measures that exeed mean + 2STD
    if sum(abs(Val_sorted(end-1:end) - mean(Val_sorted(end-21:end-2))) > 2*std(Val_sorted(end-21:end-2)))==2
        Shewart = 3; message = '2 consecutive issues: Evaluate Instrument';
    end
    % 4: range of 4SD, difference between two consecutive measures exceeds 4SD
    if abs(Val_sorted(end)-Val_sorted(end-1))>4*std(Val_sorted(end-21:end-2))
        Shewart = 4; message = 'Big Change in system: Requires Evaluation';
    end
    % 5: 4 consectuive measures exeed same limit of 1SD or -1 SD
    if sum(abs(Val_sorted(end-3:end) - mean(Val_sorted(end-23:end-4))) > 1*std(Val_sorted(end-23:end-4)))==4
        Shewart = 5; message = '4 consecutive issues: Requires Evaluation';
    end
    if sum(Val_sorted(end-9:end) < mean(Val_sorted(end-29:end-10)))==0 || sum(Val_sorted(end-9:end) < mean(Val_sorted(end-29:end-10)))==10
        Shewart = 6; message = 'Consistent Drift: Instrument Evaluation Required';
    end
    
end



end
