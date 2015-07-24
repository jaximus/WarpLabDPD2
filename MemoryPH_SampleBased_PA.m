function PA_OutSample = MemoryPH_SampleBased_PA(PAinSample,PH_f1,PH_f3,PH_f5)                                         

persistent PreviousState1
persistent PreviousState3
persistent PreviousState5

Epsi_1 = PAinSample;
if isempty(PreviousState1)
    [PH_SampleBranch1 CurrentState1] = filter(PH_f1,1,Epsi_1);
else
    [PH_SampleBranch1 CurrentState1] = filter(PH_f1,1,Epsi_1,PreviousState1);
end
PreviousState1 = CurrentState1;

Epsi_3 = PAinSample.*abs(PAinSample).^2;
if isempty(PreviousState3)
    [PH_SampleBranch3 CurrentState3] = filter(PH_f3,1,Epsi_3);
else
    [PH_SampleBranch3 CurrentState3] = filter(PH_f3,1,Epsi_3,PreviousState3);
end
PreviousState3 = CurrentState3;

Epsi_5 = PAinSample.*abs(PAinSample).^4;
if isempty(PreviousState5)
    [PH_SampleBranch5 CurrentState5] = filter(PH_f5,1,Epsi_5);
else
    [PH_SampleBranch5 CurrentState5] = filter(PH_f5,1,Epsi_5,PreviousState5);
end
PreviousState5 = CurrentState5;

PA_OutSample = PH_SampleBranch1 + PH_SampleBranch3 + PH_SampleBranch5;

