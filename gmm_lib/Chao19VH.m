function [VH,SaH,SaV]=Chao19VH(T,M,Rrup,Ztor,Vs30,Z10,FRO,FSS,FNO,Finter,Fintra,Fas,Fma,FVS30)

SaH=Chao19H(T,M,Rrup,Ztor,Vs30,Z10,FRO,FSS,FNO,Finter,Fintra,Fas,Fma,FVS30);
SaV=Chao19V(T,M,Rrup,Ztor,Vs30,Z10,FRO,FSS,FNO,Finter,Fintra,Fas,Fma,FVS30);
VH=SaV./SaH;
