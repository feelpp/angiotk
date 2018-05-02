
unit=10;//1e-1;
lc=0.05;//0.001;//0.0001;//0.1;
cptPt=1;
cptLine=1;
r=0.05;//*unit;

For i In {0:14*Pi:0.1}
  Point(cptPt) = { r*Cos(i)*unit,r*Sin(i)*unit,0.01*i*unit,lc};
  cptPt=cptPt+1;
  If ( cptPt>2 )
    Line(cptLine) = { cptPt-2,cptPt-1 };
    cptLine=cptLine+1;
  EndIf
EndFor

