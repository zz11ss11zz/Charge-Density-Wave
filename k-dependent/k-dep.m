(* ::Package:: *)

SetSystemOptions["ParallelOptions" -> "MathLinkTimeout" -> 30.]

LaunchKernels[16];  (*parallel 16 cores*)

inputkxFileName="../inputs/kx.txt";
inputkyFileName="../inputs/ky.txt";
inputEnergyFileName="../inputs/energy.txt";
inputqFileName="../inputs/111.dat";
inputw90FileName="../inputs/wannier90_geninterp.dat";

initialOutput={};

T=0.005;  (*electronic temperature*)
fermiDirac[x_]:=1/(Exp[x/T]+1);
dFD[x_]:=-Exp[x/T]/((1+Exp[x/T])^2* T);

kx=ReadList[OpenRead[inputkxFileName],Number];
ky=ReadList[OpenRead[inputkyFileName],Number];
band=ReadList[OpenRead[inputEnergyFileName],Number];
gData=Import[inputqFileName,"Table"][[All,1]];
w90Data=Import[inputw90FileName,"Table"];


w90vx=w90Data[[All,-3]];
w90vy=w90Data[[All,-2]];

nPts=Length[kx];
gData=Partition[gData,nPts];

data={};
kvector={};
For[i=1,i<=Length[ky],i++,
AppendTo[data,{kx[[i]],ky[[i]],band[[i]]}];
AppendTo[kvector,{kx[[i]],ky[[i]]}]
]
bz=ConvexHullMesh @ kvector;
periodx=Max[kx]-Min[kx];
periody=Max[ky]-Min[ky];

upperLeftky=Min[ky];
lowerRightky=Max[ky];
For[i=1,i<=Length[ky],i++,
thisvector={kx[[i]],ky[[i]]};
If[thisvector[[1]]==Min[kx]&&thisvector[[2]]>upperLeftky,
upperLeftky=thisvector[[2]];
upperLeftkvector=thisvector
];
If[thisvector[[1]]==Max[kx]&&thisvector[[2]]<lowerRightky,
lowerRightky=thisvector[[2]];
lowerRightkvector=thisvector
];
]

basisvector1={Max[kx],Max[ky]}-{Min[kx],upperLeftkvector[[2]]};
basisvector2={0,upperLeftkvector[[2]]-Min[ky]};

indexqmax=Length[kx];
indexkmax=Length[kx];

ParallelEvaluate[file=OpenWrite["outputCore"<>ToString[$KernelID]<>".txt"]];
ParallelDo[
thiskx=kx[[indexk]];
thisky=ky[[indexk]];
thiskvector={thiskx,thisky};
sus=0;
susIm=0;
g2SumEPW=0;
g2SumW90=0;
d2SumEPW=0;
d2SumW90=0;

Do[
thisqx=kx[[indexq]];
thisqy=ky[[indexq]];
thisqvector={thisqx,thisqy};
thiskplusqx=thiskx+thisqx;
thiskplusqy=thisky+thisqy;
thiskplusqvector={thiskplusqx,thiskplusqy};
If[RegionMember[bz,thiskplusqvector]==False,While[thiskplusqvector[[1]]>Max[kx],thiskplusqvector=thiskplusqvector-basisvector1];
While[thiskplusqvector[[1]]<Min[kx],thiskplusqvector=thiskplusqvector+basisvector1];
];
If[thiskplusqvector[[2]]>Max[ky],While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector-basisvector2]];
If[thiskplusqvector[[2]]<Min[ky],While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector+basisvector2]];
If[thiskplusqvector[[2]]>upperLeftky,While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector-basisvector2]];
If[thiskplusqvector[[2]]<lowerRightky,While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector+basisvector2]];

nearestkplusqvector=Nearest[kvector,thiskplusqvector][[1]];
indexNearestkplusqvector=Position[kvector,nearestkplusqvector][[1]][[1]];
thisgEPW=gData[[indexq,indexk]];
thisgxW90=w90vx[[indexNearestkplusqvector]]-w90vx[[indexk]];
thisgyW90=w90vy[[indexNearestkplusqvector]]-w90vy[[indexk]];
thisg2W90=thisgxW90^2+thisgyW90^2;
\[Epsilon]k=band[[indexk]];
\[Epsilon]kplusq=band[[indexNearestkplusqvector]];
f1=fermiDirac[\[Epsilon]k];
f2=fermiDirac[\[Epsilon]kplusq];
If[\[Epsilon]k==\[Epsilon]kplusq,thisSusRe=0;];
If[\[Epsilon]k!=\[Epsilon]kplusq,thisSusRe=(f1-f2)/(\[Epsilon]k-\[Epsilon]kplusq);];

thisSusIm=dFD[\[Epsilon]k]*dFD[\[Epsilon]kplusq];

thisd2EPW=thisgEPW^2*thisSusRe;
thisd2W90=thisg2W90*thisSusRe;

sus=sus+thisSusRe;
susIm=susIm+thisSusIm;
g2SumEPW=g2SumEPW+thisgEPW^2;
g2SumW90=g2SumW90+thisg2W90;
d2SumEPW=d2SumEPW+thisd2EPW;
d2SumW90=d2SumW90+thisd2W90;
,{indexq,1,indexqmax}];

timeStamp=DateString[{"ISOOrdinalDate","-","Hour24","-","Minute","Second"}];
WriteString[file,indexk,"\t",thiskx,"\t",thisky,"\t",CForm[g2SumEPW],"\t",CForm[g2SumW90],"\t",CForm[sus],"\t",CForm[susIm],"\t",CForm[d2SumEPW],"\t",CForm[d2SumW90],"\t",timeStamp,"\t","\n"]

,{indexk,1,indexkmax}];



ParallelEvaluate[Close@file];

outputCoreData={};
fileNumbers=Length[Kernels[]];
For[i=1,i<=fileNumbers,i++,
fileName="outputCore"<>ToString[i]<>".txt";
thisoutputCoreData=Import[fileName,"Table"];
outputCoreData=Join[outputCoreData,thisoutputCoreData]
]

timeList=Flatten[Take[outputCoreData,All,{-1}]];
outputCoreData=Drop[outputCoreData,None,-1];
outputCoreData=ToExpression[outputCoreData];

indexkList=Flatten[Take[outputCoreData,All,{1}]];
kxList=Flatten[Take[outputCoreData,All,{2}]];
kyList=Flatten[Take[outputCoreData,All,{3}]];
g2SumEPWList=Flatten[Take[outputCoreData,All,{4}]];
g2SumW90List=Flatten[Take[outputCoreData,All,{5}]];
susList=Flatten[Take[outputCoreData,All,{6}]];
susImList=Flatten[Take[outputCoreData,All,{7}]];
d2SumEPWList=Flatten[Take[outputCoreData,All,{8}]];
d2SumW90List=Flatten[Take[outputCoreData,All,{9}]];

order=Ordering[indexkList];
indexkListOrdered=indexkList[[order]];
kxListOrdered=kxList[[order]];
kyListOrdered=kyList[[order]];
g2SumEPWListOrdered=g2SumEPWList[[order]];
g2SumW90ListOrdered=g2SumW90List[[order]];
susListOrdered=susList[[order]];
susImListOrdered=susImList[[order]];
d2SumEPWListOrdered=d2SumEPWList[[order]];
d2SumW90ListOrdered=d2SumW90List[[order]];
timeListOrdered=timeList[[order]];
finalOutput=Transpose[{indexkListOrdered,kxListOrdered,kyListOrdered,g2SumEPWListOrdered,g2SumW90ListOrdered,susListOrdered,susImListOrdered,d2SumEPWListOrdered,d2SumW90ListOrdered}];
finalOutputHeader={"indexk","kx","ky","g2EPW","g2W90","susRe","susIm","d2EPW","d2W90"};
finalOutput=Prepend[finalOutput,finalOutputHeader];

Export["outputFinalIndexOfk.txt",indexkListOrdered];
Export["outputFinalkx.txt",kxListOrdered];
Export["outputFinalky.txt",kyListOrdered];
Export["outputFinalSusRe.csv",susListOrdered];
Export["outputFinalSusIm.csv",susImListOrdered];
Export["outputFinalg2SumEPW.csv",g2SumEPWListOrdered];
Export["outputFinalg2SumW90.csv",g2SumW90ListOrdered];
Export["outputFinald2SumEPW.csv",d2SumEPWListOrdered];
Export["outputFinald2SumW90.csv",d2SumW90ListOrdered];
Export["outputFinalTimeSequence.txt",timeListOrdered];
Export["finalOutput.csv",finalOutput];



