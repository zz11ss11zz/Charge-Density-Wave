(* ::Package:: *)

LaunchKernels[16]; (*parallel 16 cores*)

inputkxFileName="../inputs/kx.txt";
inputkyFileName="../inputs/ky.txt";
inputEnergyFileName="../inputs/energy.txt";

initialOutput={};

T=0.005; (*electronic temperature*)
fermiDirac[x_]:=1/(Exp[x/T]+1);



kx=ReadList[OpenRead[inputkxFileName],Number];
ky=ReadList[OpenRead[inputkyFileName],Number];
band=ReadList[OpenRead[inputEnergyFileName],Number];

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
thisqx=kx[[indexq]];
thisqy=ky[[indexq]];
thisqvector={thisqx,thisqy};
sus=0;


For[indexk=1,indexk<=indexkmax,indexk++,
thiskx=kx[[indexk]];
thisky=ky[[indexk]];
thiskvector={thiskx,thisky};

thiskplusqx=thiskx+thisqx;
thiskplusqy=thisky+thisqy;
thiskplusqvector={thiskplusqx,thiskplusqy};


If[RegionMember[bz,thiskplusqvector]==False,
While[thiskplusqvector[[1]]>Max[kx],thiskplusqvector=thiskplusqvector-basisvector1];
While[thiskplusqvector[[1]]<Min[kx],thiskplusqvector=thiskplusqvector+basisvector1];
];

If[thiskplusqvector[[2]]>Max[ky],While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector-basisvector2]];
If[thiskplusqvector[[2]]<Min[ky],While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector+basisvector2]];

If[thiskplusqvector[[2]]>upperLeftky,While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector-basisvector2]];
If[thiskplusqvector[[2]]<lowerRightky,While[RegionMember[bz,thiskplusqvector]==False,thiskplusqvector=thiskplusqvector+basisvector2]];

nearestkplusqvector=Nearest[kvector,thiskplusqvector][[1]];

indexNearestkplusqvector=Position[kvector,nearestkplusqvector][[1]][[1]];


\[Epsilon]k=band[[indexk]];
\[Epsilon]kplusq=band[[indexNearestkplusqvector]];


f1=fermiDirac[\[Epsilon]k];
f2=fermiDirac[\[Epsilon]kplusq];

If[\[Epsilon]k==\[Epsilon]kplusq,thiskterm=0;];
If[\[Epsilon]k!=\[Epsilon]kplusq,thiskterm=(f1-f2)/(\[Epsilon]k-\[Epsilon]kplusq);];

sus=sus+thiskterm;

];



timeStamp=DateString[{"ISOOrdinalDate","-","Hour24","-","Minute","Second"}];
WriteString[file,indexq,"\t",thisqx,"\t",thisqy,"\t",NumberForm[sus,20],"\t",timeStamp,"\n"]

,{indexq,1,indexqmax}];

ParallelEvaluate[Close@file];

outputCoreData={};
fileNumbers=Length[Kernels[]];
For[i=1,i<=fileNumbers,i++,
fileName="outputCore"<>ToString[i]<>".txt";
thisOutputCoreData=StringSplit[ReadList[OpenRead[fileName],String]];
outputCoreData=Join[outputCoreData,thisOutputCoreData];
]

timeList=Flatten[Take[outputCoreData,All,{-1}]];
outputCoreData=Drop[outputCoreData,None,-1];
outputCoreData=ToExpression[outputCoreData];

indexqList=Flatten[Take[outputCoreData,All,{1}]];
qxList=Flatten[Take[outputCoreData,All,{2}]];
qyList=Flatten[Take[outputCoreData,All,{3}]];
susList=Flatten[Take[outputCoreData,All,{4}]];

order=Ordering[indexqList];
indexqListOrdered=indexqList[[order]];
qxListOrdered=qxList[[order]];
qyListOrdered=qyList[[order]];
susListOrdered=susList[[order]];
timeListOrdered=timeList[[order]];

Export["outputFinalIndexOfq.txt",indexqListOrdered];
Export["outputFinalqx.txt",qxListOrdered];
Export["outputFinalqy.txt",qyListOrdered];
Export["outputFinalSus.txt",susListOrdered];
Export["outputFinalTimeSequence.txt",timeListOrdered];







