(* ::Package:: *)
#!/usr/local/bin/MathematicaScript -script (*.m: linux version*)

LaunchKernels[16];  (*parallel 16 cores*)

inputkxFileName="../inputs/kx.txt";
inputkyFileName="../inputs/ky.txt";
inputEnergyFileName="../inputs/energy.txt";

initialOutput={};

(*select smearing used in the delta function*)
criticalValue1=10^-3;
criticalValue2=10^-2;
criticalValue3=10^-1;



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
susIm0=0;
susIm1=0;
susIm2=0;
susIm3=0;



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

susImthiskterm0=Boole[Abs[\[Epsilon]k]==0]*Boole[Abs[\[Epsilon]kplusq]==0];
susImthiskterm1=Boole[Abs[\[Epsilon]k]<=criticalValue1]*Boole[Abs[\[Epsilon]kplusq]<=criticalValue1];
susImthiskterm2=Boole[Abs[\[Epsilon]k]<=criticalValue2]*Boole[Abs[\[Epsilon]kplusq]<=criticalValue2];
susImthiskterm3=Boole[Abs[\[Epsilon]k]<=criticalValue3]*Boole[Abs[\[Epsilon]kplusq]<=criticalValue3];

susIm0=susIm0+susImthiskterm0;
susIm1=susIm1+susImthiskterm1;
susIm2=susIm2+susImthiskterm2;
susIm3=susIm3+susImthiskterm3;


];



timeStamp=DateString[{"ISOOrdinalDate","-","Hour24","-","Minute","Second"}];
WriteString[file,indexq,"\t",thisqx,"\t",thisqy,"\t",susIm0,"\t",susIm1,"\t",susIm2,"\t",susIm3,"\t",timeStamp,"\n"]
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
susIm0List=Flatten[Take[outputCoreData,All,{4}]];
susIm1List=Flatten[Take[outputCoreData,All,{5}]];
susIm2List=Flatten[Take[outputCoreData,All,{6}]];
susIm3List=Flatten[Take[outputCoreData,All,{7}]];

order=Ordering[indexqList];
indexqListOrdered=indexqList[[order]];
qxListOrdered=qxList[[order]];
qyListOrdered=qyList[[order]];
susIm0ListOrdered=susIm0List[[order]];
susIm1ListOrdered=susIm1List[[order]];
susIm2ListOrdered=susIm2List[[order]];
susIm3ListOrdered=susIm3List[[order]];
timeListOrdered=timeList[[order]];

Export["outputFinalIndexOfq.txt",indexqListOrdered];
Export["outputFinalqx.txt",qxListOrdered];
Export["outputFinalqy.txt",qyListOrdered];
Export["outputFinalSusIm0.txt",susIm0ListOrdered];
Export["outputFinalSusIm1.txt",susIm1ListOrdered];
Export["outputFinalSusIm2.txt",susIm2ListOrdered];
Export["outputFinalSusIm3.txt",susIm3ListOrdered];
Export["outputFinalTimeSequence.txt",timeListOrdered];







