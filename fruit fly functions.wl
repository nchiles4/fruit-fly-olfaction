(* ::Package:: *)

(* ::Subtitle:: *)
(*Two Bloom Filter Initialization*)


(* BeginPackage["pos`"];
pos::usage=
"pos returns a list of novelty scores for a set of queries evaluated by a Fly Bloom Filter 
trained only with active molecules."

pos[actData_,inactData_][train_,test_,iterations_]:=ParallelTable[
Total[twoFilter[sampleSplit[actData,train+test,train],sampleSplit[inactData,train,train]]]/test,iterations];
EndPackage[]; *)


(* BeginPackage["neg`"];
neg::usage=
"neg returns a list of novelty scores for a set of queries evaluated by a Fly Bloom Filter 
trained only with inactive molecules."

neg[actData_,inactData_][train_,test_,iterations_]:=ParallelTable[
Total[twoFilter[sampleSplit[actData,train,train],sampleSplit[inactData,train+test,train]]]/test,iterations];
EndPackage[]; *)


(* Begin["`Private`"]; *)

(* defines random sparse binary matrix *)
M[size_:2000,numberOfBits_:6,inputLength_:2048]:=Table[
RandomSample[Join[Table[1,numberOfBits],
Table[0,size-numberOfBits]]],inputLength];

(* reduces dimensions of initial fingerprint by taking product of fp and random matrix *)
projectM[statRandM_][traindata_,m_:2000,k_:25]:=ReplacePart[
Table[1,m],#->0&/@Ordering[traindata . statRandM,-k]];
consolidateM[statRandM_?MatrixQ,training_List]:=projectM[statRandM]/@training//Total

(* unitizes output *)
flattenM[consolidateM_,n_]:=EqualTo[n]/@consolidateM//Boole

(* recalls flattening functions given training data size *)
inputFunction[statRandM_,traindata_List]:=flattenM[consolidateM[statRandM,traindata],Length[traindata]];

(* defines functions to compare outputs between trained filter and testing fingerprints *)
novelty[statRandM_,train_List][test_]:=With[
{tr=inputFunction[statRandM,train],
te=inputFunction[statRandM,{test}]},
 Mean@Pick[tr,te,0]]

(* creates two filters trained with different data, outputs 1 if higher similarity with actives, 
0 if higher similarity with inactives *)
twoFilter[splitDataActive_,splitDataInactive_]:=With[
{staticRandomMatrix=M[],
activetrain=splitDataActive[[1]],
inactivetrain=splitDataInactive[[1]],
test=Join[splitDataActive[[2]],splitDataInactive[[2]]]},
NonNegative[novelty[staticRandomMatrix,inactivetrain]/@test-
novelty[staticRandomMatrix,activetrain]/@test]//Boole]
(* EndPackage[]; *)


pos[actData_,inactData_][train_,test_,iterations_]:=ParallelTable[
Total[twoFilter[sampleSplit[actData,train+test,train],sampleSplit[inactData,train,train]]]/test,iterations]; 
neg[actData_,inactData_][train_,test_,iterations_]:=ParallelTable[
Total[twoFilter[sampleSplit[actData,train,train],sampleSplit[inactData,train+test,train]]]/test,iterations];




(* ::Subtitle:: *)
(*Splitting Sample Data*)


(* splitting data into training and testing groups *)
sampleSplit[Data_,totalSize_,n_]:=With[
{randomSample=RandomSample[Data,totalSize]},
{randomSample[[;;n]],randomSample[[(n+1);;]]}](*{train,test}*)


sampleSplit2[activeData_,inactiveData_,n_]:=With[
{a=ResourceFunction["TrainTestSplit"][activeData,"TrainingSetSize"->n],
i=ResourceFunction["TrainTestSplit"][inactiveData,"TrainingSetSize"->n]},
<|"train"->Join[a[[1]],i[[1]]],"test"->Join[a[[2]],i[[2]]]|>]


(* ::Subtitle:: *)
(*Traditional Classifier Setup*)


(* easy input function for traditional classifiers *)
evaluateModel[splitData_Association,method_String]:=With[
{model=Classify[splitData["train"],Method->method]},
ClassifierMeasurements[model,splitData["test"],{"Accuracy","Recall"}]]


booleTrue[active_]:=active-> True;
booleFalse[inactive_]:=inactive->False;
boolePos[actData_][trainplustest_]:=Map[booleTrue,RandomSample[actData,trainplustest]];
booleNeg[inactData_][trainplustest_]:=Map[booleFalse,RandomSample[inactData,trainplustest]];


(* easy input train and test size for traditional classifiers *)
stats[modelType_String, actData_,inactData_][train_,test_]:=With[
{learningdata=sampleSplit2[boolePos[actData][train+test],booleNeg[inactData][train+test],train]},
evaluateModel[learningdata,modelType]]


classStats[modelType_String,actData_,inactData_][train_,test_,iterations_]:=ParallelTable[
stats[modelType, actData,inactData][train,test],iterations]//Transpose;


(* ::Subtitle:: *)
(*Box Whisker Charts*)


accChart[bloom_,dt_,rf_,lr_,nn_]:=BoxWhiskerChart[
{Join[bloom[[1]],1-bloom[[2]]],dt[[1]],rf[[1]],lr[[1]],nn[[1]]},
ChartStyle->{Lighter[Red],Darker[Green],Darker[Cyan],Lighter[Purple],Darker[Purple]},
LabelStyle->{FontSize->14,FontFamily->"Times",Black,Bold},
PlotRange->{.1,0.9}]


recChart[bloom_,dt_,rf_,lr_,nn_]:=BoxWhiskerChart[
{bloom[[1]],
Table[dt[[2,x,2]],{x,1,25}],
Table[rf[[2,x,2]],{x,1,25}],
Table[lr[[2,x,2]],{x,1,25}],
Table[nn[[2,x,2]],{x,1,25}]},
ChartStyle->{Lighter[Red],Darker[Green],Darker[Cyan],Lighter[Purple],Darker[Purple]},
LabelStyle->{FontSize->14,
FontFamily->"Times",Black,Bold},
PlotRange->{.1,0.9}];


horiz=Plot[0.5,{x,0,6},PlotStyle->Black];


(* ::Subtitle:: *)
(*F1 Scores*)


flyConfusionMatrix[data_]:={data[[1]],1-data[[1]],data[[2]],1-data[[2]]};


confusionMatrix[data_]:={Table[data[[2,x,2]],{x,1,25}],
1-Table[data[[2,x,2]],{x,1,25}],
1-Table[data[[2,x,1]],{x,1,25}],
Table[data[[2,x,1]],{x,1,25}]};


F1[CF_]:=CF[[1]]/(CF[[1]]+.5(CF[[3]]+CF[[2]]))


F1List[fly_,dt_,rf_,lr_,nn_]:={F1[flyConfusionMatrix[fly]],F1[confusionMatrix[dt]],F1[confusionMatrix[rf]],F1[confusionMatrix[lr]],F1[confusionMatrix[nn]]};


F1Chart[allData_]:=BoxWhiskerChart[
allData,
ChartStyle->{Lighter[Red],Darker[Green],Darker[Cyan],Lighter[Purple],Darker[Purple]},
LabelStyle->{FontSize->14,FontFamily->"Times",Black,Bold},
PlotRange->{0.1,0.9}];


(* SetDirectory[]@NotebookDirectory[] *)
