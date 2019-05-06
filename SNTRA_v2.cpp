#include <stdio.h>
#include "InputData.hh"
//#include "VectorLib.cpp"
#include "Experiment.cpp"
//#include <vector>
//#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>

//#include "ExtStrings.cpp"

#include "PrimitiveShellModel.cpp"
//root-зависимые библиотеки

#include <TGraph.h>
#include <TFile.h>
#include <TSystem.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TLine.h>
#include <TMultiGraph.h>



using namespace std;

TCanvas *cc1=new TCanvas("cc1","cc1");




/*ostream& operator << (ostream &s, CoupleOfExperiments &CE)
{
	s<<"Pickup:"<<CE.Pickup.reference<<" Stripping:"<<CE.Stripping.reference<<"\n";
	s<<CE.SPE.size()<<" "<<CE.EJP.size()<<" "<<CE.n.size()<<"\n";
	for(unsigned int i=0;i<CE.EJP.size();i++)
	{
		
		//s<<CE.SPE[i]<<" "<<CE.EJP[i]<<" "<<CE.n[i]<<"\n";
	}
}*/

vector<string> ListFiles(string mask)
{
	vector<string> FileNames;
	string s;
	FILE* fp;
	char result [1000];
	fp = popen(("ls "+mask).c_str(),"r");
	fread(result,1,sizeof(result),fp);
	s=result;
	fclose (fp);
	stringstream ss(s);
	while(ss)
	{
		s.resize(0);
		ss>>s;
		//if((GetFileName(s)[0]>='.')&&(GetFileName(s)[0]<='z'))
		FileNames.push_back(s);
		//cout<<s<<"\n";
	}
	return FileNames;
}

vector<string> ListFiles(string dirname, string ext) 
{
	TSystemDirectory dir(dirname.c_str(), dirname.c_str()); 
	TList *files = dir.GetListOfFiles(); 
	vector<string> result;
	if(files) 
	{
		TSystemFile *file; TString fname; 
		TIter next(files);
		while((file=(TSystemFile*)next()))
		{
			fname = file->GetName(); 
			if(!file->IsDirectory() && fname.EndsWith(ext.c_str()))
			{ 
				result.push_back(dirname+(string)fname); 
				cout<<(string)fname<<"\n";
			} 
		} 
	}
	return result; 
}

vector<CoupleOfExperiments> CreateCouplesOfExperiments(vector<Experiment> &Pickup,vector<Experiment> &Stripping,parameters &par)
{
	vector<CoupleOfExperiments> result;
	for(unsigned int i=0;i<Pickup.size();i++)
	{
		Pickup[i].CalculateNPH(par);
		if(par.Normalize())
		{
			cout<<"Pickup_normalization1:"<<Pickup[i].GetNumberOfParticles()<<" "<<(float)par.NumberOfParticlesInUsedShell<<" "<<(float)par.NumberOfParticlesInUsedShell/Pickup[i].GetNumberOfParticles()<<"\n";
			Pickup[i].Scale((float)par.NumberOfParticlesInUsedShell/Pickup[i].GetNumberOfParticles());
			Pickup[i].CalculateNPH(par);
			cout<<"Pickup_normalization2:"<<Pickup[i].GetNumberOfParticles()<<" "<<(float)par.NumberOfParticlesInUsedShell<<"\n";
		}
	}
	for(unsigned int j=0;j<Stripping.size();j++)
	{
		Stripping[j].CalculateNPH(par);
		if(par.Normalize())
		{
			Stripping[j].Scale((float)par.NumberOfHolesInUsedShell/Stripping[j].GetNumberOfHoles());
		}
	}
	for(unsigned int i=0;i<Pickup.size();i++)
	{
		//cout<<i<<" pickup "<<Pickup[i].reference<<"\n";
		for(unsigned int j=0;j<Stripping.size();j++)
		{
			//cout<<j<<" stripping "<<Stripping[j].reference<<"\n";
			CoupleOfExperiments CE(Pickup[i],Stripping[j]);
			CE.par=par;
			result.push_back(CE);
		}
	}
	return result;
}

void ReadFilesInDirectory(string PathToFiles,vector<Experiment> &Pickup,vector<Experiment> &Stripping,string particle, int ListFilesFlag=0)
{
	vector<string> FileNames;
	if(ListFilesFlag==0)
	{
		FileNames=ListFiles(PathToFiles);
	}
	else
	{
		stringstream ss(PathToFiles);
		string path,ext;
		ss>>path;
		ss>>ext;
		FileNames=ListFiles(path,ext); 
	}
	////cout<<"size="<<FileNames.size()<<"\n";
	for(unsigned int i=0;i<FileNames.size();i++)
	{
		Experiment E;
		E.ReadInputFile(FileNames[i]);
		//cout<<E.GetType()<<"\n";
		int Z,A,ParticleType;//ParticleType-передача нейтрона (0) или протона(1)
		string type;//pickup/stripping
		GetAZ(E.Nucleus,Z,A);
		//cout<<E.Nucleus<<" "<<Z<<" "<<A<<"\n";
		ParceReaction(E.reaction,type,ParticleType);
				
		if(ParticleType==1)
		{
			E.BA1=GetSeparationEnergy(Z+1, A+1, 1,1)*1000;
			E.BA=GetSeparationEnergy(Z, A, 1,1)*1000;
			E.particle="proton";
		}
		else
		{
			E.BA1=GetSeparationEnergy(Z, A+1, 0,1)*1000;
			E.BA=GetSeparationEnergy(Z, A, 0,1)*1000;
			E.particle="neutron";
		}
		
		E.ProcessExperimentalData();
		
		if((E.particle==particle)||(particle==""))
		{
			if(E.GetType()=="pickup")
			{
				for(unsigned int j=0;j<Pickup.size();j++)
				{
					if(Pickup[j].reference==E.reference)
					{
						E.ErrorString="Non-unique reference found!\n"+E.ErrorString+"\n";
						E.reference=FileNames[i];
					}
				}
				Pickup.push_back(E);
			}
			else if(E.GetType()=="stripping")
			{
				for(unsigned int j=0;j<Stripping.size();j++)
				{
					if(Stripping[j].reference==E.reference)
					{
						E.ErrorString="Non-unique reference found!\n"+E.ErrorString+"\n";
						E.reference=FileNames[i];
					}
				}
				Stripping.push_back(E);
			}
			SplitExperiments(Pickup);
			SplitExperiments(Stripping);	
		}
	}
}



void CalculatePenaltyFunction(vector<CoupleOfExperiments> &v)
{
	float MaxEfError,MaxDeltaError;
	int NumberOfPickupStatesMax, NumberOfStrippingStatesMax, AverageNumberOfCalculatedStates;
	NumberOfPickupStatesMax=0;
	NumberOfStrippingStatesMax=0;
	AverageNumberOfCalculatedStates=0;
	float MaxNParticles=0;
	float MaxNHoles=0;
	for(int i=0;i<v.size();i++)
	{
		if(MaxEfError<v[i].Ef_error)
		{
			MaxEfError=v[i].Ef_error;
		}
		if(MaxDeltaError<v[i].Delta_error)
		{
			MaxDeltaError=v[i].Delta_error;
		}
		if(NumberOfPickupStatesMax<v[i].GetNumberOfPickupStates())//GetNlevels())
		{
			//NumberOfPickupStatesMax=v[i].Pickup.GetNlevels();
			NumberOfPickupStatesMax=v[i].GetNumberOfPickupStates();
			cout<<"v[i].GetNumberOfPickupStates():"<<v[i].GetNumberOfPickupStates()<<"\n";
		}
		if(NumberOfStrippingStatesMax<v[i].GetNumberOfStrippingStates())//GetNlevels())
		{
			NumberOfStrippingStatesMax=v[i].GetNumberOfStrippingStates();//)//GetNlevels();
			cout<<"v[i].GetNumberOfStrippingStates():"<<v[i].GetNumberOfStrippingStates()<<"\n";
		}
		AverageNumberOfCalculatedStates+=v[i].SPE.size();
		float NP=abs(v[i].par.NumberOfParticlesInUsedShell-v[i].NumberOfParticlesInUsedShell);
		float NH=abs(v[i].par.NumberOfHolesInUsedShell-v[i].NumberOfHolesInUsedShell);
		if(MaxNParticles<NP)
		{
			MaxNParticles=NP;
		}
		if(MaxNHoles<NH)
		{
			MaxNHoles=NH;
		}
	}
	if(v.size()>0)
	AverageNumberOfCalculatedStates=round(AverageNumberOfCalculatedStates/v.size());
	else
	return;
	
	for(unsigned int i=0;i<v.size();i++)
	{
		//v[i].PenaltyComponents.resize(v[i].par.UsedPenaltyFunctionComponents.size());
		//cout<<"size "<<v[i].par.UsedPenaltyFunctionComponents.size()<<"\n";
		for(unsigned int j=0;j<v[i].par.UsedPenaltyFunctionComponents.size();j++)
		{
			//cout<<"comp "<<(int)v[i].par.UsedPenaltyFunctionComponents[j]<<"\n";
			if(v[i].par.UsedPenaltyFunctionComponents[j]==1)
			{
				v[i].PenaltyComponents.push_back(abs(1-Average(v[i].ParticlesAndHolesSum)));
			}
			else if(v[i].par.UsedPenaltyFunctionComponents[j]==2)
			{
				v[i].PenaltyComponents.push_back(1-((float)v[i].GetNumberOfPickupStates()/NumberOfPickupStatesMax));
			}
			else if(v[i].par.UsedPenaltyFunctionComponents[j]==3)
			{
				v[i].PenaltyComponents.push_back(1-((float)v[i].GetNumberOfStrippingStates()/NumberOfStrippingStatesMax));
			}
			
			else if(v[i].par.UsedPenaltyFunctionComponents[j]==4)
			{
				if((MaxEfError!=0)&&(MaxDeltaError!=0)&&(v[i].FitStatus==4000)&&(v[i].OccStatus==0))
				{
					v[i].PenaltyComponents.push_back(v[i].Ef_error/MaxEfError);
				}
				else
				{
					v[i].PenaltyComponents.push_back(MaxEfError);
				}
			}
			else if(v[i].par.UsedPenaltyFunctionComponents[j]==5)
			{
				if((MaxEfError!=0)&&(MaxDeltaError!=0)&&(v[i].FitStatus==4000)&&(v[i].OccStatus==0))
				{
					v[i].PenaltyComponents.push_back(v[i].Delta_error/MaxDeltaError);
				}
				else
				{
					v[i].PenaltyComponents.push_back(MaxDeltaError);
				}
			}
			else if(v[i].par.UsedPenaltyFunctionComponents[j]==6)
			{
				float N=abs(v[i].par.NumberOfParticlesInUsedShell-v[i].NumberOfParticlesInUsedShell)/MaxNParticles;
				cout<<"v[i].par.NumberOfParticlesInUsedShell:"<<v[i].par.NumberOfParticlesInUsedShell<<" v[i].NumberOfParticlesInUsedShell:"<<v[i].NumberOfParticlesInUsedShell<<" MaxNParticles:"<<MaxNParticles<<"\n";
				v[i].PenaltyComponents.push_back(N);
			}
			else if(v[i].par.UsedPenaltyFunctionComponents[j]==7)
			{
				float N=abs(v[i].par.NumberOfHolesInUsedShell-v[i].NumberOfHolesInUsedShell)/MaxNHoles;
				cout<<"v[i].par.NumberOfHolesInUsedShell:"<<v[i].par.NumberOfHolesInUsedShell<<" v[i].NumberOfHolesInUsedShell:"<<v[i].NumberOfHolesInUsedShell<<" MaxNHoles:"<<MaxNHoles<<"\n";
				v[i].PenaltyComponents.push_back(N);
			}
			//cout<<"v[i].par.weights[j]:"<<v[i].par.weights[j]<<"\n";
		}
		
		//cout<<v[i].Pickup.reference<<" "<<v[i].Stripping.reference<<"\n";
		float weights_sum=0;
		for(unsigned int j=0;j<v[i].PenaltyComponents.size();j++)
		{
			//cout<<"p["<<j<<"]="<<v[i].PenaltyComponents[j]<<" "<<Average(v[i].ParticlesAndHolesSum)<<"\n";
			v[i].PenaltyComponents[j]=v[i].PenaltyComponents[j]*v[i].par.weights[j];
			weights_sum+=v[i].par.weights[j];
			v[i].penalty+=v[i].PenaltyComponents[j];
		}
		v[i].penalty=v[i].penalty/weights_sum;
	}
	
}

void PrintCalculationResult(vector<CoupleOfExperiments> v,string OutputFileName)
{
	ofstream OutputTextFile((OutputFileName+".txt").c_str());
	cc1->Print((OutputFileName+".pdf[").c_str(),"pdf");
	for(unsigned int i=0;i<v.size();i++)
	{
		SpectroscopicFactorHistogram HistPickup=v[i].Pickup.BuildSpectroscopicFactorHistogram();
		SpectroscopicFactorHistogram HistStrip=v[i].Stripping.BuildSpectroscopicFactorHistogram();
		cc1->Clear();
		cc1->Divide(3,2);
		cc1->cd(1);
		//gPad->SetLogy(1);
		HistPickup.PrintSpectroscopicFactorHistogram();
		cc1->cd(2);
		TMultiGraph mgr;
		v[i].occupancies.SetTitle("Occupancy;E,keV;v^2");
		mgr.Add(&v[i].Pickup_occupancies);
		mgr.Add(&v[i].Stripping_occupancies);
		mgr.Add(&v[i].Both_occupancies);
		mgr.Draw("ap");
		v[i].occupancies.Draw("p");
		v[i].BCS.Draw("l same");
		cc1->cd(3);
		TH1F PenaltyComponents=v[i].BuildPenaltyComponentsHistogram();
		gPad->SetLogy(1);
		PenaltyComponents.Draw();
		cc1->cd(4);
		//gPad->SetLogy(1);
		HistStrip.PrintSpectroscopicFactorHistogram();
		string TextOutput=v[i].ResultsInTextForm(1);
		stringstream s(TextOutput);
		OutputTextFile<<TextOutput<<"\n";
		
		if(v[i].Pickup.ErrorString.size()>0)
		OutputTextFile<<"Errors (pickup): "<<v[i].Pickup.ErrorString<<"\n";
		if(v[i].Stripping.ErrorString.size()>0)
		OutputTextFile<<"Errors (stripping): "<<v[i].Stripping.ErrorString<<"\n";
		OutputTextFile<<"Penalty:\n";
		for(unsigned int j=0;j<v[i].PenaltyComponents.size();j++)
		{
			OutputTextFile<<v[i].par.GetComponentName(j)<<": "<<v[i].PenaltyComponents[j]<<"\n";
		}
		
		cc1->cd(5);
		TGraph* gr=new TGraph();//"h1","Calculated shell scheme;1 ;E, keV",10,0,1);
		gr->SetPoint(0,0,0);
		gr->SetMinimum(GetMinimum(v[i].SPE)-1000);
		gr->SetMaximum(GetMaximum(v[i].SPE)+1000);
		
		gr->Draw("ap");
		for(unsigned int j=0;j<v[i].SPE.size();j++)
		{
			TLine line;
			TText txt;
			int color=v[i].SP_centroids[j].GetColor();
			line.SetLineColor(color);
			line.DrawLine(0.1,v[i].SPE[j],0.7,v[i].SPE[j]);
			txt.SetTextColor(color);
			txt.DrawText(0.8,v[i].SPE[j], NLJToString(v[i].SP_centroids[j].n,v[i].SP_centroids[j].l,v[i].SP_centroids[j].JP).c_str());
		}
		
		cc1->cd(6);
		v[i].DrawResultsInTextForm();
		cc1->Print((OutputFileName+".pdf").c_str(),"pdf");
		gPad->Update();
		
	}
	cc1->Print((OutputFileName+".pdf]").c_str(),"pdf");
}

void ArrangeByPenalty(vector<CoupleOfExperiments> &v)
{
	for(unsigned int i=0;i<v.size();i++)
	{
		int NumberOfExcanges=0;
		for(unsigned int j=0;j<v.size()-i-1;j++)
		{
			if(v[j].penalty>v[j+1].penalty)
			{
				CoupleOfExperiments tmp=v[j];
				v[j]=v[j+1];
				v[j+1]=tmp;
				NumberOfExcanges++;
			}
		}
		if(NumberOfExcanges==0)
		{
			return;
		}
	}
}

//int main()
void SNTRA_v2(string PathToFiles, string particle="", int ListFilesFlag=0)
{
	vector<Experiment> Pickup;
	vector<Experiment> Stripping;
	
	ReadFilesInDirectory(PathToFiles,Pickup,Stripping,particle,ListFilesFlag);
	parameters par;
	stringstream s1(PathToFiles);
	string ParFileName;
	s1>>ParFileName;
	par.ReadParameters(ParFileName+"parameters.par");
	vector<CoupleOfExperiments> CE=CreateCouplesOfExperiments(Pickup,Stripping,par);
	//TCanvas *cc1=new TCanvas("cc1","cc1");
	for(unsigned int i=0;i<CE.size();i++)
	{
		CE[i].CalcSPE_and_OCC(cc1);
		//cout<<CE[i].SPE.size()<<" "<<CE[i].EJP.size()<<" "<<CE[i].n.size()<<"\n";
		
		for(int j=0;j<CE[i].SPE.size();j++)
		{
			//cout<<CE[i].SPE[j]<<" "<<CE[i].n[j]<<" "<<CE[i].EJP[j]<<"\n"; 
		}
		
		////cout<<CE[i]<<"\n";
	}
	string OutputFileName;
	if((Pickup.size()>0)&&(Stripping.size()>0))
	{
		OutputFileName=Pickup[0].Nucleus+"_"+Pickup[0].particle;		
	}
	else
	{
		return ;
	}
	CalculatePenaltyFunction(CE);
	ArrangeByPenalty(CE);
	PrintCalculationResult(CE,OutputFileName);
}

int main(int argc, char** argv)
{
	string path=argv[1];
	string ext=argv[2];
	cout<<path+" "+ext<<"\n";
	SNTRA_v2(path+" "+ext,"",1);
}
