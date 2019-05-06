#include "Experiment.hh"
#include "TFitResult.h"
void SetTGraphLimits(TGraph &gr,float xmin,float xmax,float ymin, float ymax)
{
	gr.Draw("AP");
	gr.GetXaxis()->SetLimits(xmin*0.8,xmax*1.2);
	gr.SetMinimum(ymin*0.8);
	gr.SetMaximum(ymax*1.2);
	gr.SetMarkerStyle(8);
}

void SpectroscopicFactorHistogram::PrintSpectroscopicFactorHistogram()
{
	Legend=new TLegend(0.7,0.7,0.9,0.9);
	for(unsigned int j=0;j<Histograms.size();j++)
	{
		if(j==0)
		{
			Histograms[j].SetTitle((Reference+";E,kev;G").c_str());
			Histograms[j].Draw("histo");
			Legend->AddEntry(&Histograms[j],NLJToString(n[j],L[j], JP[j]).c_str(),"f");
			Histograms[j].GetYaxis()->SetRangeUser(10e-9,maximum*1.5);
		}
		else
		{
			Legend->AddEntry(&Histograms[j],NLJToString(n[j],L[j], JP[j]).c_str(),"f");
			Histograms[j].Draw("histo same");
		}
		Legend->Draw();
	}
}
StateParameters::StateParameters()
{
	
}
StateParameters::StateParameters(int n,int l,double JP,string c_flag)
{
	this->n=n;
	this->l=l;
	this->JP=JP;
	if(c_flag=="pickup")
	{
		couple_flag=1;
	}
	else if(c_flag=="stripping")
	{
		couple_flag=2;
	}
	else if(c_flag=="both")
	{
		couple_flag=3;
	}
	else
	{
		couple_flag=0;
	}
}
unsigned char StateParameters::GetColor()
{
	if(couple_flag==3)
	{
		return 1;
	}
	if(couple_flag==2)
	{
		return 2;
	}
	if(couple_flag==1)
	{
		return 4;
	}
}
bool StateParameters::CompareQN(StateParameters &s)
{
	if((n==s.n)&&(l==s.l)&&(JP==s.JP))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool StateParameters::CompareQN(int n,int l,float JP)
{
	if((n==this->n)&&(l==this->l)&&(abs(JP)==abs(this->JP)))
	{
		return true;
	}
	else
	{
		return false;
	}
}
string StateParameters::GetType()
{
	if(couple_flag==3)
	{
		return "both";
	}
	if(couple_flag==2)
	{
		return "stripping";
	}
	if(couple_flag==1)
	{
		return "pickup";
	}
}
bool parameters::Normalize()
{
	if(NormalizeFlag==1)
	{
		return true;
	}
	else
	{
		return false;
	}
}
void parameters::ReadParameters(string filename)
{

	ifstream ifs(filename.c_str());
	string line;
	while(getline(ifs,line))
	{
		stringstream s(line);
		string tmp;
		s>>tmp;
		if(tmp=="UseIncompleteCouples:")
		{
			s>>tmp;
			if(tmp=="all")
			{
				IncompleteCouplesFlag=1;
			}
			else if(tmp=="pickup")
			{
				IncompleteCouplesFlag=2;
			}
			else if(tmp=="stripping")
			{
				IncompleteCouplesFlag=3;
			}
			else if(tmp=="no")
			{
				IncompleteCouplesFlag=4;
			}
		}
		else if(tmp=="NormalizeToNPH")
		{
			NormalizeFlag=1;
		}
		else if(tmp=="UsedPenaltyFunctionComponents:")
		{
			string tmp2;
			while(s)
			{
				tmp2.resize(0);
				s>>tmp2;
				if(tmp2=="a_ij")
				{
					UsedPenaltyFunctionComponents.push_back(1);
				}
				if(tmp2=="NPickupMax")
				{
					UsedPenaltyFunctionComponents.push_back(2);
				}
				if(tmp2=="NStrippingMax")
				{
					UsedPenaltyFunctionComponents.push_back(3);
				}
				if(tmp2=="EF_err")
				{
					UsedPenaltyFunctionComponents.push_back(4);
				}
				if(tmp2=="Delta_err")
				{
					UsedPenaltyFunctionComponents.push_back(5);
				}
				if(tmp2=="NumberOfParticlesInUsedShell")
				{
					UsedPenaltyFunctionComponents.push_back(6);
					UsedPenaltyFunctionComponents.push_back(7);
				}
			}
		}
		else if(tmp=="SubShellsUsedForOccupancyFit:")
		{
			string tmp2;
			while(s)
			{
				s>>tmp2;
				int n,l;
				float JP;
				if(StringToNLJ(tmp2,n,l,JP))
				{
					StateParameters sp(n,l,JP,"both");
					SubShellsUsedForOccupancyFit.push_back(sp);
				}
			}
		}
		else if(tmp=="SubShellsUsedForPenaltyFunction:")
		{
			string tmp2;
			while(s)
			{
				s>>tmp2;
				int n,l;
				float JP;
				if(StringToNLJ(tmp2,n,l,JP))
				{
					StateParameters sp(n,l,JP,"both");
					SubShellsUsedForPenaltyFunction.push_back(sp);
				}
			}
		}
		else if(tmp=="Weights:")
		{
			while(s)
			{
				float w;
				s>>w;
				weights.push_back(w);
			}
			
		}
		else if(tmp=="NParticles:")
		{
			s>>NumberOfParticlesInUsedShell;
		}
		else if(tmp=="NHoles:")
		{
			s>>NumberOfHolesInUsedShell;
		}
	}
	cout<<"IncompleteCouplesFlag "<<(int)IncompleteCouplesFlag<<"\n";
	for(int i=0;i<UsedPenaltyFunctionComponents.size();i++)
	{
		cout<<i<<"PenaltyFunctionComponents: "<<(int)UsedPenaltyFunctionComponents[i]<<"\n";
	}
	float W_norm=0;
	if(weights.size()>=UsedPenaltyFunctionComponents.size())
	{
		for(unsigned int i=0;i<UsedPenaltyFunctionComponents.size();i++)//нормировка весов
		{
			if(W_norm<weights[i])
			{
				W_norm=weights[i];
			}
		}
	}
	else
	{
		weights.resize(UsedPenaltyFunctionComponents.size());
		for(unsigned int i=0;i<weights.size();i++)//нормировка весов
		{
			weights[i]=1;
		}
	}
	if(W_norm!=0)
	{
		for(unsigned int i=0;i<weights.size();i++)//нормировка весов
		{
			weights[i]=weights[i]/W_norm;
		}
	}
	
}
bool parameters::CheckStateParameters(StateParameters &s)
{
	if((IncompleteCouplesFlag==1)&&(s.couple_flag==1||s.couple_flag==2||s.couple_flag==3))
	{
		return true;
	}
	else if((IncompleteCouplesFlag==2)&&(s.couple_flag==1||s.couple_flag==3))
	{
		return true;
	}
	else if((IncompleteCouplesFlag==3)&&(s.couple_flag==2||s.couple_flag==3))
	{
		return true;
	}
	else if((IncompleteCouplesFlag==4)&&(s.couple_flag==3))
	{
		return true;
	}
	else
	{
		return false;
	}
}
bool parameters::CheckOccupancy(StateParameters &s)
{
	if(SubShellsUsedForOccupancyFit.size()==0)
	{
		return true;//если состояния, для которых фитируется заселенность, не определены, используем все состояния
	}
	for(unsigned int i=0;i<SubShellsUsedForOccupancyFit.size();i++)
	{
		if(SubShellsUsedForOccupancyFit[i].CompareQN(s))
		{
			return true;
		}
	}
	return false;
}
bool parameters::CheckSatesForPenaltyFunction(StateParameters &s)
{
	if(SubShellsUsedForPenaltyFunction.size()==0)
	{
		return true;//если состояния, для которых считается штраф.функция, не определены, используем все состояния
	}
	for(unsigned int i=0;i<SubShellsUsedForPenaltyFunction.size();i++)
	{
		if(SubShellsUsedForPenaltyFunction[i].CompareQN(s))
		{
			return true;
		}
	}
	return false;
}

bool parameters::CheckSatesForPenaltyFunction(State &s)
{
	if(SubShellsUsedForPenaltyFunction.size()==0)
	{
		return true;//если состояния, для которых считается штраф.функция, не определены, используем все состояния
	}
	for(unsigned int i=0;i<SubShellsUsedForPenaltyFunction.size();i++)
	{
		if(SubShellsUsedForPenaltyFunction[i].CompareQN(s.GetN(),s.GetL(),s.GetJP()))
		{
			return true;
		}
	}
	return false;
}

string parameters::PrintUsedSubShells()
{
	stringstream s;
	for(unsigned int i=0;i<SubShellsUsedForOccupancyFit.size();i++)
	{
		s<<"s:"<<SubShellsUsedForOccupancyFit[i].n<<"("<<SubShellsUsedForOccupancyFit[i].l<<")"<<SubShellsUsedForOccupancyFit[i].JP<<" ";
	}
}
string parameters::GetComponentName(unsigned int iterator)
{
	if(iterator<UsedPenaltyFunctionComponents.size())
	{
		if(UsedPenaltyFunctionComponents[iterator]==1)
		{
			return "a_ij";
		}
		if(UsedPenaltyFunctionComponents[iterator]==2)
		{
			return "NPickupMax";
		}
		if(UsedPenaltyFunctionComponents[iterator]==3)
		{
			return "NStrippingMax";
		}
		if(UsedPenaltyFunctionComponents[iterator]==4)
		{
			return "EF_err";
		}
		if(UsedPenaltyFunctionComponents[iterator]==5)
		{
			return "Delta_err";
		}
		if(UsedPenaltyFunctionComponents[iterator]==6)
		{
			return "dNParticles";
		}
		if(UsedPenaltyFunctionComponents[iterator]==7)
		{
			return "dNHoles";
		}
	}
	return "unknown";
}
int Experiment::GetColor(int L, float JP)
{
	if(L==2)
	{
		if(abs(JP)==2.5)
		{
			return kGreen+2;
		}
		if(abs(JP)==1.5)
		{
			return kGreen+1;
		}
	}
	if(L==0)
	{
		return kGreen; 
	}
	if(L==3)
	{
		if(abs(JP)==3.5)
		{
			return kRed+2;
		}
		if(abs(JP)==2.5)
		{
			return kRed+1;
		}
	}
	if(L==1)
	{
		if(abs(JP)==1.5)
		{
			return kRed-6;
		}
		if(abs(JP)==0.5)
		{
			return kRed-7;
		} 
	}
	if(L==4)
	{
		if(abs(JP)==4.5)
		{
			return kRed-10;
		}
		if(abs(JP)==3.5)
		{
			return kRed-9;
		}
	}
}
Experiment::Experiment()
{
	E_iterator=0; n_iterator=1; L_iterator=2; JP_iterator=3; SF_iterator=4;
}
string Experiment::GetType()
{
	if(type==1)
	{
		return "stripping";
	}
	if(type==0)
	{
		return "pickup";
	}
}
void Experiment::ReadInputFile(string filename)//просто чтение файла с данными. Сначала поиск ключевых слов, если их нет, то попытка считать строку как состояние, наблюдаемое в эксперименте
{
	ifstream ifs(filename.c_str());
	string line;
	while(getline(ifs,line))
	{
		//cout<<line<<"\n";
		std::stringstream sstr(line);
		string line_tmp;
		sstr>>line_tmp;
		if((line_tmp=="Reference:")||(line_tmp=="Author:"))
		{
			sstr>>reference;
		}
		if(line_tmp=="Nucleus:")
		{
			sstr>>Nucleus;
		}
		if(line_tmp=="Type:")
		{
			string typeR;
			sstr>>typeR;
			if(typeR=="stripping")
			{
				type=1;
			}
			if(typeR=="pickup")
			{
				type=0;
			}
		}
		else if(line_tmp=="Reaction:")
		{
			sstr>>reaction;
		}
		else if(line_tmp=="ProjectileEnergy:")
		{
			sstr>>ProjectileEnergy;
		}
		else if(line_tmp=="BA:")
		{
			sstr>>BA;
			if(BA<100)
			{
				BA=BA*1000;
			}
		}
		else if(line_tmp=="BA1:")
		{
			sstr>>BA1;
			if(BA1<100)
			{
				BA1=BA1*1000;
			}
		}
		else if(line_tmp=="JP0:")
		{
			sstr>>JP0;
		}
		else
		{
			if((line_tmp[0]>='0')&&(line_tmp[0]<='9')&&((line_tmp[0]!='*')))
			{
				State s_tmp(line);
				s_tmp.UseFlag=1;
				if(s_tmp.Good())
				{	
					if((s_tmp.JP.size()>1)||(s_tmp.L.size()>1)||(s_tmp.n.size()>1)||(s_tmp.SpectroscopicFactor.size()>1))
					{
						IndexesOfMultipleStates.push_back(States.size());
						cout<<"Multiple"<<IndexesOfMultipleStates.size();
					}
					States.push_back(s_tmp);
					States[States.size()-1].type=type;
					States[States.size()-1].JP0=JP0;
				}
			}
			else if((line_tmp[1]>='0')&&(line_tmp[1]<='9')&&((line_tmp[0]=='*')))
			{
				line_tmp.erase(0,1);
				
				State s_tmp(line);
				s_tmp.UseFlag=0;
				cout<<s_tmp.Good()<<"*\n";
				if(s_tmp.Good())
				{	
					States.push_back(s_tmp);
					States[States.size()-1].type=type;
					States[States.size()-1].JP0=JP0;
				}					
			}
			else
			{
				int iterator=0;
				stringstream str_stream(line);
				while(str_stream)
				{
					string tmp_string;
					str_stream>>tmp_string;
					if(tmp_string=="E")
					{
						E_iterator=iterator;
					}
					if(tmp_string=="n")
					{
						n_iterator=iterator;
					}
					if(tmp_string=="l")
					{
						L_iterator=iterator;
					}
					if(tmp_string=="JP")
					{
						JP_iterator=iterator;
					}
					if(tmp_string=="C2S")
					{
						SF_iterator=iterator;
					}
					iterator++;
					if(iterator>4)
					{
						break;
					}
				}
			}
		}
	}
	for(unsigned int i=0;i<States.size();i++)
	{
		if(States[i].Good()==0)
		{
			ErrorString+="Bad state found:"+to_string(i)+"\n";
		}
	}
}

void Experiment::ProcessExperimentalData()//надо будет добавить перебор n и j
{
	cout<<"Exp:"<<reference<<"\n";
	SSD.Calculate(States);
}
void Experiment::Scale(float kNorm,string option)
{
	for(unsigned int i=0;i<States.size();i++)
	{
		States[i].Scale(kNorm,option);
	}
	ProcessExperimentalData();
}
void Experiment::CalculateNPH(parameters par)
{
	vector<StateParameters> states;
	NPH=0;
	for(int i=0;i<SSD.size();i++)
	{
		StateParameters s(SSD.States[i].n,SSD.States[i].L,SSD.States[i].JP,"");
		if(par.CheckSatesForPenaltyFunction(s))
		{
			states.push_back(s);
		}		
	}
	for(unsigned int i=0;i<states.size();i++)
	{
		NPH+=SSD.GetState(states[i].n,states[i].l,states[i].JP).SumG;
	}
}
float Experiment::GetNumberOfHoles()
{
	if(type==1)
	{
		return NPH;
	}
	else
	{
		return -1;
	}
}
float Experiment::GetNumberOfParticles()
{
	if(type==0)
	{
		return NPH;
	}
	else
	{
		return -1;
	}
}
double Experiment::GetCentroid(int n,int l,double JP_inp)//Возвращает центроид для данных JP
{
	return SSD.GetState(n,l,JP_inp).C;
}
double Experiment::GetSumSF(int n,int l,double JP_inp)//Возвращает сумму СФ для данных JP
{
	return SSD.GetState(n,l,JP_inp).SumG;
}
double Experiment::GetSumSF(StateParameters &s)
{
	return SSD.GetState(s.n,s.l,s.JP).SumG;
}
double Experiment::GetCentroid(StateParameters &s)
{
	return SSD.GetState(s.n,s.l,s.JP).C;
}
int Experiment::GetNlevels()//возвращает число уровней, для которых вычислены центроиды
{
	return SSD.size();
}
int Experiment::size()//возвращает число состояний, зарегистрированных в эксперименте
{
	return States.size();
}
int Experiment::SSSsize()
{
	return SSD.size();
}
SummarizedSpectroscopicState&  Experiment::operator [] (int index)
{
	if(index<SSD.States.size())
	{
		return SSD.States[index];
	}
	else
	{
		cout<<reference<<" SummarizedSpectroscopicState index error: index="<<index<<" States.size()="<<SSD.States.size()<<"\n";
		SummarizedSpectroscopicState s1;
		return s1;
	}
}
//vector<TH1F> BuildSpectroscopicFactorHistogram(float &maximum)
SpectroscopicFactorHistogram  Experiment::BuildSpectroscopicFactorHistogram()
{
	SpectroscopicFactorHistogram SFHistograms;
	//Histograms.resize(0);
	////cout<<"SSD.size(): "<<SSD.size()<<"\n";
	SFHistograms.maximum=0;
	SFHistograms.Reference=reference;
	for(unsigned int i=0;i<SSD.size();i++)
	{
		string name=reference+sprintf("_%d_%d_%d",SSD.States[i].n,SSD.States[i].L,(int)(SSD.States[i].JP*2));
		TH1F histogram(name.c_str(),name.c_str(),100,0,10000);
		for(unsigned int j=0;j<States.size();j++)
		{
			if((States[j].n[0]==SSD.States[i].n)&&(States[j].L[0]==SSD.States[i].L)&&(States[j].JP[0]==SSD.States[i].JP))
			histogram.SetBinContent(histogram.GetXaxis()->FindFixBin(States[j].Energy),States[j].G());
		}
		histogram.SetLineColor(GetColor(SSD.States[i].L, SSD.States[i].JP));
		histogram.SetFillColor(GetColor(SSD.States[i].L, SSD.States[i].JP));
		if(histogram.GetMaximum()>SFHistograms.maximum)
		{
			SFHistograms.maximum=histogram.GetMaximum();
		} 
		SFHistograms.Histograms.push_back(histogram);
		SFHistograms.n.push_back(SSD.States[i].n);
		SFHistograms.L.push_back(SSD.States[i].L);
		SFHistograms.JP.push_back(SSD.States[i].JP);
	}
	return SFHistograms;
}	
TH1F CoupleOfExperiments::BuildPenaltyComponentsHistogram()
{
	TH1F result("PFC","Penalty function components",PenaltyComponents.size()+1,0,PenaltyComponents.size()+1);
	result.GetYaxis()->SetRangeUser(10e-12,1);
	for(unsigned int i=0;i<PenaltyComponents.size();i++)
	{
		result.SetBinContent(i+1,PenaltyComponents[i]);
		result.GetXaxis()->SetBinLabel(i+1,par.GetComponentName(i).c_str());
	}
	return result;
}

//	TFormula *PenaltyFormula;

/*	void ConstructFormula()
{
	PenaltyFunction.replace(PenaltyFunction.find("a_{ij}"),6,"[0]");//поиск подстроки 
	PenaltyFunction.replace(PenaltyFunction.find("N^+"),3,"[1]");
	PenaltyFunction.replace(PenaltyFunction.find("N^-"),3,"[2]");
	PenaltyFunction.replace(PenaltyFunction.find("\sigma(E_f)"),11,"[3]");
	PenaltyFunction.replace(PenaltyFunction.find("\sigma(\Delta)"),14,"[4]");
	PenaltyFunction.replace(PenaltyFunction.find("N^+_{max}"),9,"[5]");
	PenaltyFunction.replace(PenaltyFunction.find("N^-_{max}"),9,"[6]");
	PenaltyFunction.replace(PenaltyFunction.find("\sigma(E_f)_{max}"),17,"[7]");
	PenaltyFunction.replace(PenaltyFunction.find("\sigma(\Delta)_{max}"),20,"[8]");
	PenaltyFormula=new TFormula("Penalty",PenaltyFunction.c_str());
	
	PenaltyFormula->SetParameter(0,a);
	PenaltyFormula->SetParameter(1,NumberOfStrippingStates);
	PenaltyFormula->SetParameter(2,NumberOfPickupStates);
	PenaltyFormula->SetParameter(3,Ef_error);
	PenaltyFormula->SetParameter(4,Delta_error);
	PenaltyFormula->SetParameter(5,NumberOfStrippingStatesMax);
	PenaltyFormula->SetParameter(6,NumberOfPickupStatesMax);
	PenaltyFormula->SetParameter(7,Ef_error_max);
	PenaltyFormula->SetParameter(8,Delta_error_max);
	penalty=PenaltyFormula->Eval(0);
}
*/
CoupleOfExperiments::CoupleOfExperiments(Experiment &InpPickup,Experiment &InpStripping)//конструктор, аргументы которого представляют из себя эксперименты по подхвату и срыву
{
	Pickup=InpPickup;
	Stripping=InpStripping;
	FitStatus=0;
	OccStatus=0;
	NumberOfParticlesInUsedShell=0;
	NumberOfHolesInUsedShell=0;
}
void CoupleOfExperiments::GenerateCommonNJPList()
{
	for(int i=0;i<Pickup.SSSsize();i++)
	{
		unsigned char flag=0;
		for(int j=0;j<Stripping.SSSsize();j++)
		{
			if(Pickup[i].Compare(Stripping[j]))
			{
				StateParameters s(Pickup[i].n,Pickup[i].L,Pickup[i].JP,"both");
				SP.push_back(s);
				flag=1;
			}
		}
		if((flag==0)&&((par.IncompleteCouplesFlag==2)||(par.IncompleteCouplesFlag==1)))
		{
			StateParameters s(Pickup[i].n,Pickup[i].L,Pickup[i].JP,"pickup");
			SP.push_back(s);
		}
		
	}//не особо эффективный алгоритм поиска отсутствия совпадений, но сейчас скорость не особо принципиальна
	for(int i=0;i<Stripping.SSSsize();i++)
	{
		unsigned char flag=0;
		for(int j=0;j<Pickup.SSSsize();j++)
		{
			if(Pickup[j].Compare(Stripping[i]))
			{
				flag=1;
			}
		}
		if((flag==0)&&((par.IncompleteCouplesFlag==3)||(par.IncompleteCouplesFlag==1)))
		{
			StateParameters s(Stripping[i].n,Stripping[i].L,Stripping[i].JP,"stripping");
			SP.push_back(s);
		}		
	}
}
void CoupleOfExperiments::CalcSPE_and_OCC(TCanvas *cc1)
{
	//cout<<"\n Generate \n";
	GenerateCommonNJPList();
	//cout<<"\n Generated \n";
	vector<double> OccupanciesForBCSFit;//отдельные векторы заселенностей для аппроксимации БКШ
	vector<double> EnergiesForBCSFit;
	for(int i=0;i<SP.size();i++)
	{
		double C_pickup=Pickup.GetCentroid(SP[i]);
		double C_stripping=Stripping.GetCentroid(SP[i]);	
	
		if((C_stripping!=-1)&&(C_pickup!=-1)&&(!isnan(C_stripping))&&(!isnan(C_pickup)))//индусский fix, потом проверить, что генерирует nan
		{
			double E_pickup=-Pickup.BA-C_pickup;//Диплом Марковой М.Л., ф-ла 4
			double E_stripping=-Stripping.BA1+C_stripping;//Диплом Марковой М.Л., ф-ла 5
			//cout<<"E_pickup="<<E_pickup<<" E_stripping="<<E_stripping<<"\n";
			float G_plus=0,G_minus=0;
			
			cout<<"p:"<<Pickup.reference<<" s:"<<Stripping.reference<<"\n";
			
			G_plus=Stripping.GetSumSF(SP[i]);
			G_minus=Pickup.GetSumSF(SP[i]);
			cout<<"BA:"<<Pickup.BA<<" BA1:"<<Stripping.BA1<<"\n";
			cout<<"G_plus:"<<G_plus<<" G_minus:"<<G_minus<<" E_pickup:"<<E_pickup<<" E_stripping:"<<E_stripping<<"\n";
			double SPE_tmp=(G_minus*E_pickup+G_plus*E_stripping)/(G_minus+G_plus);
			
			cout<<NLJToString(SP[i].n,SP[i].l,SP[i].JP)<<" pickup_c: "<<C_pickup<<" stripping_c:"<<C_stripping<<" "<<Pickup.GetSumSF(SP[i])*E_pickup+Stripping.GetSumSF(SP[i])*E_stripping<<" "<<Pickup.GetSumSF(SP[i])+Stripping.GetSumSF(SP[i])<<"\n";
			cout<<NLJToString(SP[i].n,SP[i].l,SP[i].JP)<<" pickup_sum: "<<Pickup.GetSumSF(SP[i])<<" pickup_E:"<<E_pickup<<" stripping_sum:"<<Stripping.GetSumSF(SP[i])<<" stripping_E: "<<E_stripping<<" E:"<<SPE_tmp<<"\n";
			
			SPE.push_back(SPE_tmp);//Диплом Марковой М.Л., ф-ла 17
			double OCC_tmp=(Pickup.GetSumSF(SP[i])-Stripping.GetSumSF(SP[i])+2*abs(SP[i].JP)+1)/(2*(2*abs(SP[i].JP)+1));
			if(OCC_tmp<0)
			{
				OccStatus=1;
			}
			OCC.push_back(OCC_tmp);//Диплом Марковой М.Л., ф-ла 18
			if(par.CheckSatesForPenaltyFunction(SP[i]))
			{
				ParticlesAndHolesSum.push_back((Pickup.GetSumSF(SP[i])+Stripping.GetSumSF(SP[i]))/(2*abs(SP[i].JP)+1));
				NumberOfParticlesInUsedShell+=Pickup.GetSumSF(SP[i]);
				NumberOfHolesInUsedShell+=Stripping.GetSumSF(SP[i]);
			}
			SP_centroids.push_back(SP[i]);
			
			if(par.CheckOccupancy(SP[i]))
			{
				OccupanciesForBCSFit.push_back(OCC_tmp);//отдельные векторы заселенностей для аппроксимации БКШ
				EnergiesForBCSFit.push_back(SPE_tmp);
			}
		}
		else
		{
			SPE.push_back(0);
			OCC.push_back(0);
		}
	}
	occupancies=TGraph(OccupanciesForBCSFit.size(),&EnergiesForBCSFit[0],&OccupanciesForBCSFit[0]);
	BCS=TF1("BCS","0.5*(1-(x-[0])/(sqrt((x-[0])^2+[1]^2)))",-50000,0);
	BCS.SetParameter(0,-8000);
	BCS.SetParameter(1,15000);
	float min_E,max_E,min_OCC,max_OCC;
	for(unsigned int i=0;i<OCC.size();i++)
	{
		if(SP_centroids[i].GetType()=="pickup")
		{
			Pickup_occupancies.SetPoint(Pickup_occupancies.GetN(),SPE[i],OCC[i]);
		}
		else if(SP_centroids[i].GetType()=="stripping")
		{
			Stripping_occupancies.SetPoint(Stripping_occupancies.GetN(),SPE[i],OCC[i]);
		}
		else if(SP_centroids[i].GetType()=="both")
		{
			Both_occupancies.SetPoint(Both_occupancies.GetN(),SPE[i],OCC[i]);
		}
		if(min_E>SPE[i])
		{
			min_E=SPE[i];
		}
		if(max_E<SPE[i])
		{
			max_E=SPE[i];
		}
		if(min_OCC>OCC[i])
		{
			min_OCC=OCC[i];
		}
		if(max_OCC<OCC[i])
		{
			max_OCC=OCC[i];
		}
	}
	
	SetTGraphLimits(Pickup_occupancies,min_E,max_E,min_OCC,max_OCC);
	Pickup_occupancies.SetMarkerColor(4);
	SetTGraphLimits(Stripping_occupancies,min_E,max_E,min_OCC,max_OCC);
	Stripping_occupancies.SetMarkerColor(2);
	SetTGraphLimits(Both_occupancies,min_E,max_E,min_OCC,max_OCC);
	Both_occupancies.SetMarkerColor(1);
	occupancies.Draw("AP");	
	occupancies.SetMarkerStyle(28);
	occupancies.SetMarkerSize(2);
	FitStatus=occupancies.Fit(&BCS,"M");
	BCS.Draw("l same");
	
	cc1->Print("BCS.pdf","pdf");
	Ef=BCS.GetParameter(0);
	Ef_error=BCS.GetParError(0);
	Delta=BCS.GetParameter(1);
	Delta_error=BCS.GetParError(1);
}
string CoupleOfExperiments::ResultsInTextForm(char verbose_level)
{
	stringstream s;
	if(verbose_level==0)
	{
		s<<"Experiment: "<<Pickup.reference<<" ("<<Pickup.size()<<") "<<Stripping.reference<<" ("<<Stripping.size()<<") \n";
	}
	else if(verbose_level==1)
	{
		s<<"Experiment(pickup): "<<Pickup.reference<<" ("<<Pickup.size()<<")\n";
		s<<Pickup.ChangesLog<<"\n";
		s<<"Experiment(stripping): "<<Stripping.reference<<" ("<<Stripping.size()<<")\n";
		s<<Stripping.ChangesLog<<"\n";
	}
	
	s<<Pickup.particle<<" transfer\n";
	s<<Pickup.particle[0]<<" separation energy A:"<<Pickup.BA<<", A+1: "<<Pickup.BA1<<"\n";
	s<<"E_F: "<<Ef<<" #pm "<<Ef_error<<"  keV \n #Delta: "<<Delta<<" #pm "<<Delta_error<<" keV\n";
	s<<"penalty: "<<penalty<<"\n";
	s<<"SPE,keV nlj OCC #frac{G^{+}+G^{-}}{2J+1}\n";
	for(unsigned int i=0;i<SPE.size();i++)
	{
		s<<SPE[i]<<" "<<NLJToString(SP_centroids[i].n,SP_centroids[i].l,SP_centroids[i].JP)<<" "<<OCC[i]<<" "<<ParticlesAndHolesSum[i]<<"\n";
	}
	return s.str();
}
void CoupleOfExperiments::DrawResultsInTextForm()
{
	stringstream s(ResultsInTextForm());
	string LatexLineTmp;
	TLatex latex;
	float x=0.05,y=0.9, step=-0.08;
	while(getline(s,LatexLineTmp))
	{
		latex.DrawLatexNDC(x,y,LatexLineTmp.c_str());
		y+=step;
	}
}
int CoupleOfExperiments::GetNumberOfPickupStates()
{
	int value=0;
	for(unsigned int i=0;i<Pickup.States.size();i++)
	{
		if(par.CheckSatesForPenaltyFunction(Pickup.States[i]))
		{
			value++;
		}		
	}
	return value;
}
int CoupleOfExperiments::GetNumberOfStrippingStates()
{
	int value=0;
	for(unsigned int i=0;i<Stripping.States.size();i++)
	{
		if(par.CheckSatesForPenaltyFunction(Stripping.States[i]))
		{
			value++;
		}		
	}
	return value;
}
