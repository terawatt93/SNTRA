#include "State.hh"

double State::G(int JP_Index,int SF_Index)
{
	if(type==0)//pickup
	{
		return SpectroscopicFactor[SF_Index];
	}
	else if(type==1)//stripping
	{
		if(JP_Index<JP.size())
		{
			return (2*abs(JP[JP_Index])+1)/(2*abs(JP0)+1)*SpectroscopicFactor[SF_Index];
		}
		else
		{
			//cout<<"Error in function G(int JP_Index): wrong JP_Index="<<JP_Index<<" when JP.size()="<<JP.size()<<"\n";
		}
	}
	return 0;
}
double State::GetJP(int JP_Index)
{
	////cout<<"JP_Index="<<JP_Index<<"\n";
	if(JP.size()>JP_Index)
	{
		return JP[JP_Index];
	}
	else
	{
		return -1;
	}
}
double State::GetN(int N_Index)
{
	if(n.size()>N_Index)
	{
		return n[N_Index];
	}
	else
	{
		return -1;
	}
}
int State::GetL(int L_index)
{
	return L[L_index];
}
State::State()
{
	JP0=0;
	UseFlag=0;
}
State::State(string InputString,int E_iterator, int n_iterator, int L_iterator, int JP_iterator, int SF_iterator)
{
	////cout<<"iterators:"<<E_iterator<<" "<<L_iterator<<" "<<SF_iterator<<" "<<n_iterator<<" "<<JP_iterator<<"\n";
	JP0=0;
	stringstream str_stream(InputString);		
	if(InputString.size()>0)
	{
		int iterator=0;
		while(str_stream)
		{
			////cout<<iterator<<"\n";
			if(iterator==E_iterator)
			{
				str_stream>>Energy;
				//cout<<"E="<<Energy<<" ";
				iterator++;
			}
			if(iterator==L_iterator)
			{
				string buffer;
				string buffer2;
				str_stream>>buffer;
				for(int i=0;i<buffer.size();i++)
				{
					if(buffer[i]!=';')
					{
						buffer2+=buffer[i];
					}
					else
					{
						L.push_back(atof(buffer2.c_str()));
						buffer2.resize(0);					
					}
				}
				//if(L.size()==0)
				//{
				//	//cout<<"bufn="<<buffer2<<"\n";
					L.push_back(atof(buffer2.c_str()));
				//}
				iterator++;
			}
			if(iterator==SF_iterator)
			{
				string buffer;
				string buffer2;
				str_stream>>buffer;
				for(int i=0;i<buffer.size();i++)
				{
					if(buffer[i]!=';')
					{
						buffer2+=buffer[i];
					}
					else
					{
						SpectroscopicFactor.push_back(atof(buffer2.c_str()));
						buffer2.resize(0);						
					}
				}
				//if(SpectroscopicFactor.size()==0)
				//{
				//	//cout<<"bufn="<<buffer2<<"\n";
					SpectroscopicFactor.push_back(atof(buffer2.c_str()));
				//}
				iterator++;
			}
			if(iterator==n_iterator)
			{
				string buffer;
				string buffer2;
				str_stream>>buffer;
				for(int i=0;i<buffer.size();i++)
				{
					if(buffer[i]!=';')
					{
						buffer2+=buffer[i];
					}
					else
					{
						n.push_back(atoi(buffer2.c_str()));
						buffer2.resize(0);					
					}
				}
				//if(n.size()==0)
				//{
				//	//cout<<"bufn="<<buffer2<<"\n";
					n.push_back(atof(buffer2.c_str()));
				//}
				
				iterator++;
			}
			if(iterator==JP_iterator)
			{
				string buffer;
				string buffer2;
				str_stream>>buffer;
				for(int i=0;i<buffer.size();i++)
				{
					//cout<<buffer[i];
					if(buffer[i]!=';')
					{
						buffer2+=buffer[i];
						//cout<<buffer[i];
					}
					else
					{
						cout<<"!";
						JP.push_back(atof(buffer2.c_str()));
						buffer2.resize(0);					
					}
				}
				//if(JP.size()==0)
				//{
					////cout<<"bufJP="<<buffer2<<"\n";
					JP.push_back(atof(buffer2.c_str()));
				//}
				iterator++;
				
			}
			if(iterator>4)
			{
				break;
			}	
		}
	}
}
int State::Good()
{
	if((JP.size()>0)&&(n.size()>0))
	{
		for(unsigned int i=0;i<n.size();i++)
		{
			for(unsigned int j=0;j<L.size();j++)
			{
				for(unsigned int k=0;k<JP.size();k++)
				{
					if(abs(L[j]-JP[k])>0.5)
					{
						return 0;
					}
				}
			}
		}		
		return 1;
	}
	else
	{
		return 0;
	}
}
void State::Scale(float kNorm, string option)
{
	if(SpectroscopicFactor.size()>0)
	{
		SpectroscopicFactor[0]=SpectroscopicFactor[0]*kNorm;
	}
}
SummarizedSpectroscopicState::SummarizedSpectroscopicState()
{
	JP=0; C=0; SumG=0; n=0;
}
SummarizedSpectroscopicState::SummarizedSpectroscopicState(int n_value,int L_value,double JP_value)
{
	JP=JP_value; C=0; SumG=0; n=n_value; L=L_value;
}
SummarizedSpectroscopicState::SummarizedSpectroscopicState(int n_value,int L_value,double JP_value, double C_value,double SumG_value)
{
	JP=JP_value; C=C_value; SumG=SumG_value; n=n_value; L=L_value;
}
const bool SummarizedSpectroscopicState::Compare(SummarizedSpectroscopicState &S)
{
	////cout<<JP<<" "<<S.JP<<" "<<n<<" "<<S.n<<"\n";
	if(SumG<0)
	{
		return 0;
	}
	if((JP==S.JP)&&(n==S.n)&&(L==S.L))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
const bool SummarizedSpectroscopicState::Compare(State &S, int N_index, int JP_index)
{
	if((JP==S.JP[JP_index])&&(n==S.n[N_index])&&(L==S.L[0])&&(S.UseFlag==1))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
const bool SummarizedSpectroscopicState::Compare(int n_value, int L_value, double JP_value)
{
	////cout<<JP<<" "<<JP_value<<" "<<n<<" "<<n_value<<"\n";
	
	if((JP==JP_value)&&(n==n_value)&&(L==L_value))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
const bool SummarizedSpectroscopicState::operator == (SummarizedSpectroscopicState &v2)
{
	if((SumG<0)||(v2.SumG<0))
	{
		return 0;
	}
	if((n==v2.n)&&(JP==v2.JP)&&(L==v2.L))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
void SummarizedSpectroscopicData::CreateJPList(State &s)
{
	char exsistance_flag=0;
	for(int i=0;i<States.size();i++)
	{
		if(States[i].Compare(s)==1)
		{
			exsistance_flag=1;
		}
	}
	if((exsistance_flag==0)&&(s.UseFlag==1))
	{
		SummarizedSpectroscopicState state_tmp(s.GetN(),s.GetL(),s.GetJP());
		////cout<<"CreateJPList"<<s.GetN()<<" "<<s.GetL()<<" "<<s.GetJP()<<"\n";
		States.push_back(state_tmp);
	}
}
void SummarizedSpectroscopicData::Calculate(vector<State> &states)
{
	States.resize(0);
	for(int i=0;i<states.size();i++)
	{
		CreateJPList(states[i]);
	}
	for(int i=0;i<states.size();i++)//вычисление сумм СФ и центроидов
	{
		for(int j=0;j<States.size();j++)
		{
			if(States[j].Compare(states[i])==1)
			{
				States[j].SumG+=states[i].G();
				cout<<states[i].GetJP()<<" G:"<<states[i].G()<<"\n";
				States[j].C+=states[i].Energy*states[i].G();
				cout<<NLJToString(states[i].n[0],states[i].L[0],states[i].JP[0])<<" "<<states[i].Energy<<" "<<states[i].G()<<"\n";
			}
		}
	}
	for(int i=0;i<States.size();i++)
	{
		States[i].C=States[i].C/States[i].SumG;
		cout<<"centroid:"<<NLJToString(States[i].n,States[i].L,States[i].JP)<<" "<<States[i].C<<" G:"<<States[i].SumG<<"\n";
	}
}
SummarizedSpectroscopicState SummarizedSpectroscopicData::GetState(int n, int l, double JP)
{
	////cout<<"GetState:"<<n<<" "<<JP<<"\n";
	for(int i=0;i<States.size();i++)
	{
		////cout<<"States[i].Compare(n,JP):"<<States[i].Compare(n,JP)<<"\n";
		if(States[i].Compare(n,l,JP)==1)
		{
			return States[i];
		}
	}
	SummarizedSpectroscopicState s;
	return s;
	
}
int SummarizedSpectroscopicData::size()
{
	return States.size();
}
