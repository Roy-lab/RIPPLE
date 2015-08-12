#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "Dataset.H"
#include "Distance.H"
#include "Framework.H"

Framework::Framework()
{
	filterZeros=true;
}

Framework::~Framework()
{
}

//Read the enhancer-promoter interactions
int 
Framework::readPairs(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strchr(buffer,'-')==NULL && strchr(buffer,'\t')==NULL)
		{
			continue;
		}
		char* tok=strchr(buffer,'\t');
		if(tok!=NULL)
		{
			*tok='\0';
		}
		char* e=buffer;
		char* p=tok+1;
		char* tok2=strchr(p,'\t');
		if(tok2!=NULL)
		{
			*tok2='\0';
		}
		char* primerpref=strrchr(e,'|');
		if(primerpref!=NULL)
		{
			e=primerpref+1;
		}
		
		int i=0;
		while(e[i]!='\0')
		{
			if(e[i]=='-' || e[i]==':' )
			{
				e[i]='_';
			}
			i++;
		}
		primerpref=strrchr(p,'|');
		if(primerpref!=NULL)
		{
			p=primerpref+1;
		}
		i=0;
		while(p[i]!='\0')
		{
			if(p[i]=='-' || p[i]==':')
			{
				p[i]='_';
			}
			i++;
		}
		string eKey(e);
		if(regionSet.find(eKey)==regionSet.end())
		{
			Framework::Region* eregion=new Framework::Region;
			regionSet[eKey]=eregion;	
			string chrom;	
			int start;
			int end;
			tok=strtok(e,"_");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					eregion->chromosome.append(tok);
				}
				else if(tokCnt==1)
				{	
					eregion->begin=atoi(tok);
				}
				else if(tokCnt==2)
				{
					eregion->end=atoi(tok);
				}
				tok=strtok(NULL,"_");
				tokCnt++;
			}
		}
		string pKey(p);
		if(regionSet.find(pKey)==regionSet.end())
		{
			Framework::Region* pregion=new Framework::Region;
			regionSet[pKey]=pregion;	
			tok=strtok(p,"_");
			int tokCnt=0;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					pregion->chromosome.append(tok);
				}
				else if(tokCnt==1)
				{	
					pregion->begin=atoi(tok);
				}
				else if(tokCnt==2)
				{
					pregion->end=atoi(tok);
				}
				tok=strtok(NULL,"_");
				tokCnt++;
			}
			regionSet[pKey]=pregion;
		}
		map<string,int>* promSet=NULL;
		if(pairSet.find(eKey)==pairSet.end())
		{
			promSet=new map<string,int>;
			pairSet[eKey]=promSet;
		}	
		else
		{
			promSet=pairSet[eKey];
		}
		(*promSet)[pKey]=0;
	}
	inFile.close();
	return 0;
}

int
Framework::setFilterZeros(bool flag)
{
	filterZeros=flag;
	return 0;
}

int
Framework::setCorrelation(bool flag)
{
	correlation=flag;
	return 0;
}

int 
Framework::readFeatures(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int fID=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}	
		if(strstr(buffer,"#")!=NULL)
		{
			continue;
		}
		string featName;
		char fileName[1024];
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				featName.append(tok);
			}
			else if(tokCnt==1)
			{
				strcpy(fileName,tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		cout <<"Reading " << fID <<" "<< featName <<" ";
		Dataset* d=new Dataset;
		if(datasets.find(featName)==datasets.end())
		{
			//featureOrder.push_back(featName);
			datasets[featName]=d;
			d->readDataSet(fileName);
		}
		else
		{
			d=datasets[featName];
			d->readDataSet(fileName);
			
		}
		fID++;	
	}
	inFile.close();
	for(map<string,Dataset*>::iterator dIter=datasets.begin();dIter!=datasets.end();dIter++)
	{
		featureOrder.push_back(dIter->first);
	}
	return 0;
}

int
Framework::readExpression(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string region;
		double exp;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				region.append(tok);
			}
			else
			{
				exp=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(geneexp.find(region)==geneexp.end())
		{
			geneexp[region]=exp;
		}
		else
		{
			double cval=geneexp[region];
			if(cval<exp)
			{
				geneexp[region]=exp;
			}
		}
	}
	inFile.close();
	return 0;
}

double
Framework::getExp(string& region)
{
	double exp=0;
	if(geneexp.find(region)!=geneexp.end())
	{
		exp=geneexp[region];
	}
	return exp;
}
int 
Framework::generateFeatureFiles(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f];
	}
	//oFile <<"\tDist\tCorrelation";
	if(correlation)
	{
		oFile <<"\tCorrelation";
	}
	if(geneexp.size()>0)
	{
		oFile <<"\tExp";
	}
	oFile << "\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	for(int f=0;f<featureOrder.size();f++)
	{
		eFile <<"\t" << featureOrder[f];
	}
	eFile << endl;
	ofstream pFile(aFName_P);
	for(int f=0;f<featureOrder.size();f++)
	{
		pFile <<"\t" << featureOrder[f];
	}
	pFile << endl;

	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0)
			{
				cout <<"Stop here " <<endl;
			}
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				if(strcmp(p.c_str(),"chr7_116434729_116454408")==0)
				{
					cout <<"Stop here " <<endl;
				}
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0 && strcmp(p.c_str(),"chr7_116434729_116454408")==0)
			{
				//showFeatures(p,pfeatures);
				cout <<"Stop here"<< endl;
			}
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_e=0.01;
				double featVal_p=0.01;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=(*efeatures)[featureOrder[f]]+0.01;
					founde=true;
				}
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				if(founde ||foundp)
				{
					if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
					{
						featureFrequency[featureOrder[f]]=1;
					}
					else
					{
						featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
					}
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{	
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				oFile <<"\t" <<featuresForPair[f];
			}
			int dist=getDistance(e,p);
			if(correlation)
			{
				double cc=getCCFeature(pfeatures,efeatures);
				oFile <<"\t"  <<cc; 
			}
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"<<exp;
			}
			oFile << "\t" << classLabel<< endl;
			featuresForPair.clear();
		}
	}
	oFile.close();
	//Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
	//although the user can decide
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
		
	return 0;
}


int 
Framework::generateFeatureFiles_ConcatProd(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f];
	}
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_E";
	}
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_P";
	}
	//oFile <<"\tDist\tCorrelation";
	oFile <<"\tCorrelation";
	oFile << "\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	for(int f=0;f<featureOrder.size();f++)
	{
		eFile <<"\t" << featureOrder[f];
	}
	eFile << endl;
	ofstream pFile(aFName_P);
	for(int f=0;f<featureOrder.size();f++)
	{
		pFile <<"\t" << featureOrder[f];
	}
	pFile << endl;
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0)
			{
				cout <<"Stop here " <<endl;
			}
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				if(strcmp(p.c_str(),"chr7_116434729_116454408")==0)
				{
					cout <<"Stop here " <<endl;
				}
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0 && strcmp(p.c_str(),"chr7_116434729_116454408")==0)
			{
				//showFeatures(p,pfeatures);
				cout <<"Stop here"<< endl;
			}
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_e=0.01;
				double featVal_p=0.01;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=(*efeatures)[featureOrder[f]]+0.01;
					founde=true;
				}
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				if(founde ||foundp)
				{
					if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
					{
						featureFrequency[featureOrder[f]]=1;
					}
					else
					{
						featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
					}
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{	
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				oFile <<"\t" <<featuresForPair[f];
			}
			//Now show the enhancer features
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_e=0.01;
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=(*efeatures)[featureOrder[f]]+0.01;
				}
				oFile <<"\t" << featVal_e;
			}
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_p=0.01;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
				}
				oFile <<"\t" << featVal_p;
			}
			int dist=getDistance(e,p);
			double cc=getCCFeature(pfeatures,efeatures);
			//oFile <<"\t" << dist<< "\t" <<cc << "\t" << classLabel<< endl;
			oFile <<"\t" <<cc << "\t" << classLabel<< endl;


			//Now show the promoter features
			//oFile << "\t" << classLabel<< endl;
			featuresForPair.clear();
		}
	}
	oFile.close();
	//Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
	//although the user can decide
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
		
	return 0;
}


int 
Framework::generateFeatureFiles_Concat(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_E";
	}
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_P";
	}
	//oFile <<"\tDist\tCorrelation";
	if(correlation)
	{
		oFile <<"\tCorrelation";
	}
	if(geneexp.size()>0)
	{
		oFile <<"\tExp";
	}
	oFile << "\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	for(int f=0;f<featureOrder.size();f++)
	{
		eFile <<"\t" << featureOrder[f];
	}
	eFile << endl;
	ofstream pFile(aFName_P);
	for(int f=0;f<featureOrder.size();f++)
	{
		pFile <<"\t" << featureOrder[f];
	}
	pFile << endl;
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0)
			{
				cout <<"Stop here " <<endl;
			}
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				if(strcmp(p.c_str(),"chr7_116434729_116454408")==0)
				{
					cout <<"Stop here " <<endl;
				}
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0 && strcmp(p.c_str(),"chr7_116434729_116454408")==0)
			{
				//showFeatures(p,pfeatures);
				cout <<"Stop here"<< endl;
			}
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_e=0.01;
				double featVal_p=0.01;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=(*efeatures)[featureOrder[f]]+0.01;
					founde=true;
				}
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				if(founde ||foundp)
				{
					if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
					{
						featureFrequency[featureOrder[f]]=1;
					}
					else
					{
						featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
					}
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{	
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			oFile <<e<<"-"<<p;
			//Now show the enhancer features
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_e=0.01;
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=(*efeatures)[featureOrder[f]]+0.01;
				}
				oFile <<"\t" << featVal_e;
			}
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_p=0.01;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
				}
				oFile <<"\t" << featVal_p;
			}

			int dist=getDistance(e,p);
			if(correlation)
			{
				double cc=getCCFeature(pfeatures,efeatures);
				oFile <<"\t"  <<cc; 
			}
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"<<exp;
			}
			//oFile <<"\t" << dist<< "\t" <<cc << "\t" << classLabel<< endl;
			//oFile <<"\t" <<cc << "\t" << classLabel<< endl;
			oFile << "\t" << classLabel<< endl;
			featuresForPair.clear();

		}
	}
	oFile.close();
	//Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
	//although the user can decide
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
		
	return 0;

}

int 
Framework::generateFeatureFiles_Binary(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f];
	}
	//oFile <<"\tDist\tCorrelation";
	if(correlation)
	{
		oFile <<"\tCorrelation";
	}
	if(geneexp.size()>0)
	{
		cout <<"Showing expression! "<< endl;
		oFile <<"\tExp";
	}
	oFile << "\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	for(int f=0;f<featureOrder.size();f++)
	{
		eFile <<"\t" << featureOrder[f];
	}
	eFile << endl;
	ofstream pFile(aFName_P);
	for(int f=0;f<featureOrder.size();f++)
	{
		pFile <<"\t" << featureOrder[f];
	}
	pFile << endl;
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0)
			{
				cout <<"Stop here " <<endl;
			}
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				if(strcmp(p.c_str(),"chr7_116434729_116454408")==0)
				{
					cout <<"Stop here " <<endl;
				}
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0 && strcmp(p.c_str(),"chr7_116434729_116454408")==0)
			{
				//showFeatures(p,pfeatures);
				cout <<"Stop here"<< endl;
			}
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_e=0;
				int featVal_p=0;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=1;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=1;
					founde=true;
				}
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				if(founde ||foundp)
				{
					if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
					{
						featureFrequency[featureOrder[f]]=1;
					}
					else
					{
						featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
					}
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			//Filter zeros only if explicitly asked by user. 
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{		
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				oFile <<"\t" <<featuresForPair[f];
			}
			int dist=getDistance(e,p);
			if(correlation)
			{
				double cc=getCCFeature(pfeatures,efeatures);
				//double cc=getCCFeature_Bin(pfeatures,efeatures);
				oFile <<"\t"  <<cc; 
			}
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"<<exp;
			}
			oFile << "\t" << classLabel<< endl;
			featuresForPair.clear();
		}
	}
	oFile.close();
	//Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
	//although the user can decide
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
	return 0;
}


int 
Framework::generateFeatureFiles_ConcatProd_Binary(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f];
	}
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_E";
	}
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_P";
	}
	//oFile <<"\tDist\tCorrelation";
	oFile <<"\tCorrelation";
	oFile << "\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	for(int f=0;f<featureOrder.size();f++)
	{
		eFile <<"\t" << featureOrder[f];
	}
	eFile << endl;
	ofstream pFile(aFName_P);
	for(int f=0;f<featureOrder.size();f++)
	{
		pFile <<"\t" << featureOrder[f];
	}
	pFile << endl;
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0)
			{
				cout <<"Stop here " <<endl;
			}
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				if(strcmp(p.c_str(),"chr7_116434729_116454408")==0)
				{
					cout <<"Stop here " <<endl;
				}
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0 && strcmp(p.c_str(),"chr7_116434729_116454408")==0)
			{
				//showFeatures(p,pfeatures);
				cout <<"Stop here"<< endl;
			}
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_e=0;
				int featVal_p=0;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=1;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=1;
					founde=true;
				}
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				if(founde ||foundp)
				{
					if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
					{
						featureFrequency[featureOrder[f]]=1;
					}
					else
					{
						featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
					}
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			//Filter zeros only if explicitly asked by user. 
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{		
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				oFile <<"\t" <<featuresForPair[f];
			}
			//Now show the enhancer features
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_e=0;
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=1;
				}
				oFile <<"\t" << featVal_e;
			}
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_p=0;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=1;
				}
				oFile <<"\t" << featVal_p;
			}
			int dist=getDistance(e,p);
			double cc=getCCFeature(pfeatures,efeatures);
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"  <<cc << "\t" <<exp <<"\t" << classLabel<< endl;
			}
			else
			{
				oFile <<"\t"  <<cc << "\t" << classLabel<< endl;
			}
			//oFile << "\t" << classLabel<< endl;
			featuresForPair.clear();
		}
	}
	oFile.close();
	//Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
	//although the user can decide
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
	return 0;
}



int 
Framework::generateFeatureFiles_Concat_Binary(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_E";
	}
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f]<<"_P";
	}
	//oFile <<"\tDist\tCorrelation";
	if(correlation)
	{
		oFile <<"\tCorrelation";
	}
	if(geneexp.size()>0)
	{
		cout <<"Showing expression! "<< endl;
		oFile <<"\tExp";
	}
	oFile << "\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	for(int f=0;f<featureOrder.size();f++)
	{
		eFile <<"\t" << featureOrder[f];
	}
	eFile << endl;
	ofstream pFile(aFName_P);
	for(int f=0;f<featureOrder.size();f++)
	{
		pFile <<"\t" << featureOrder[f];
	}
	pFile << endl;
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if((strcmp(e.c_str(),"chr1_1209695_1210344")==0) && (strcmp(p.c_str(),"chr1_1208103_1213103")==0))
			{
				cout <<"Stop here " <<endl;
			}
		
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				if(strcmp(p.c_str(),"chr7_116434729_116454408")==0)
				{
					cout <<"Stop here " <<endl;
				}
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			if(strcmp(e.c_str(),"chr7_116232703_116234744")==0 && strcmp(p.c_str(),"chr7_116434729_116454408")==0)
			{
				//showFeatures(p,pfeatures);
				cout <<"Stop here"<< endl;
			}
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_e=0;
				int featVal_p=0;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=1;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=1;
					founde=true;
				}
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				if(founde ||foundp)
				{
					if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
					{
						featureFrequency[featureOrder[f]]=1;
					}
					else
					{
						featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
					}
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			//Filter zeros only if explicitly asked by user. 
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{		
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			oFile <<e<<"-"<<p;
			//Now show the enhancer features
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_e=0;
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=1;
				}
				oFile <<"\t" << featVal_e;
			}
			for(int f=0;f<featureOrder.size();f++)
			{
				int featVal_p=0;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=1;
				}
				oFile <<"\t" << featVal_p;
			}
			featuresForPair.clear();
			int dist=getDistance(e,p);
			if(correlation)
			{
				double cc=getCCFeature(pfeatures,efeatures);
				//double cc=getCCFeature_Bin(pfeatures,efeatures);
				oFile <<"\t"  <<cc; 
			}
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"<<exp;
			}
			//oFile <<"\t" << dist<< "\t" <<cc << "\t" << classLabel<< endl;
			//oFile << "\t" <<cc << "\t" << classLabel<< endl;
			oFile << "\t" << classLabel<< endl;
		}
	}
	oFile.close();
	//Do the feature frequency check. If there are just zeros, get rid of them and reobtain features. Also, this makes sense for the positive class
	//although the user can decide
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
	return 0;
}

double
Framework::getCCFeature(map<string,double>* pfeatures,map<string,double>* efeatures)
{	
	double cc=0;
	vector<double> pvect;
	vector<double> evect;
	Distance dist;
	for(int i=0;i<featureOrder.size();i++)
	{
		if(pfeatures->find(featureOrder[i])==pfeatures->end())
		{
			pvect.push_back(0.001);
		}
		else
		{
			pvect.push_back((*pfeatures)[featureOrder[i]]);
		}
		if(efeatures->find(featureOrder[i])==efeatures->end())
		{
			evect.push_back(0.001);
		}
		else
		{
			evect.push_back((*efeatures)[featureOrder[i]]);
		}
	}
	cc=dist.computeCC(pvect,evect);
	pvect.clear();
	evect.clear();
	return cc;
}


double
Framework::getCCFeature_Bin(map<string,double>* pfeatures,map<string,double>* efeatures)
{	
	double cc=0;
	double hit=0;
	Distance dist;
	for(map<string,double>::iterator pIter=pfeatures->begin();pIter!=pfeatures->end();pIter++)
	{
		if(efeatures->find(pIter->first)!=efeatures->end())
		{
			hit=hit+1;
		}
	}
	/*for(map<string,double>::iterator pIter=efeatures->begin();pIter!=efeatures->end();pIter++)
	{
		if(pfeatures->find(pIter->first)==pfeatures->end())
		{
			cc=cc+1;
		}
	}*/
	cc=hit/(((double)(pfeatures->size()+efeatures->size()))-hit);
	//cc=cc/((double)featureOrder.size());
	//cc=1-cc;
	return cc;
}


int 
Framework::generateFeatureFilesDist(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		oFile <<"\t" << featureOrder[f];
	}
	oFile << "\tDist\tCorrelation\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	ofstream pFile(aFName_P);
	Distance corr;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			/*if(strcmp(p.c_str(),"chr11_116737744_116757175")==0)
			{
				showFeatures(p,pfeatures);
			}*/
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			vector<double> p_vect;
			vector<double> e_vect;
			for(int f=0;f<featureOrder.size();f++)
			{
				double featVal_e=0.01;
				double featVal_p=0.01;
				bool founde=false;
				bool foundp=false;
				if(pfeatures->find(featureOrder[f])!=pfeatures->end())
				{
					featVal_p=(*pfeatures)[featureOrder[f]]+0.01;
					foundp=true;
				}
				if(efeatures->find(featureOrder[f])!=efeatures->end())
				{
					featVal_e=(*efeatures)[featureOrder[f]]+0.01;
					founde=true;
				}
				p_vect.push_back(featVal_p);
				e_vect.push_back(featVal_e);
				if(!founde)
				{
					miss_e++;
				}
				if(!foundp)	
				{
					miss_p++;
				}
				double f_ep=featVal_e*featVal_p;
				featuresForPair[f]=f_ep;
			}	
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{
					p_vect.clear();
					e_vect.clear();
					continue;
				}
			}
			int dist=getDistance(e,p);
			if(dist>1000000)
			{
				p_vect.clear();
				e_vect.clear();
				continue;
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				oFile <<"\t" <<featuresForPair[f];
			}
			double cc=corr.computeCC(p_vect,e_vect);
			oFile <<"\t" << dist<< "\t" <<cc << "\t" << classLabel<< endl;
			featuresForPair.clear();
			p_vect.clear();
			e_vect.clear();
		}
	}	
	oFile.close();
	return 0;
}


int 
Framework::generateFeatureFilesOuterprod(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		for(int g=0;g<featureOrder.size();g++)
		{
			oFile <<"\t" << featureOrder[f]<<"_E_"<<featureOrder[g]<<"_P";
		}
	}
	if(correlation)
	{
		oFile <<"\tCorrelation";
	}
	if(geneexp.size()>0)
	{
		cout <<"Showing expression! "<< endl;
		oFile <<"\tExp";
	}
	oFile << "\tClass" <<endl;
	//oFile << "\tDist\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	ofstream pFile(aFName_P);
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			/*if(strcmp(p.c_str(),"chr11_116737744_116757175")==0)
			{
				showFeatures(p,pfeatures);
			}*/
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				for(int g=0;g<featureOrder.size();g++)
				{
					double featVal_e=0.01;
					double featVal_p=0.01;
					bool founde=false;
					bool foundp=false;
					if(pfeatures->find(featureOrder[g])!=pfeatures->end())
					{
						featVal_p=(*pfeatures)[featureOrder[g]]+0.01;
						foundp=true;
					}
					if(efeatures->find(featureOrder[f])!=efeatures->end())
					{
						featVal_e=(*efeatures)[featureOrder[f]]+0.01;
						founde=true;
					}
					if(!founde)
					{
						miss_e++;
					}
					if(!foundp)	
					{
						miss_p++;
					}
					double f_ep=featVal_e*featVal_p;
					int pos=(f*featureOrder.size())+g;
					featuresForPair[pos]=f_ep;
					if(founde ||foundp)
					{
						if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
						{
							featureFrequency[featureOrder[f]]=1;
						}
						else
						{
							featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
						}
					}
				}
			}	
			//Filter zeros only if explicitly asked by user. 
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{		
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			int dist=getDistance(e,p);
			if(dist>1000000)
			{
				continue;
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				for(int g=0;g<featureOrder.size();g++)
				{
					int pos=(f*featureOrder.size())+g;
					oFile <<"\t" <<featuresForPair[pos];
				}
			}
			if(correlation)
			{
				double cc=getCCFeature(pfeatures,efeatures);
				oFile <<"\t"  <<cc; 
			}
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"<<exp;
			}
			//oFile <<"\t" << dist<< "\t" << classLabel<< endl;
			oFile <<"\t" << classLabel<< endl;
			featuresForPair.clear();
		}
	}	
	oFile.close();
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
	return 0;

}


int 
Framework::generateFeatureFilesOuterprod_Binary(const char* aFName, const char* classLabel)
{
	ofstream oFile(aFName);
	for(int f=0;f<featureOrder.size();f++)
	{
		for(int g=0;g<featureOrder.size();g++)
		{
			oFile <<"\t" << featureOrder[f]<<"_E_"<<featureOrder[g]<<"_P";
		}
	}
	if(correlation)
	{
		oFile <<"\tCorrelation";
	}
	if(geneexp.size()>0)
	{
		cout <<"Showing expression! "<< endl;
		oFile <<"\tExp";
	}
	oFile << "\tClass" <<endl;
	//oFile << "\tDist\tClass" <<endl;
	char aFName_P[1024];
	char aFName_E[1024];
	sprintf(aFName_E,"%s.enhancer",aFName);
	sprintf(aFName_P,"%s.promoter",aFName);
	ofstream eFile(aFName_E);
	ofstream pFile(aFName_P);
	map<string,int> featureFrequency;
	for(map<string,map<string,int>*>::iterator eIter=pairSet.begin();eIter!=pairSet.end();eIter++)
	{
		string& e=(string&)eIter->first;
		map<string,double>* efeatures=NULL;
		if(regionFeatures.find(e)==regionFeatures.end())
		{
			efeatures=new map<string,double>;
			regionFeatures[e]=efeatures;
			populateFeaturesForRegion(efeatures,e);
			/*eFile << e;
			for(int i=0;i<featureOrder.size();i++)
			{
				if(efeatures->find(featureOrder[i])==efeatures->end())
				{
					eFile <<"\t-1";
				}
				else
				{
					eFile <<"\t"<<(*efeatures)[featureOrder[i]];
				}
				eFile << endl;
			}*/
			showFeatures(e,efeatures,eFile);
		}
		else
		{
			efeatures=regionFeatures[e];
		}
		map<string,int>* pSet=eIter->second;
		/*if(strcmp(e.c_str(),"chr11_116617499_116620420")==0)
		{
			showFeatures(e,efeatures);
		}*/
		for(map<string,int>::iterator pIter=pSet->begin();pIter!=pSet->end();pIter++)
		{
			string& p=(string&)pIter->first;
			map<string,double>* pfeatures=NULL;
			if(regionFeatures.find(p)==regionFeatures.end())
			{
				pfeatures=new map<string,double>;
				regionFeatures[p]=pfeatures;
				populateFeaturesForRegion(pfeatures,p);
				/*pFile << p;
				for(int i=0;i<featureOrder.size();i++)
				{
					if(efeatures->find(featureOrder[i])==efeatures->end())
					{
						eFile <<"\t-1";
					}
					else
					{
						eFile <<"\t"<<(*efeatures)[featureOrder[i]];
					}	
					eFile << endl;
				}*/
				showFeatures(p,pfeatures,pFile);
			}
			else
			{
				pfeatures=regionFeatures[p];
			}
			/*if(strcmp(p.c_str(),"chr11_116737744_116757175")==0)
			{
				showFeatures(p,pfeatures);
			}*/
			map<int,double> featuresForPair;
			int miss_e=0;
			int miss_p=0;
			for(int f=0;f<featureOrder.size();f++)
			{
				for(int g=0;g<featureOrder.size();g++)
				{
					double featVal_e=0;
					double featVal_p=0;
					bool founde=false;
					bool foundp=false;
					if(pfeatures->find(featureOrder[g])!=pfeatures->end())
					{
						featVal_p=1;
						foundp=true;
					}
					if(efeatures->find(featureOrder[f])!=efeatures->end())
					{
						featVal_e=1;
						founde=true;
					}
					if(!founde)
					{
						miss_e++;
					}
					if(!foundp)	
					{
						miss_p++;
					}
					double f_ep=featVal_e*featVal_p;
					int pos=(f*featureOrder.size())+g;
					featuresForPair[pos]=f_ep;
					if(founde ||foundp)
					{
						if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
						{
							featureFrequency[featureOrder[f]]=1;
						}
						else
						{
							featureFrequency[featureOrder[f]]=featureFrequency[featureOrder[f]]+1;
						}
					}
				}
			}	
			//Filter zeros only if explicitly asked by user. 
			if(filterZeros)
			{
				if(miss_e==featureOrder.size() || miss_p==featureOrder.size())
				{		
					cout <<"Skipping" << e<<"-" << p << endl;
					continue;
				}
			}
			int dist=getDistance(e,p);
			if(dist>1000000)
			{
				continue;
			}
			oFile <<e<<"-"<<p;
			for(int f=0;f<featureOrder.size();f++)
			{
				for(int g=0;g<featureOrder.size();g++)
				{
					int pos=(f*featureOrder.size())+g;
					oFile <<"\t" <<featuresForPair[pos];
				}
			}
			if(correlation)
			{
				double cc=getCCFeature(pfeatures,efeatures);
				oFile <<"\t"  <<cc; 
			}
			if(geneexp.size()>0)
			{
				double exp=getExp(p);
				oFile <<"\t"<<exp;
			}
			//oFile <<"\t" << dist<< "\t" << classLabel<< endl;
			oFile <<"\t" << classLabel<< endl;
			featuresForPair.clear();
		}
	}	
	oFile.close();
	for(int f=0;f<featureOrder.size();f++)
	{
		if(featureFrequency.find(featureOrder[f])==featureFrequency.end())
		{
			cout <<"Get rid of " << featureOrder[f] << endl;	
		}
	}
	return 0;

}


int
Framework::populateFeaturesForRegion(map<string,double>* populateMe, string& key)
{
	if(regionSet.find(key)==regionSet.end())
	{
		cout <<"No region with name " << key<< endl;
		exit(0);
	}
	Region* r=regionSet[key];
	for(map<string,Dataset*>::iterator dIter=datasets.begin();dIter!=datasets.end();dIter++)
	{
		Dataset* d=dIter->second;
		double fval=d->getFeatureVal(r->chromosome,r->begin,r->end);
		if(fval<0)
		{
			continue;
		}
		(*populateMe)[dIter->first]=fval;
	}	
	return 0;
}

int
Framework::showFeatures(string& fName,map<string,double>* fVect,ofstream& oFile)
{
	oFile<< fName;
	for(int i=0;i<featureOrder.size();i++)
	{
		if(fVect->find(featureOrder[i])==fVect->end())
		{
			oFile<<"\t-1";
		}
		else
		{
			oFile<<"\t" <<(*fVect)[featureOrder[i]];
		}
	}
	oFile<< endl;
	return 0;
}

int
Framework::getDistance(string& ekey, string& pkey)
{
	Region* e=regionSet[ekey];
	Region* p=regionSet[pkey];
	int dist=0;
	int first1=e->begin;
	int first2=e->end;
	int second1=p->begin;
	int second2=p->end;
	if(first1>second1)
	{
		first1=p->begin;
		first2=p->end;
		second1=e->begin;
		second2=e->end;
	}
	if(second1<first2) 
	{
		//there is overlap
		dist=0;
	}
	else
	{
		dist=second1-first2;
	}
	return dist;
}

int
main(int argc, const char** argv)
{
	if(argc!=10)
	{
		cout <<"Usage: getFeatures regionPairs featurefile outputfile label leavezeros[yes|no] featuretype[prod|concat|prod_concat] exp[file|null] correlation[yes|no] binarize[binary|continous]" << endl;
		return 0;
	}
	Framework fw;
	fw.readPairs(argv[1]);
	fw.readFeatures(argv[2]);
	if(strcmp(argv[5],"yes")==0)
	{
		fw.setFilterZeros(true);
	}
	else
	{
		fw.setFilterZeros(false);
	}
	if(strcmp(argv[7],"null")!=0)
	{
		fw.readExpression(argv[7]);
	}
	if(strcmp(argv[8],"yes")==0)
	{
		fw.setCorrelation(true);
	}
	else
	{	
		fw.setCorrelation(false);
	}
	bool binary=false;
	if(strcmp(argv[9],"binary")==0)
	{
		binary=true;
	}
	if(strcmp(argv[6],"prod")==0)
	{	
		if(binary)
		{
			fw.generateFeatureFiles_Binary(argv[3],argv[4]);
		}
		else
		{
			fw.generateFeatureFiles(argv[3],argv[4]);
		}
	}
	else if(strcmp(argv[6],"prod_concat")==0)
	{
		//fw.generateFeatureFiles_ConcatProd_Binary(argv[3],argv[4]);
		//fw.generateFeatureFiles_ConcatProd(argv[3],argv[4]);
		cout <<"Not handled" <<endl;
	}
	else if(strcmp(argv[6],"concat")==0)
	{
		if(binary)	
		{
			fw.generateFeatureFiles_Concat_Binary(argv[3],argv[4]);
		}
		else
		{
			fw.generateFeatureFiles_Concat(argv[3],argv[4]);
		}
	}
	else if(strcmp(argv[6],"outerprod")==0)
	{
		if(binary)
		{
			cout <<"Not handled yet" << endl;
		}
		else
		{
			fw.generateFeatureFilesOuterprod(argv[3],argv[4]);
		}
	}
	//fw.generateFeatureFilesDist(argv[3],argv[4]);
	return 0;
}
