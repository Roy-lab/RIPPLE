#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "Dataset.H"

Dataset::Dataset()
{
}

Dataset::~Dataset()
{
}

int 
Dataset::readDataSet(const char* aFName)
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
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string chrom;
		int begin;
		int end;
		double qval=-1;
		double pval=-1;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				chrom.append(tok);
			}
			else if(tokCnt==1)
			{
				begin=atoi(tok);
			}
			else if(tokCnt==2)
			{
				end=atoi(tok);	
			}
			else if(tokCnt==7)
			{
				pval=atof(tok);
			}
			else if(tokCnt==8)
			{
				qval=atof(tok);
				if(qval<0)
				{
					qval=-1;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}	
		Dataset::Peak* r=new Dataset::Peak;	
		r->start=begin;
		r->end=end;
		if(qval>0)
		{
			r->lqval=qval;
		}
		else
		{
			r->lqval=pval;
		}
		map<int,Peak*>* rSet=NULL;
		if(peakSet.find(chrom)==peakSet.end())
		{
			rSet=new map<int,Peak*>;
			peakSet[chrom]=rSet;
		}
		else
		{
			rSet=peakSet[chrom];
		}
		if(rSet->find(begin)==rSet->end())
		{
			(*rSet)[begin]=r;
		}
		else
		{
			//If for some reason the other dataset has a hit at he same peak location, then use the greater value
			Peak* oldpeak=(*rSet)[begin];
			if(r->lqval>oldpeak->lqval)
			{
				(*rSet)[begin]=r;
			}
		}
	}
	inFile.close();
	cout <<"Done reading "<< peakSet.size() << " different chromosomes"<<endl;
	return 0;
}

//Return the value with the greatest overlap
double
Dataset::getFeatureVal(string& chrom, int begin, int end)
{
	if(peakSet.find(chrom)==peakSet.end())
	{
		return -1;
	}	
	map<int,Peak*>* peaks=peakSet[chrom];
	double maxVal=-1;
	double peakcnt=0;
	for(map<int,Peak*>::iterator pIter=peaks->begin();pIter!=peaks->end();pIter++)
	{
		Peak* peak=pIter->second;
		int begin1=begin;
		int end1=end;
		int begin2=peak->start;
		int end2=peak->end;
		if(begin1>begin2)
		{
			begin2=begin;
			end2=end;
			begin1=peak->start;
			end1=peak->end;
		}
		if(end1<begin2)
		{
			continue;
		}
		if(maxVal<peak->lqval)
		{
			maxVal=peak->lqval;
			peakcnt++;
		}
	}
	return maxVal;
	//return peakcnt;
}
