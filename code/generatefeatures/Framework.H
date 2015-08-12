#ifndef _FRAMEWORK_
#define _FRAMEWORK_
#include <map>
#include <string>
#include <vector>
#include <fstream>
using namespace std;
class Dataset;
class Framework
{
	public:
	Framework();
	~Framework();
	int readPairs(const char*);
	//Must include features in required order
	int readFeatures(const char*);
	int readExpression(const char*);
	int generateFeatureFiles(const char*,const char*);
	int generateFeatureFiles_ConcatProd(const char* aFName, const char* classLabel);
	int generateFeatureFiles_Concat(const char* aFName, const char* classLabel);
	int generateFeatureFiles_ConcatProd_Binary(const char* aFName, const char* classLabel);
	int generateFeatureFiles_Concat_Binary(const char* aFName, const char* classLabel);
	int generateFeatureFiles_Binary(const char*, const char*);
	int generateFeatureFilesDist(const char*,const char*);
	int generateFeatureFilesOuterprod(const char*, const char*);
	int generateFeatureFilesOuterprod_Binary(const char*, const char*);
	int setFilterZeros(bool);
	int setCorrelation(bool);
	double getCCFeature(map<string,double>*, map<string,double>*);
	double getCCFeature_Bin(map<string,double>*, map<string,double>*);
	double getExp(string&);
	struct Region
	{
		string chromosome;
		int begin;
		int end;	
	};
	private:
	int showFeatures(string&,map<string,double>*,ofstream&);
	int populateFeaturesForRegion(map<string,double>*,string&);
	map<string,map<string,double>*> regionFeatures;
	map<string,Dataset*> datasets;
	vector<string> featureOrder;
	map<string,map<string,int>*> pairSet;
	map<string,Region*> regionSet;
	int getDistance(string&,string&);
	bool filterZeros;
	map<string,double> geneexp;
	bool correlation;
};

#endif
