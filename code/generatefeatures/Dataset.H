#ifndef _DATASET
#define _DATASET

#include <map>
#include <string>
using namespace std;
class Dataset
{
	public:
	Dataset();
	~Dataset();
	int readDataSet(const char*);
	double getFeatureVal(string& chrom, int begin, int end);
	struct Peak
	{
		int start;
		int end;
		double lqval;
		double signal;
		int peakCnt;
	};
	private:
	map<string,map<int,Peak*>*> peakSet;
};
#endif
