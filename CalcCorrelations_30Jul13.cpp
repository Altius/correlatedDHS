#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

struct DHSdata {
  string chrom;
  string beg;
  string end;
  string stringVec;
  vector<double> vec;
  double correlation;
};

void GetMeanVarSD(const vector<double>& vec, double& mean, double& variance, double& SD);
void GetMeanVarSD(const vector<double>& vec, double& mean, double& variance, double& SD)
{
  double N(static_cast<double>(vec.size())), sum(0.);
  if (!N)
    {
      cerr << "Error:  GetMeanVarSD() received an empty vector." << endl;
      exit(1);
    }
  if (1. == N)
    {
      mean = vec[0];
      variance = SD = 0;
      return;
    }
  for (int i = 0; i < vec.size(); i++)
    sum += vec[i];
  mean = sum/N;
  sum = 0.;
  for (int i = 0; i < vec.size(); i++)
    sum += (vec[i] - mean)*(vec[i] - mean);
  variance = sum/(N - 1.);
  SD = sqrt(variance);
}

// compute Pearson's correlation coefficient; return 0 if SD(a) = 0 or SD(b) = 0                                                                                                            
double correlation(const vector<double>& a, const vector<double>& b);
double correlation(const vector<double>& a, const vector<double>& b)
{
  if (a.empty() || b.empty())
    {
      cerr << "Error:  correlation() received an empty vector." << endl;
      exit(1);
    }
  if (a.size() != b.size())
    {
      cerr << "Error:  correlation() received vectors of unequal size ("
           << a.size() << " and " << b.size() << ")." << endl;
      exit(1);
    }
  double mean_a, var_a, SD_a, mean_b, var_b, SD_b, sum(0.), N(static_cast<double>(a.size()));
  GetMeanVarSD(a, mean_a, var_a, SD_a);
  GetMeanVarSD(b, mean_b, var_b, SD_b);
  if (!SD_a || !SD_b)
    {
      /*
      cerr << "Warning:  Attempted to compute correlation with a vector whose SD = 0, between ("
           << a[0];
      for (int j = 1; j < a.size(); j++)
        cerr << ", " << a[j];
      cerr << ") and (" << b[0];
      for (int j = 1; j < a.size(); j++)
        cerr << ", " << b[j];
      cerr << ")." << endl;
*/
      return 0;
    }
  for (int i = 0; i < a.size(); i++)
    sum += (a[i] - mean_a)*(b[i] - mean_b);
  return sum/((N-1.)*SD_a*SD_b);
}

bool vecFromVecString(char *pString, vector<double>& vec);
bool vecFromVecString(char *pString, vector<double>& vec)
{
  char *p = pString;
  vec.clear();
  if (!(p = strtok(pString,",")))
    {
      cerr << "Error:  Failed to find comma-delimited list of numbers in \""
	   << pString << "\"." << endl;
      return false;
    }
  vec.push_back(atof(p));
  while (p = strtok(NULL,","))
    vec.push_back(atof(p));

  return true;
}

bool doEverything(ifstream& ifs, ofstream& ofs);
bool doEverything(ifstream& ifs, ofstream& ofs)
{
  const long BUFSIZE(900000L);//BUFSIZE(550000L);
  char buf[BUFSIZE], *p;
  long linenum(0), fieldnum;
  string promChrom, promBeg, promEnd, promGene, promStringVec;
  vector<double> promVec;
  DHSdata d;
  vector<DHSdata> dvec;

  while (ifs.getline(buf,BUFSIZE))
    {
      linenum++;
      fieldnum = 1;
      if (!(p = strtok(buf,"\t")))
        {
        MissingField:
          cerr << "Error:  Failed to find field " << fieldnum << " on line " << linenum
               << " of the input file." << endl << endl;
          return false;
        }
      promChrom = string(p);
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
        goto MissingField;
      promBeg = string(p);
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
        goto MissingField;
      promEnd = string(p);
      fieldnum++;
      if (!(p = strtok(NULL,"\t")))
        goto MissingField;
      promGene = string(p);
      fieldnum++;
      if (!(p = strtok(NULL,"|")))
        goto MissingField;
      promStringVec = string(p);
      fieldnum++;
      dvec.clear();
      while (p = strtok(NULL,"\t"))
	{
	  d.chrom = string(p);
	  fieldnum++;
	  if (!(p = strtok(NULL,"\t")))
	    goto MissingField;
	  d.beg = string(p);
	  fieldnum++;
	  if (!(p = strtok(NULL,"\t")))
	    goto MissingField;
	  d.end = string(p);
	  fieldnum++;
	  if (!(p = strtok(NULL,";")))
	    goto MissingField;
	  d.stringVec = string(p);
	  fieldnum++;
	  dvec.push_back(d);
	}
      if (dvec.empty())
	continue;

      p = const_cast<char *>(promStringVec.c_str());
      if (!vecFromVecString(p, promVec))
	{
	  cerr << "Error:  On line " << linenum << ", the promoter\'s tag counts are absent, "
	       << "or something else is up.  Exiting..." << endl << endl;
	  return false;
	}
      for (int i = 0; i < dvec.size(); i++)
	{
	  p = const_cast<char *>(dvec[i].stringVec.c_str());
	  if (!vecFromVecString(p, dvec[i].vec))
	    {
	      cerr << "Error:  On line " << linenum << ", DHS number " << i+1
		   << " in the set of DHSs (" << dvec[i].chrom << ':' << dvec[i].beg << '-' << dvec[i].end
		   << "), the string of tag counts is absent, "
		   << "or something else is up.  Exiting..." << endl << endl;
	      return false;
	    }
	  if (dvec[i].vec.size() != promVec.size())
	    {
	      cerr << "Error:  On line " << linenum << ", DHS number " << i+1
		   << " in the set of DHSs (" << dvec[i].chrom << ':' << dvec[i].beg << '-' << dvec[i].end
		   << "), the DHS set of tag counts contains a different number of cell types ("
		   << dvec[i].vec.size() << ") than the promoter set does ("
		   << promVec.size() << ").  Exiting..." << endl << endl;
	      return false;
	    }
	  dvec[i].correlation = correlation(promVec, dvec[i].vec);

	  ofs << promChrom << '\t' << promBeg << '\t' << promEnd << '\t' << promGene << '\t'
	      << dvec[i].chrom << '\t' << dvec[i].beg << '\t' << dvec[i].end << '\t' << dvec[i].correlation
	      << endl;
	}

    }

  return true;
}

int main(int argc, const char* argv[])
{
  if (3 != argc)
    {
      cerr << "Usage:  " << argv[0] << " infile outfile\n"
	   << "\twhere infile is a bed file with the promoter DHS and gene name and comma-delimited tag sums in field 1\n"
	   << "\t(followed immediately by \'|\'), immediately followed by a non-empty semicolon-delimited list\n"
	   << "\tof distal DHSs and their comma-delimited tag sums.\n"
	   << "\toutfile is a .bed8 file with no \'|\' characters in it,\n"
	   << "\tjust the promoter DHS, gene name, distal DHS, correlation,\n"
	   << "\twith one unique pair per line (note:  duplicate pairs with different gene names will occur)."
	   << endl << endl;
      return -1;
    }
  ifstream ifs(argv[1]);
  if (!ifs)
    {
      cerr << "Error:  Unable to open file \"" << argv[1] << "\" for read." << endl << endl;
      return -1;
    }
  ofstream ofs(argv[2]);
  if (!ofs)
    {
      cerr << "Error:  Unable to open file \"" << argv[2] << "\" for write." << endl << endl;
      return -1;
    }
  if (!doEverything(ifs, ofs))
    return -1;

  return 0;
}

