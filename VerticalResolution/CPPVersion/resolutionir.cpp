#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <functional>
#include <complex>
#include <algorithm>
using namespace std;



double stddev(std::vector<double> v)
{
	if (v.size()  <= 2)
		return 0.;
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double m =  sum / v.size();

	double accum = 0.0;
	for(std::vector<double>::iterator d = v.begin(); d != v.end(); ++d)
	{
		accum += (*d - m) * (*d - m);
	}

	double stdev = sqrt(accum / (v.size()-1));
	return stdev;
}

void convolve(const vector<double>& vIrc, const vector<double>& vCo, vector<double> & ir)
{
	ir.resize(vIrc.size(), 0.);
	int nc = (vCo.size() - 1) / 2;
	for(unsigned j = nc; j < vIrc.size() - nc; ++j)
	{
		double tmp = 0;
		for(int a = - nc; a < nc + 1; ++a)
		{
			tmp += vIrc[j + a] * vCo[nc + a];
		}
		ir[j] = tmp;
	}
}


void NDAAC_ResolIR(double vZsampling, const vector<double>& vCoef, vector<double> vIRinp, unsigned &rStatus, double &rDzfwhm, vector<double>& rIRout)
{

	unsigned ncoef = vCoef.size();
	unsigned nk = vIRinp.size();
	rStatus = 0;
	double missval = -99;
	vector<double> irout(vIRinp.begin(), vIRinp.end());
	vector<double> irconvol(vIRinp.begin(), vIRinp.end());
	vector<double> ir(nk, 0.);
	int nc = (ncoef - 1) / 2;
	if((ncoef % 2) == 0)
	{
		cout<<"ERROR: the number of coefficients is even"<<endl;
		return;
	}
	if ((ncoef == 1) and vCoef[0] != 1)
	{
		rStatus = 0;
		cout<<"ERROR: we have a non-filter with a non unit coefficient"<<endl;
		return;
	}
	if ((ncoef == 1) and vCoef[0] == 1)
	{
		rStatus = 1;
	}else
	{
		rStatus = 42;
		for(int ic = 1; ic < nc + 1; ++ic)
		{
			if(abs(vCoef[nc - ic]) < 1E-20 && abs(vCoef[nc + ic]) < 1E-20)
			{
				cout<<"zero  "<<endl;
			}else if (abs(vCoef[nc - ic] - vCoef[nc + ic]) < 1E-3 * abs(vCoef[nc + ic]))
			{
				if (rStatus == 2)
				{
					rStatus = 0;
					cout<<"ERROR: we have a mix of symmetric and antisymmetric"<<endl;
					return;
				}
				rStatus = 1;
			}else if (abs(vCoef[nc - ic] + vCoef[nc + ic]) < 1E-3 * abs(vCoef[nc + ic]) and vCoef[nc] == 0)
			{
				if (rStatus == 1)
				{
					rStatus = 0;
					cout<<"ERROR: we have a mix of symmetric and antisymmetric"<<endl;
					return;
				}
				rStatus = 2;
			}else
			{
				rStatus = 0;
				cout<<"ERROR: we have a mix of symmetric and antisymmetric (impossible parameters)"<<endl;
				return;
			}

		}	
	}
	if (accumulate(irconvol.begin(), irconvol.end(),0.) == 0 and stddev(irconvol) == 0)
	{
		irconvol[nk / 2] = 1;
	}
	if (rStatus == 2)
	{
		cout<<"Status 2"<<endl;
		vector<double> irtemp(irconvol.size(), 0.);
		irtemp[0] =irconvol[0];
		double maxirtemp = irtemp[0];
		for(unsigned k=0; k < nk; ++k)
		{
			//irtemp[k] = sum(irconvol[0:k + 1]);
			irtemp[k] = accumulate(irconvol.begin(), irconvol.begin() + k + 1, 0); 
			if (irtemp[k] > maxirtemp)
			{
				maxirtemp = irtemp[k];
			}
		}
		for(unsigned k=0; k< irconvol.size(); ++k)
		{
			irconvol[k] = irtemp[k] / maxirtemp;
		}
	}


	convolve(irconvol, vCoef, ir);

	double maxir = *(std::max_element(ir.begin(), ir.end()));
	if(maxir < 0)
	{
		for(unsigned i = 0; i < ir.size(); ++i)
		{
			ir[i] = abs(ir[i]);
		}
		maxir = *(std::max_element(ir.begin(), ir.end()));
	}
	rIRout.resize(ir.size(), 0.);
	copy(ir.begin(), ir.end(), rIRout.begin());
	for(unsigned i = 0; i < ir.size(); ++i)
	{
		irout[i] = (ir[i]) / maxir;
	}

	unsigned k1 = 0;
	double x1 = 0;
	double x2 = 0;
	double irc = 0.5;// pow(10, -6 / 20.)


	for(unsigned jf = 0; jf < nk - 1; ++jf)
	{
		if((irout[jf + 1] - irc) * (irout[jf] - irc) <= 0)
		{
			if(abs(irout[jf] - irc) < 1E-9)
			{
				k1 = jf;
				x1 = k1;
			}else if (abs(irout[jf + 1] - irc) < 1E-9)
			{

				k1 = jf + 1;
				x1 = k1;
			}else
			{
				double weight = abs(irout[jf] - irc) / abs(irout[jf + 1] - irout[jf]);
				k1 = jf;
				x1 = jf + weight;
			}
			break;
		}
	}

	for(unsigned jf = k1 + 1; jf < nk - 1; ++jf)
	{
		if((irout[jf + 1] - irc) * (irout[jf] - irc) <= 0)
		{
			if(abs(irout[jf + 1] - irc) < 1E-9)
			{
				x2 = jf + 1;
			}else if (abs(irout[jf] - irc) < 1E-9)
			{

				x2 = jf;
			}else
			{
				double weight = abs(irout[jf] - irc) / abs(irout[jf + 1] - irout[jf]);
				x2 = jf + weight;
			}
			break;
		}
	}


	rDzfwhm = vZsampling * (x2 - x1); 
}







int main(int argc, char** argv)
{

	vector<double> impulse(2400, 0.); // vector of 2400 values initialized at 1

//vector<double> filter(9, 1/9.); // vector of 9 values initialized at 1/9.
	vector<double> filter(3,0); // vector of 9 values initialized at 1/9.
	filter[0] = 0.333;
	filter[1] = 0.333;
	filter[2] = 0.333;

	unsigned status = 0;
	double rDz = 0.;
	vector<double> irout;

	NDAAC_ResolIR(100, filter, impulse, status, rDz, irout);

	cout<<"The status after the computation is "<<status<< "  and the calculated resolution is "<<rDz<<endl;

	cout<<"The output function is "<<endl;

	for(unsigned i = 0; i < irout.size(); ++i)
	{
		cout<<irout[i]<<"\t";
		//cout<<impulse[i]<<"\t";
	}
	cout<<endl;
	return 0;
}
