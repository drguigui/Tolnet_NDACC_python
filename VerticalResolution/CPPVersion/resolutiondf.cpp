#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <functional>
#include <complex>
using namespace std;


/**
 * Computes the standard deviation of an input vector
 *
 *
 */
double stddev(std::vector<double> v)
{
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double m =  sum / v.size();

	double accum = 0.0;
//	std::for_each (std::begin(v), std::end(v), [&](const double d) {a
	for(std::vector<double>::iterator d = v.begin(); d != v.end(); ++d)
	{
		accum += (*d - m) * (*d - m);
	}

	double stdev = sqrt(accum / (v.size()-1));
	return stdev;
}



void NDAAC_ResolDF(double vZsampling, vector<double> vCoef, vector<double> vHfinp, unsigned &rStatus, double &rDzfwhm, vector<double>& rHfout)
{

	unsigned ncoef = vCoef.size();
	unsigned nf = vHfinp.size();
	rStatus = 0;
	double missval = -99;
	vector<double> hfinp(vHfinp.begin(), vHfinp.end());
	vector<double> hf(nf, 0.);
	vector<double> hfout(nf, 0.);
	double dzcutf = missval;
	double fc = missval;
	double fcinp = missval;
	unsigned nc = (ncoef - 1) / 2;
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
		for(unsigned ic = 1; ic < nc + 1; ++ic)
		{
			if( (abs(vCoef[nc - ic]) < 1E-20) && (abs(vCoef[nc + ic]) < 1E-20))
			{
	//			cout<<"zero and "<<rStatus<<endl;
				cout<<"zero"<<endl;
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
				cout<<"ERROR: we have a mix of symmetric and antisymmetric"<<endl;
				return;
			}

		}	
	}

	vector<double> f;
	for(unsigned i = 1; i < nf + 1; ++i)
	{
		f.push_back(double(i) * 0.5 / double(nf));
	}

	for(unsigned jf = 0; jf < nf - 1; ++jf)
	{
		if ((hfinp[jf + 1] - 0.5) * (hfinp[jf] - 0.5) <= 0)
		{
			if (abs(hfinp[jf] - 0.5) < 1E-9)
			{
				fcinp = f[jf];
			}else if (abs(hfinp[jf + 1] - 0.5) < 1E-9)
			{
				fcinp = f[jf + 1];
			}else
			{
				double weight = abs(hfinp[jf] - 0.5) / abs(hfinp[jf + 1] - hfinp[jf]);
				fcinp = f[jf] + weight * abs(f[jf+1] - f[jf]);
			}
			break;
		}
	}
	if (accumulate(hfinp.begin(), hfinp.end(),0.) == nf and stddev(hfinp) == 0)
	{
		fcinp = 0.5;
	}

	vector<complex<double> > hfcomp(nf, 0.);

	for(unsigned jf = 0; jf < nf ; ++ jf)
	{
		complex<double> summ =0.;
		for(int ic = -int(nc); ic < int(nc) + 1; ++ic)
		{
			complex<double> bla(cos(2 * M_PI * f[jf] * ic), - sin(2 * M_PI * f[jf] * ic));
			summ += vCoef[nc - ic] * bla;
		}
		hfcomp[jf] = summ;
	}


	if(rStatus == 2)
	{
		for(unsigned i = 0; i < hfcomp.size(); ++i)
		{
			double d = std::imag(hfinp[i] * hfcomp[i]) / (2 * M_PI * f[i]);
			hfout[i] = d;
		}
	}else
	{
		for(unsigned i = 0; i < hfcomp.size(); ++i)
		{
			double d = std::real(hfinp[i] * hfcomp[i]);
			hfout[i] = d;
		}
	}


	double hfc = 0.5; // pow(10, -6 / 20.)

	for(unsigned jf = 0; jf < nf - 1; ++jf)
	{
		if((hfout[jf + 1] - hfc) * (hfout[jf] - hfc) <= 0)
		{
			if(abs(hfout[jf] - hfc) < 1E-9)
			{
				fc = f[jf];
			}else if (abs(hfout[jf + 1] - hfc) < 1E-9)
			{
				fc = f[jf + 1];
			}else
			{
				double weight = abs(hfout[jf] - hfc) / abs(hfout[jf + 1] - hfout[jf]);
				fc = f[jf] + weight * abs(f[jf+1] - f[jf]);
			}
			break;
		}
	}

	if(ncoef == 1)
	{
		fc = fcinp;
	}
	if (fc > 0)
	{
		rDzfwhm = vZsampling / (2 * fc);
	}
	rHfout.resize(hfout.size());
	std::copy(hfout.begin(), hfout.end(), rHfout.begin());
}







int main(int argc, char** argv)
{

	vector<double> impulse(2400, 1.); // vector of 2400 values initialized at 1

	vector<double> filter(25, 1/25.); // vector of 9 values initialized at 1/9.

	vector<double> filter2(3,0); // vector of 9 values initialized at 1/9.
	filter2[0] = -1/2.;
	filter2[1] = 0;
	filter2[2] = 1/2.;
	unsigned status = 0;
	double rDz = 0.;
	vector<double> irout, irout2;

	NDAAC_ResolDF(1, filter, impulse, status, rDz, irout);

	cout<<"The status after the computation is "<<status<< "  and the calculated resolution is "<<rDz<<endl;



	NDAAC_ResolDF(1, filter2, irout, status, rDz, irout2);
	cout<<"The status after the computation is "<<status<< "  and the calculated resolution is "<<rDz<<endl;
	cout<<"The output function is "<<endl;




	for(unsigned i = 0; i < irout2.size(); ++i)
	{
		cout<<irout2[i]<<"\t";
		//cout<<impulse[i]<<"\t";
	}
	cout<<endl;
	return 0;
}
