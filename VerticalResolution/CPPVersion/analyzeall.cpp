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

void CheckIR(vector<double> vCoef)
{
	vector<double> impulse(2400, 0.); // vector of 2400 values initialized at 1
	unsigned status = 0;
	double rDz = 0.;
	vector<double> irout;
	NDAAC_ResolIR(100, vCoef, impulse, status, rDz, irout);
	cout<<"dzIR = "<< rDz<<endl;
}


void CheckDF(vector<double> vCoef)
{
	vector<double> impulse(2400, 1.); // vector of 2400 values initialized at 1
	unsigned status = 0;
	double rDz = 0.;
	vector<double> irout, irout2;
	NDAAC_ResolDF(100, vCoef, impulse, status, rDz, irout);
	cout<<"dzDF = "<< rDz<<endl;
}


void PrintVector(vector<double> vCoef)
{
	cout<<"Coeff = ";

	for(unsigned i = 0 ; i < vCoef.size(); ++i)
	{
		cout<<vCoef[i]<<",";
	}
	cout<<endl;
}


void check0()
{
	vector<double> weight;
	weight.push_back(0.333333333333);
	weight.push_back(0.333333333333);
	weight.push_back(0.333333333333);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check1(){
	vector<double> weight;
	weight.push_back(0.2);
	weight.push_back(0.2);
	weight.push_back(0.2);
	weight.push_back(0.2);
	weight.push_back(0.2);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check2(){
	vector<double> weight;
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	weight.push_back(0.111111111111);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check3(){
	vector<double> weight;
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	weight.push_back(0.0588235294118);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check4(){
	vector<double> weight;
	weight.push_back(0.0);
	weight.push_back(1.0);
	weight.push_back(0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check5(){
	vector<double> weight;
	weight.push_back(0.0);
	weight.push_back(0.25);
	weight.push_back(0.5);
	weight.push_back(0.25);
	weight.push_back(0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check6(){
	vector<double> weight;
	weight.push_back(0.0);
	weight.push_back(0.036611652);
	weight.push_back(0.125);
	weight.push_back(0.21338835);
	weight.push_back(0.25);
	weight.push_back(0.21338835);
	weight.push_back(0.125);
	weight.push_back(0.036611652);
	weight.push_back(0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check7(){
	vector<double> weight;
	weight.push_back(0.0);
	weight.push_back(0.0047575292);
	weight.push_back(0.018305826);
	weight.push_back(0.038582285);
	weight.push_back(0.0625);
	weight.push_back(0.086417715);
	weight.push_back(0.10669417);
	weight.push_back(0.12024247);
	weight.push_back(0.125);
	weight.push_back(0.12024247);
	weight.push_back(0.10669417);
	weight.push_back(0.086417715);
	weight.push_back(0.0625);
	weight.push_back(0.038582285);
	weight.push_back(0.018305826);
	weight.push_back(0.0047575292);
	weight.push_back(0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check8(){
	vector<double> weight;
	weight.push_back(0.0);
	weight.push_back(1.0);
	weight.push_back(0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check9(){
	vector<double> weight;
	weight.push_back(-0.0857142857143);
	weight.push_back(0.342857142857);
	weight.push_back(0.485714285714);
	weight.push_back(0.342857142857);
	weight.push_back(-0.0857142857143);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check10(){
	vector<double> weight;
	weight.push_back(-0.0909090909091);
	weight.push_back(0.0606060606061);
	weight.push_back(0.168831168831);
	weight.push_back(0.233766233766);
	weight.push_back(0.255411255411);
	weight.push_back(0.233766233766);
	weight.push_back(0.168831168831);
	weight.push_back(0.0606060606061);
	weight.push_back(-0.0909090909091);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check11(){
	vector<double> weight;
	weight.push_back(-0.0650154798762);
	weight.push_back(-0.0185758513932);
	weight.push_back(0.0216718266254);
	weight.push_back(0.0557275541796);
	weight.push_back(0.0835913312693);
	weight.push_back(0.105263157895);
	weight.push_back(0.120743034056);
	weight.push_back(0.130030959752);
	weight.push_back(0.133126934985);
	weight.push_back(0.130030959752);
	weight.push_back(0.120743034056);
	weight.push_back(0.105263157895);
	weight.push_back(0.0835913312693);
	weight.push_back(0.0557275541796);
	weight.push_back(0.0216718266254);
	weight.push_back(-0.0185758513932);
	weight.push_back(-0.0650154798762);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check12(){
	vector<double> weight;
	weight.push_back(0);
	weight.push_back(1);
	weight.push_back(0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check13(){
	vector<double> weight;
	weight.push_back(-0.0);
	weight.push_back(0.20689655);
	weight.push_back(0.5862069);
	weight.push_back(0.20689655);
	weight.push_back(-0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check14(){
	vector<double> weight;
	weight.push_back(-0.0);
	weight.push_back(0.010552849);
	weight.push_back(0.10036839);
	weight.push_back(0.2372394);
	weight.push_back(0.30367873);
	weight.push_back(0.2372394);
	weight.push_back(0.10036839);
	weight.push_back(0.010552849);
	weight.push_back(-0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check15(){
	vector<double> weight;
	weight.push_back(-0.0);
	weight.push_back(-0.00082412208);
	weight.push_back(0.0036995271);
	weight.push_back(0.020050227);
	weight.push_back(0.048719477);
	weight.push_back(0.084828255);
	weight.push_back(0.12013351);
	weight.push_back(0.14580285);
	weight.push_back(0.15518056);
	weight.push_back(0.14580285);
	weight.push_back(0.12013351);
	weight.push_back(0.084828255);
	weight.push_back(0.048719477);
	weight.push_back(0.020050227);
	weight.push_back(0.0036995271);
	weight.push_back(-0.00082412208);
	weight.push_back(-0.0);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check16(){
	vector<double> weight;
	weight.push_back(-0.5);
	weight.push_back(0.0);
	weight.push_back(0.5);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check17(){
	vector<double> weight;
	weight.push_back(-0.2);
	weight.push_back(-0.1);
	weight.push_back(0.0);
	weight.push_back(0.1);
	weight.push_back(0.2);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check18(){
	vector<double> weight;
	weight.push_back(-0.0666666666667);
	weight.push_back(-0.05);
	weight.push_back(-0.0333333333333);
	weight.push_back(-0.0166666666667);
	weight.push_back(0.0);
	weight.push_back(0.0166666666667);
	weight.push_back(0.0333333333333);
	weight.push_back(0.05);
	weight.push_back(0.0666666666667);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check19(){
	vector<double> weight;
	weight.push_back(-0.0196078431373);
	weight.push_back(-0.0171568627451);
	weight.push_back(-0.0147058823529);
	weight.push_back(-0.0122549019608);
	weight.push_back(-0.00980392156863);
	weight.push_back(-0.00735294117647);
	weight.push_back(-0.00490196078431);
	weight.push_back(-0.00245098039216);
	weight.push_back(0.0);
	weight.push_back(0.00245098039216);
	weight.push_back(0.00490196078431);
	weight.push_back(0.00735294117647);
	weight.push_back(0.00980392156863);
	weight.push_back(0.0122549019608);
	weight.push_back(0.0147058823529);
	weight.push_back(0.0171568627451);
	weight.push_back(0.0196078431373);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check20(){
	vector<double> weight;
	weight.push_back(0.0833333333333);
	weight.push_back(-0.666666666667);
	weight.push_back(0.0);
	weight.push_back(0.666666666667);
	weight.push_back(-0.0833333333333);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check21(){
	vector<double> weight;
	weight.push_back(0.0723905723906);
	weight.push_back(-0.119528619529);
	weight.push_back(-0.162457912458);
	weight.push_back(-0.106060606061);
	weight.push_back(0.0);
	weight.push_back(0.106060606061);
	weight.push_back(0.162457912458);
	weight.push_back(0.119528619529);
	weight.push_back(-0.0723905723906);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}
void check22(){
	vector<double> weight;
	weight.push_back(0.0321637426901);
	weight.push_back(-0.00421396628827);
	weight.push_back(-0.0276487788098);
	weight.push_back(-0.0399896800826);
	weight.push_back(-0.0430856553148);
	weight.push_back(-0.0387856897145);
	weight.push_back(-0.0289387684899);
	weight.push_back(-0.015393876849);
	weight.push_back(0.0);
	weight.push_back(0.015393876849);
	weight.push_back(0.0289387684899);
	weight.push_back(0.0387856897145);
	weight.push_back(0.0430856553148);
	weight.push_back(0.0399896800826);
	weight.push_back(0.0276487788098);
	weight.push_back(0.00421396628827);
	weight.push_back(-0.0321637426901);
	PrintVector(weight);
	CheckIR(weight);
	CheckDF(weight);
}

int main(int argc, char** argv)
{
	cout<<"We start at degree 0 "<<endl;

	cout<<"3 points"<<endl;
	check0();
	cout<<"5 points"<<endl;
	check1();
	cout<<"9 points"<<endl;
	check2();
	cout<<"17 points"<<endl;
	check3();

	cout<<"We are at degree 0, hanning "<<endl;

	cout<<"3 points"<<endl;
	check4();
	cout<<"5 points"<<endl;
	check5();
	cout<<"9 points"<<endl;
	check6();
	cout<<"17 points"<<endl;
	check7();
	cout<<"We are at degree 2 "<<endl;

	cout<<"3 points"<<endl;
	check8();
	cout<<"5 points"<<endl;
	check9();
	cout<<"9 points"<<endl;
	check10();
	cout<<"17 points"<<endl;
	check11();


	cout<<"We are at degree 2, hanning "<<endl;

	cout<<"3 points"<<endl;
	check12();
	cout<<"5 points"<<endl;
	check13();
	cout<<"9 points"<<endl;
	check14();
	cout<<"17 points"<<endl;
	check15();

	cout<<"We are at degree 1, derivative "<<endl;

	cout<<"3 points"<<endl;
	check16();
	cout<<"5 points"<<endl;
	check17();
	cout<<"9 points"<<endl;
	check18();
	cout<<"17 points"<<endl;
	check19();

	cout<<"We are at degree 2, derivative "<<endl;

	cout<<"5 points"<<endl;
	check20();
	cout<<"9 points"<<endl;
	check21();
	cout<<"17 points"<<endl;
	check22();






}

