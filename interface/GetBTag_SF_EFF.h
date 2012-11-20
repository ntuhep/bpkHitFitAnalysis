#ifndef GETBTAG_SF_EFF
#define GETBTAG_SF_EFF

#include "TMath.h"
#include "MyAna/bprimeKit/interface/format.h"
#include <iostream>
//class JetInfoBranches;
float _btag_sfb_ptmin[14] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
float _btag_sfb_ptmax[14] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};
float _sfb_error[14] = {
 0.0311456,
 0.0303825,
 0.0209488,
 0.0216987,
 0.0227149,
 0.0260294,
 0.0205766,
 0.0227065,
 0.0260481,
 0.0278001,
 0.0295361,
 0.0306555,
 0.0367805,
 0.0527368 };

class BTag_Utility {
	public :
	BTag_Utility(JetInfoBranches jet, double sigma);
	double GetSF(int ijet);
	double GetEFF(int ijet);
	double Btag_SF;
	double Btag_SF_ERR;
	double Btag_EFF;

	private :
	JetInfoBranches    _jet;
	double _sigma;
	float btag_sfb_ptmin[14];
	float btag_sfb_ptmax[14];
	float sfb_error[14];
	
};

BTag_Utility::BTag_Utility(JetInfoBranches jet, double sigma)
	:_jet(jet),
	_sigma(sigma)
{
	for(int i=0; i < 14; i++){
		btag_sfb_ptmin[i] = _btag_sfb_ptmin[i];
		btag_sfb_ptmax[i] = _btag_sfb_ptmax[i];
		sfb_error[i] = _sfb_error[i];	
	}
}

double BTag_Utility::GetSF(int ijet){
	double scaling_factor = -999;
	double sf_error = -999;
	int flavor = abs(_jet.GenFlavor[ijet]);
	double x = _jet.Pt[ijet];
	double temp_x = x;
	double eta = fabs(_jet.Eta[ijet]);
	if(flavor == 0){flavor = 1;}
	if(flavor == 5 || flavor == 4){

		if(x > 670) x = 670;
		if(x < 30) x = 30;// will result in a peak around 0.9368
		scaling_factor = 0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));

		int i;
		for(i = 0; i < 14; i++){
			if (btag_sfb_ptmin[i] <= temp_x && temp_x < btag_sfb_ptmax[i]) break;
		}
		if(i != 14){
			sf_error = sfb_error[i];
		}	
		if(i == 14){
			if(temp_x == 670) sf_error = sfb_error[13];
			if(temp_x > 670) sf_error = 2*sfb_error[13];
			if(temp_x < 30) sf_error = sfb_error[0] + 0.12;
		}
		if(flavor == 4) sf_error = 2*sf_error;
		sf_error = sf_error*_sigma;
	}
		
	else if( flavor > 0 && flavor != 5 && flavor != 4){
		if(x > 670) x = 670;
		if(x < 20) x = 20;
		if (eta >= 0 && eta <= 0.8){
			scaling_factor = ((1.06182+(0.000617034*x))+(-1.5732e-06*(x*x)))+(3.02909e-10*(x*(x*x)));
			if(_sigma < 0) sf_error = ((0.962627+(0.000448344*x))+(-1.25579e-06*(x*x)))+(4.82283e-10*(x*(x*x)));
			if(_sigma >= 0) sf_error = ((1.12368+(0.00124806*x))+(-3.9032e-06*(x*x)))+(2.80083e-09*(x*(x*x)));
		}
		if (eta > 0.8 && eta <= 1.6){
			scaling_factor = ((1.111+(-9.64191e-06*x))+(1.80811e-07*(x*x)))+(-5.44868e-10*(x*(x*x)));//lots of low pt 20 ~ 40 jet result in a peak around 1.1108
			if(_sigma < 0) sf_error = ((1.02055+(-0.000378856*x))+(1.49029e-06*(x*x)))+(-1.74966e-09*(x*(x*x)));
			if(_sigma >= 0) sf_error = ((1.20146+(0.000359543*x))+(-1.12866e-06*(x*x)))+(6.59918e-10*(x*(x*x)));
		}
		if (eta > 1.6 /*&& eta <= 2.4*/){
			scaling_factor = ((1.08498+(-0.000701422*x))+(3.43612e-06*(x*x)))+(-4.11794e-09*(x*(x*x)));
			if(_sigma < 0) sf_error = ((0.983476+(-0.000607242*x))+(3.17997e-06*(x*x)))+(-4.01242e-09*(x*(x*x)));
			if(_sigma >= 0) sf_error = ((1.18654+(-0.000795808*x))+(3.69226e-06*(x*x)))+(-4.22347e-09*(x*(x*x)));
		}
		sf_error = (sf_error-scaling_factor)*fabs(_sigma);
		if (temp_x > 670) {
			sf_error = 2*sf_error; 
		}
		if (temp_x < 20) {
			sf_error = 2*sf_error;
			//std::cout << "No SF Light Quark available for Jet Pt < 20, using Pt = 20 SF with twice error instead " << std::endl;
		}
	}
	Btag_SF = scaling_factor;
	Btag_SF_ERR = sf_error;
	return scaling_factor;
}

double BTag_Utility::GetEFF(int ijet){
	double efficiency = -999;
	int flavor = abs(_jet.GenFlavor[ijet]);
	double eta = fabs(_jet.Eta[ijet]);
	//double x = _jet.CombinedSVBJetTags[ijet];
	double x = _jet.Pt[ijet];
	if (flavor == 0) {flavor = 1;}
	if (flavor == 5) {
		x = 0.679;
		efficiency = -6.41591823466*x*x*x*x +  11.5812173893*x*x*x +  -6.94406444589*x*x +  1.13278339944*x +  0.889359753365;
		//if(_jet.Pt[ijet] > 200) std::cout << " Jet Pt > 200, still applying the Efficiency, Be aware!" << std::endl;
		//if(_jet.Pt[ijet] < 30) std::cout << " Jet Pt < 30, still applying the Efficiency, Be aware!" << std::endl;
		//if(eta > 2.4) std::cout << "Jet Eta > 2.4, Sapplying B Quark Efficiency, Be awared!" << std::endl;
	}
	if (flavor == 4) {
		x = 0.679;
		efficiency = -1.73338329789*x*x*x*x +  1.26161794785*x*x*x +  0.784721653518*x*x +  -1.03328577451*x +  1.04305075822;
		//if(_jet.Pt[ijet] > 200) std::cout << " Jet Pt > 200, still applying the Efficiency, Be aware!" << std::endl;
		//if(_jet.Pt[ijet] < 30) std::cout << " Jet Pt < 30, still applying the Efficiency, Be aware!" << std::endl;
		//if(eta > 2.4) std::cout << "Jet Eta > 2.4, Sapplying C Quark Efficiency, Be awared!" << std::endl;
	}
	else if (flavor > 0 && flavor != 5 && flavor != 4) {
		if (_jet.Pt[ijet] > 670) {
			x = 670;
			//std::cout << "Jet Pt > 670, applying the Light Quark Efficiency as if pt = 670, Be awared!" << std::endl;
		}
		if (_jet.Pt[ijet] < 20) {
			x = 20;
			//std::cout << "Jet Pt < 20, applying the Light Quark Efficiency as if pt = 20, Be awared!" << std::endl;
		}
		if (eta >= 0 && eta <= 0.8){
			efficiency = (0.00967751+(2.54564e-05*x))+(-6.92256e-10*(x*x));
		}
		if (eta > 0.8 && eta <= 1.6){
			efficiency = (0.0113428+(5.18983e-05*x))+(-2.59881e-08*(x*x));
		}
		if (eta > 1.6 && eta <= 2.4){
			efficiency = (0.013595+(0.000104538*x))+(-1.36087e-08*(x*x));
		}
		if (eta > 2.4) {
			efficiency = (0.013595+(0.000104538*x))+(-1.36087e-08*(x*x));
			//std::cout << "Jet Eta > 2.4, applying Light Quark Efficiency as if eta within 1.6~2.4, Be awared!" << std::endl; 
		}	
	}
	Btag_EFF = efficiency;
	return efficiency;
}
#endif
