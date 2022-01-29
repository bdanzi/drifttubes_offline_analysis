//plot funzione trasferimento 
#include "funcUtils.h"


void transFunc(double min_f, double max_f, int npt,double omega_0=50e+6, double width=1e+6, int type=0){
  wavefft white(2*npt);
  double step_f = (max_f-min_f)/((double)npt);
  for(int ipt=0;ipt<npt;++ipt){
    white.real[ipt]=1.0;
    white.img[ipt]=0.0;
    white.ampl[ipt]=1.0;
    white.phi[ipt]=0.0;
    white.omega[ipt]=TMath::TwoPi()*(min_f+step_f*ipt);  
  }
  double *realFltFFT = new double[2*npt];
  double *imgFltFFT = new double[2*npt];
  if (type<4) { filterNotch(white, realFltFFT, imgFltFFT, omega_0, width, type); }
  else if (type==4) { filterWaveRC(white,realFltFFT,imgFltFFT,1.0/omega_0); }
  else if (type==5) { filterWaveCR(white,realFltFFT,imgFltFFT,1.0/omega_0); }
  double *ampl_tf = new double [npt];
  double *phi_tf = new double [npt];
  for (int ip=0; ip<npt; ++ip) {
		ampl_tf[ip]=TMath::Sqrt(realFltFFT[ip]*realFltFFT[ip]+imgFltFFT[ip]*imgFltFFT[ip]);
    phi_tf[ip]=(imgFltFFT[ip]||realFltFFT[ip])?TMath::ATan2(imgFltFFT[ip],realFltFFT[ip]):0.0;           
  }	
  TCanvas *c1 = new TCanvas("transfer function_ampl","TF");		
  c1->Divide(1,2);
  c1->SetLogx();
  c1->cd(1)->SetLogy();
  TGraph *function = new TGraph(npt, &white.omega[0], ampl_tf);	
  function->GetXaxis()->SetTitle("#omega [Hz]");
  function->SetTitle("amplitude");
  function->GetYaxis()->SetTitleOffset(1.4);
  function->GetYaxis()->SetTitle("Amplitude");  
  function->Draw("AL");
  c1->cd(2);	
  TGraph *function_phi = new TGraph(npt, &white.omega[0], phi_tf);	
  function_phi->GetXaxis()->SetTitle("#omega [Hz]");
  function_phi->SetTitle("phi");
  function_phi->GetYaxis()->SetTitleOffset(1.4);
  function_phi->GetYaxis()->SetTitle("phi");  
  function_phi->Draw("AL");

}
