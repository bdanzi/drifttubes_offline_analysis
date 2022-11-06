//#include <Riostream.h>
#include <iostream>
#include <TF1.h>
#include <TMath.h>

Float_t * IsignalFromCharges (Float_t *ChargesSequence, Float_t *timeSequence, Int_t Nel, Int_t &outORCharge, Float_t ranget=2.0e-6, Float_t Rwire= 0.00125, Float_t Rtube= 2., Float_t Volt = 2000.0, Float_t mobility = 10.4e-6){

  Float_t timeres=800.0e-12;
  Int_t Npoint = TMath::Nint(2.0*ranget/timeres);
  Float_t *sign= new Float_t[Npoint];
  Float_t logRat = TMath::Log(Rtube/Rwire);
  Float_t to = pow(Rwire,2)*logRat/(2.*mobility*Volt);
  Float_t tmax = (pow(Rtube,2)-pow(Rwire,2))*logRat/(2.*mobility*Volt);
  to*=.000001;  //from micros to s
  tmax*=.000001;  //from micros to s

  TF1 *peakshape = new TF1("peakshape","[0]/([1]*([2]+x))",0.,tmax);
  peakshape->SetParameter(1,2.0*logRat);
  peakshape->SetParameter(2,to);


  for(Int_t i=0; i<Npoint; i++) {
    sign[i]=0.;
  }
  Int_t binCharge = 0;
  Float_t signTime;
  Int_t CutNpoint = (Int_t) (0.5*Npoint-1);
  Int_t jt;
  outORCharge = 0;

  for (Int_t charge=0; charge<Nel; charge++){
    binCharge = (Int_t) ((timeSequence[charge]*1.e-6)/timeres);
    if (binCharge>=CutNpoint){
      outORCharge++;
      std::cout<<timeSequence[charge]<<" "<<binCharge<<" ----- "<<Npoint <<" ------\n";
      continue;
    }
    peakshape->SetParameter(0,ChargesSequence[charge]*TMath::Qe());
    jt=0;
    for (Int_t i=binCharge; i<Npoint; i++) {
      signTime=(((Float_t)jt)+.5)*timeres;
      sign[i]+=peakshape->Eval(signTime)*1000.0;//im mA //1000000.0;  //in microA
      jt++;
    }
  }
  return sign;
}
