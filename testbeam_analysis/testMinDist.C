
void testMinDist(float sigma=1.0, float mean=0.0, int nTest=1000, int nBins=1000) {

	TH1F *hTest=new TH1F("hTest","Mean of Gaus Max/Min dist",200, 0, TMath::Nint(10*sigma));
	TRandom3 *rndm= new TRandom3();
	for (int itry=0; itry<nTest; ++itry) {
		rndm->SetSeed();
//		cout<<"itry "<<itry<<" seed "<<rndm->GetSeed()<<endl;
		float max=mean-1e6*sigma;
		for (int ibn=0; ibn<nBins; ++ibn) {
			float tmpval=rndm->Gaus(mean,sigma);
			if (tmpval>max) { max=tmpval; }
		}
		hTest->Fill(max-mean);
	}
	hTest->Draw();
}
