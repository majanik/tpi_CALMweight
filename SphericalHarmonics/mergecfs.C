void mergecfs(const char *in1, const char *in2, const char *out, const char *pair="PipPp")
{
  TFile *inf1 = new TFile(in1);
  TFile *inf2 = new TFile(in2);
  TFile *outf = new TFile(out, "RECREATE");

  char fname[200];

  sprintf(fname, "NumReYlm%i%i%sNonIdYlms", 0, 0, pair);
  cout << "Reading " << fname << endl;
      
  TH1D *nm1 = inf1->Get(fname);
  TH1D *nm2 = inf2->Get(fname);

  for (int el=0; el<=3; el++) {
    for (int em=0; em<el*2+1; em++) {
      sprintf(fname, "CfnReYlm%i%i%sNonIdYlms", el, em, pair);
      cout << "Reading " << fname << endl;
      
      TH1D *cf1 = inf1->Get(fname);
      TH1D *cf2 = inf2->Get(fname);

      sprintf(fname, "CovCfc%sNonIdYlms", pair);
      TH3D *cv1 = inf1->Get(fname);
      TH3D *cv2 = inf2->Get(fname);


      TH1D *cfo = new TH1D(*cf1);
      for (int iter=1; iter<=cfo->GetNbinsX(); iter++) {

	double w1 = cf1->GetBinError(iter);
	if (w1>0.0) {
	  w1 *= w1;
	  w1 = 1/w1;
	}
	else 
	  w1 = 0;
	double w2 = cf2->GetBinError(iter);
	if (w2 > 0.0) {
	  w2 *= w2;
	  w2 = 1/w2;
	}
	else
	  w2 = 0;

	if ((w1>0.0) && (w2>0.0)) {
	  cfo->SetBinContent(iter, (cf1->GetBinContent(iter)*w1 + cf2->GetBinContent(iter)*w2)/(w1+w2));
	  cfo->SetBinError(iter, TMath::Sqrt(1.0/(w1+w2)));
	}
	else {
	  cfo->SetBinContent(iter,0.0);
	  cfo->SetBinError(iter, 0.0);
	}
      }
      outf->cd();
      cfo->Write();
      
      sprintf(fname, "CfnImYlm%i%i%sNonIdYlms", el, em, pair);
      cout << "Reading " << fname << endl;
      
      cf1 = (TH1D *) inf1->Get(fname);
      cf2 = (TH1D *) inf2->Get(fname);
      
      cfo = new TH1D(*cf1);
      for (int iter=1; iter<=cfo->GetNbinsX(); iter++) {

	double w1 = cf1->GetBinError(iter);
	if (w1>0.0) {
	  w1 *= w1;
	  w1 = 1/w1;
	}
	else 
	  w1 = 0;
	double w2 = cf2->GetBinError(iter);
	if (w2 > 0.0) {
	  w2 *= w2;
	  w2 = 1/w2;
	}
	else
	  w2 = 0;

	if ((w1>0.0) && (w2>0.0)) {
	  cfo->SetBinContent(iter, (cf1->GetBinContent(iter)*w1 + cf2->GetBinContent(iter)*w2)/(w1+w2));
	  cfo->SetBinError(iter, TMath::Sqrt(1.0/(w1+w2)));
	}
	else {
	  cfo->SetBinContent(iter,0.0);
	  cfo->SetBinError(iter, 0.0);
	}
      }
      outf->cd();
      cfo->Write();
      
    }
  }
  
}
