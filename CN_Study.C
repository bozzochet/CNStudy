TTree* OpenTree(float CNs[2][9][16], float Signals[2][9][1024]) {
  //////////////////////////////////////////////////////////
  //   This file has been automatically generated 
    //     (Fri Dec  1 15:03:30 2023 by ROOT version6.26/06)
    //   from TTree t3/My calibration tree
    //   found on file: run_010855_ONLYCAL.root
    //////////////////////////////////////////////////////////
    
    TTree* t3 = NULL;
  
  //Reset ROOT and connect tree file
  gROOT->Reset();
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("run_010855_ONLYCAL.root");
  if (!f) {
    f = new TFile("run_010855_ONLYCAL.root");
  }
  f->GetObject("t3",t3);
    
  // Set branch addresses.
  t3->SetBranchAddress("CNs",CNs);
  t3->SetBranchAddress("Signals",Signals);
  
  //     This is the loop skeleton
  //       To read only selected branches, Insert statements like:
  // t3->SetBranchStatus("*",0);  // disable all branches
  // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname

  return t3;
}

void shiftwholeleft(float* array){

  //  printf("Shift whole left\n");

  for (int ii=0; ii<63; ii++) {
    array[ii] = array[ii+1];
  }
  
  return;
}

void shiftwholeright(float* array){

  //  printf("Shift whole right\n");

  for (int ii=63; ii>0; ii--) {
    array[ii] = array[ii-1];
  }
  
  return;
}

void shifthalfleft(float* array){

  for (int ii=0; ii<31; ii++) {
    array[ii] = array[ii+1];
  }
  
  return;
}

void shifthalfright(float* array){

  for (int ii=63; ii>32; ii--) {
    array[ii] = array[ii-1];
  }
  
  return;
}

void oddswap(float* array){

  for (int ii=0; ii<63; ii+=2) {
    if (array[ii]>array[ii+1]) {
      float saved = array[ii];
      array[ii] = array[ii+1];
      array[ii+1] = saved;
    }
  }
  
  return;
}

void evenswap(float* array){
  
  for (int ii=1; ii<62; ii+=2) {
    if (array[ii]>array[ii+1]) {
      float saved = array[ii];
      array[ii] = array[ii+1];
      array[ii+1] = saved;
    }
  }
  
  return;
}

void print_array(TString middle, float* array){
  
  printf("###############################%s###################################\n", middle.Data());
  for (int ii=0; ii<8; ii++){
    for (int jj=0; jj<8; jj++) {
      printf("a[%02d] = %f, ", 8*ii+jj, array[8*ii+jj]);
    }
    printf("\n");
  }
  
  return;
}

float MedianAlgoSimulation(float signal, bool restart, bool verbose=false){

  static float array[64];
  static int stepsdone = 0;
  if (restart) {
    stepsdone=0;
    for (int ii=0; ii<32; ii++) {
      array[ii] = -(pow(2, 14)-1-ii);
      array[ii+32] = (pow(2, 14)-1-ii);
    }
  }
  
  //insertion
  if (stepsdone%2) {
    shifthalfleft(array);
    if (verbose) print_array(" A ", array);
    array[31] = signal;
    if (verbose) print_array(" B ", array);
    shiftwholeright(array);
  }
  else {
    shifthalfright(array);
    if (verbose) print_array(" A ", array);
    array[32] = signal;
    if (verbose) print_array(" B ", array);
    shiftwholeleft(array);
  }
  stepsdone++;
  
  // it important to do even before odd:
  // - [30] with [31]
  // - [32] with [33]
  // and only afterwords [31] with [32]
  // that for median (the mean of these two) is even useless
  //even swap
  evenswap(array);

  //odd swap
  oddswap(array);
  
  if (verbose) print_array("###", array);
  float median = 0.5*(array[31] + array[32]);
  //  printf("median = %f\n", median);
  
  return median;
}

void CN_Toy(){

  double mean_ped = 1900;
  double var_ped = 100;
  double mean_sigma = 7.87;
  double var_CN = 10.0;
  int ADC_float = 4; //floating numbers are in units of 1/ADC_float ADC

  TRandom* rand = new TRandom3();

  //----------------------------------------------

  TH1F* i_ped = new TH1F("i_ped", "i_ped; ch; Entries", 1024, 0, 1024);
  TH1F* ped = new TH1F("ped", "ped; ch; ADC", 1024, 0, 1024);
  ped->SetLineColor(kRed+2);
  ped->SetMarkerColor(kRed+2);
  TH1F* delta_ped = new TH1F("delta_ped", "delta_ped; ch; #Delta_{ped} (%)", 1024, 0, 1024);
  delta_ped->SetLineColor(kRed+2);
  delta_ped->SetMarkerColor(kRed+2);
  TH1F* ped_s = new TH1F("ped_s", "ped_s; ch; ADC", 1024, 0, 1024);
  ped_s->SetLineColor(kBlue+1);
  ped_s->SetMarkerColor(kBlue+1);
  TH1F* delta_ped_s = new TH1F("delta_ped_s", "delta_ped_s; ch;  #Delta_{ped} (%)", 1024, 0, 1024);
  delta_ped_s->SetLineColor(kBlue+1);
  delta_ped_s->SetMarkerColor(kBlue+1);
  TH1F* i_sig = new TH1F("i_sig", "i_sig; ch; Entries", 1024, 0, 1024);
  TH1F* sig = new TH1F("sig", "sig; ch; ADC", 1024, 0, 1024);
  sig->SetLineColor(kRed+2);
  sig->SetMarkerColor(kRed+2);
  TH1F* delta_sig = new TH1F("delta_sig", "delta_sig; ch; #Delta_{CN} (%)", 1024, 0, 1024);
  delta_sig->SetLineColor(kRed+2);
  delta_sig->SetMarkerColor(kRed+2);
  TH1F* rsig = new TH1F("rsig", "rsig; ch; ADC", 1024, 0, 1024);
  rsig->SetLineColor(kRed+2);
  rsig->SetMarkerColor(kRed+2);

  //----------------------------------------------

  TH1F* CN = new TH1F("CN", "CN; CN (ADC); Entries", 1000, -10*var_CN, 10*var_CN);
  TH1F* CN_SW = new TH1F("CN_SW", "CN_SW; CN SW median (ADC); Entries", 1000, -10*var_CN, 10*var_CN);
  CN_SW->SetLineColor(kBlue+1);
  CN_SW->SetMarkerColor(kBlue+1);
  TH1F* CN_SW_mean = new TH1F("CN_SW_mean", "CN_SW_mean; CN SW mean (ADC); Entries", 1000, -10*var_CN, 10*var_CN);
  CN_SW_mean->SetLineColor(kGreen+2);
  CN_SW_mean->SetMarkerColor(kGreen+2);
  TH1F* CN_FPGA = new TH1F("CN_FPGA", "CN_FPGA; CN FPGA median (ADC); Entries", 1000, -10*var_CN, 10*var_CN);
  CN_FPGA->SetLineColor(kRed+2);
  CN_FPGA->SetMarkerColor(kRed+2);
  TH1F* diff = new TH1F("diff", "diff; (SW median - injected) / injected (%); Entries", 1000, -500.0/var_CN, 500.0/var_CN);
  TH1F* diff_mean = new TH1F("diff_mean", "diff_meas; (SW mean - injected) / injected (%); Entries", 1000, -500.0/var_CN, 500.0/var_CN);
  diff_mean->SetLineColor(kGreen+2);
  diff_mean->SetMarkerColor(kGreen+2);
  TH1F* diff_meas = new TH1F("diff_meas", "diff_meas; (SW median - FPGA median) / SW median (%); Entries", 1000, -500.0/var_CN, 500.0/var_CN);
  diff_meas->SetLineColor(kRed+2);
  diff_meas->SetMarkerColor(kRed+2);
  
  //----------------------------------------------

  float i_pedestal[1024];
  std::vector<float> vsignal[1024];
  std::vector<float> vCN[16];
  long int nentries = 100000;
  for (long int cc = 0; cc < 1024; cc++) {
    i_pedestal[cc] = rand->Gaus(mean_ped, var_ped);
    i_ped->SetBinContent(cc+1, i_pedestal[cc]);
    i_ped->SetBinError(cc+1, 0.0);
    i_sig->SetBinContent(cc+1, mean_sigma);
    i_sig->SetBinError(cc+1, 0.0);
  }
  for (Long64_t i=0; i<nentries; i++) {
    for (long int vv=0; vv<16; vv++){
      double thisCN = rand->Gaus(0, var_CN);
      CN->Fill(thisCN);
      vCN[vv].push_back(thisCN);
      for (long int cc = 64*vv; cc < 64*(vv+1); cc++) {
	double sig = rand->Gaus(i_pedestal[cc], mean_sigma);
	vsignal[cc].push_back(std::round(sig+thisCN));
      }
    }
  }
  
  //----------------------------------------------

  float pedestal[1024];
  float sigmaraw[1024];
  for (long int cc = 0; cc < 1024; cc++) {
    auto beginItr = std::begin(vsignal[cc]);
    auto endItr = std::end(vsignal[cc]);
    
    auto nEv = std::distance(beginItr, endItr);
    
    pedestal[cc] = std::accumulate(beginItr, endItr, 0.0f)/float(nEv);
    pedestal[cc] = std::floor(ADC_float*pedestal[cc])/ADC_float;
    ped->SetBinContent(cc+1, pedestal[cc]);
    ped->SetBinError(cc+1, 0.0);
    if (fabs((pedestal[cc] -i_pedestal[cc])/i_pedestal[cc])>1.0E-3) {
      printf("Pedestal[%lu] = %f, %f\n", cc, pedestal[cc], i_pedestal[cc]);
    }
    delta_ped->SetBinContent(cc+1, 100.0*(pedestal[cc]-i_pedestal[cc])/i_pedestal[cc]);
    delta_ped->SetBinError(cc+1, 0.0);    
    
    sigmaraw[cc] =
      std::sqrt(std::accumulate(beginItr, endItr, 0.0f,
				[&](float acc, float curr) {
				  return acc + (curr - pedestal[cc]) *
				    (curr - pedestal[cc]);
				}) / float(nEv));
    //    printf("Sigma_raw[%lu] = %f\n", cc, sigmaraw[cc]);
    sigmaraw[cc] = std::floor(ADC_float*sigmaraw[cc])/ADC_float;
    rsig->SetBinContent(cc+1, sigmaraw[cc]);
    rsig->SetBinError(cc+1, 0.0);
  }
  
  //----------------------------------------------

  std::vector<float> vsignal_s[1024];
  
  for (Long64_t i=0; i<nentries;i++) {
    //  for (Long64_t i=0; i<1;i++) {
    
    for (long int vv=0; vv<16; vv++){
      std::vector<float> vmedian;

      float mean = 0;
      float median = 0;
      float FPGAmedian = 0;
      
      for (long int cc = 64*vv; cc < 64*(vv+1); cc++) {
	vmedian.push_back(vsignal[cc].at(i)-pedestal[cc]);
	
	//--------------------------------------------------------------------------
	bool restart =  false;
	if (cc==64*vv) {//first time
	  restart = true;
	}
	//      printf("Channel = %lu\n", cc);
	//      FPGAmedian = MedianAlgoSimulation(cc, restart, true); //FOR TEST
	FPGAmedian = MedianAlgoSimulation(vsignal[cc].at(i)-pedestal[cc], restart);//at each call the median is returned
      }
      FPGAmedian = std::floor(ADC_float*FPGAmedian)/ADC_float;

      mean = std::accumulate(vmedian.begin(), vmedian.end(), 0.0f)/float(64);
      mean = std::floor(ADC_float*mean)/ADC_float;
      
      std::sort(begin(vmedian), end(vmedian));
      median = 0.5 * (vmedian[(vmedian.size() / 2) - 1] + vmedian[vmedian.size() / 2]);
      median = std::floor(ADC_float*median)/ADC_float;
      
      //--------------------------------------------------------------------------
      diff->Fill(100.0*(median - vCN[vv].at(i))/vCN[vv].at(i));
      diff_mean->Fill(100.0*(mean - vCN[vv].at(i))/vCN[vv].at(i));
      diff_meas->Fill(100.0*(median - FPGAmedian)/median);
      CN_SW->Fill(median);
      CN_SW_mean->Fill(mean);
      CN_FPGA->Fill(FPGAmedian);
      
      for (long int cc = 64*vv; cc < 64*(vv+1); cc++) {
	float to_insert = vsignal[cc].at(i)-FPGAmedian;
	to_insert = std::floor(ADC_float*to_insert)/ADC_float;
	vsignal_s[cc].push_back(to_insert);
      }

    }
    
  }

  //----------------------------------------------

  float pedestal_s[1024];
  float sigma[1024];
  
  for (long int cc = 0; cc < 1024; cc++) {
    auto beginItr = std::begin(vsignal_s[cc]);
    auto endItr = std::end(vsignal_s[cc]);
    
    auto nEv = std::distance(beginItr, endItr);
    
    pedestal_s[cc] = std::accumulate(beginItr, endItr, 0.0f)/float(nEv);
    pedestal_s[cc] = std::floor(ADC_float*pedestal_s[cc])/ADC_float;
    ped_s->SetBinContent(cc+1, pedestal_s[cc]);
    ped_s->SetBinError(cc+1, 0.0);    
    if (fabs((pedestal_s[cc] -i_pedestal[cc])/i_pedestal[cc])>1.0E-3) {
      printf("Pedestal_s[%lu] = %f, %f\n", cc, pedestal_s[cc], i_pedestal[cc]);
    }
    delta_ped_s->SetBinContent(cc+1, 100.0*(pedestal_s[cc]-i_pedestal[cc])/i_pedestal[cc]);
    delta_ped_s->SetBinError(cc+1, 0.0);    
    
    sigma[cc] =
      std::sqrt(std::accumulate(beginItr, endItr, 0.0f,
				[&](float acc, float curr) {
				  return acc + (curr - pedestal_s[cc]) *
				    (curr - pedestal_s[cc]);
				}) / float(nEv));
    //    printf("Sigma[%lu] = %f\n", cc, sigma[cc]);
    sigma[cc] = std::floor(ADC_float*sigma[cc])/ADC_float;
    sig->SetBinContent(cc+1, sigma[cc]);
    sig->SetBinError(cc+1, 0.0);
    delta_sig->SetBinContent(cc+1, 100.0*(sigma[cc]-mean_sigma)/mean_sigma);
    delta_sig->SetBinError(cc+1, 0.0);
  }
  
  //----------------------------------------------
  
  TCanvas* c0 = new TCanvas();
  c0->cd();
  c0->Divide(2,2);
  c0->cd(1);
  i_ped->Draw("hist");
  ped->Draw("samehist");
  ped_s->Draw("samehist");
  c0->cd(2);
  rsig->Draw("hist");
  c0->cd(3);
  i_sig->Draw("hist");
  sig->Draw("samehist");

  TCanvas* c0_bis = new TCanvas();
  c0_bis->cd();
  c0_bis->Divide(2,2);
  c0_bis->cd(1);
  delta_ped->Draw("samehist");
  delta_ped_s->Draw("samehist");
  c0_bis->cd(3);
  delta_sig->Draw("hist");
  
  TCanvas* c1 = new TCanvas();
  c1->cd();
  diff_mean->Draw();
  diff_meas->Draw("same");
  diff->Draw("same");


  TCanvas* c2 = new TCanvas();
  c2->cd();
  CN->Draw();
  CN_SW->Draw("same");
  CN_SW_mean->Draw("same");
  CN_FPGA->Draw("same");
  
  return;
}

void CN_Verification(){
  
  //Declaration of leaves types
  float        CNs[2][9][16];
  float        Signals[2][9][1024];

  TTree* t3 = OpenTree(CNs, Signals);

  Long64_t nentries = t3->GetEntries();

  //----------------------------------------------
  
  float pedestal[1024];
  float sigmaraw[1024];
  std::vector<float> vsignal[1024];
  Long64_t nbytes = 0;
  for (Long64_t i=0; i<nentries;i++) {
    nbytes += t3->GetEntry(i);
    
    long int jj=0;//LINF number
    long int tt=0;//LEF number
    
    for (long int cc = 0; cc < 1024; cc++) {
      if (Signals[jj][tt][cc] > 0.0) {
	vsignal[cc].push_back(Signals[jj][tt][cc]);
      }
    }
  }

  //----------------------------------------------
  
  for (long int cc = 0; cc < 1024; cc++) {
    auto beginItr = std::begin(vsignal[cc]);
    auto endItr = std::end(vsignal[cc]);
    
    auto nEv = std::distance(beginItr, endItr);

    pedestal[cc] = std::accumulate(beginItr, endItr, 0.0f)/float(nEv);
    //    printf("Pedestal[%lu] = %f\n", cc, pedestal[cc]);
    
    sigmaraw[cc] =
      std::sqrt(std::accumulate(beginItr, endItr, 0.0f,
				[&](float acc, float curr) {
				  return acc + (curr - pedestal[cc]) *
				    (curr - pedestal[cc]);
				}) / float(nEv));
    //    printf("Sigma_raw[%lu] = %f\n", cc, sigmaraw[cc]);
  }

  //----------------------------------------------

  TH1F* diff = new TH1F("diff", "diff; SW median - FPGA median (ADC); Entries", 1000, -0.1, 0.1); 
  
  //----------------------------------------------
  
  nbytes = 0;
  for (Long64_t i=0; i<nentries;i++) {
    //  for (Long64_t i=0; i<1;i++) {
    nbytes += t3->GetEntry(i);
    
    long int jj=0;//LINF number
    long int tt=0;//LEF number
    long int vv=2;//VA number
    
    std::vector<float> vmedian;
    float CN = CNs[jj][tt][vv];

    float median = 0;
    float FPGAmedian = 0;
    
    for (long int cc = 64*vv; cc < 64*(vv+1); cc++) {
      vmedian.push_back(Signals[jj][tt][cc]-pedestal[cc]);
      
      //--------------------------------------------------------------------------
      bool restart =  false;
      if (cc==64*vv) {//first time
	restart = true;
      }
      //      printf("Channel = %lu\n", cc);
      //      FPGAmedian = MedianAlgoSimulation(cc, restart, true); //FOR TEST
      FPGAmedian = MedianAlgoSimulation(Signals[jj][tt][cc]-pedestal[cc], restart);//at each call the median is returned
    }
    
    std::sort(begin(vmedian), end(vmedian));
    median = 0.5 * (vmedian[(vmedian.size() / 2) - 1] + vmedian[vmedian.size() / 2]);

    //--------------------------------------------------------------------------
    //    printf("CN = %f, median = %f, diff = %E\n", CN, median, CN-median);//for some reason the diffs are 10E-2 - 10E-3, even if the code to compute seems really the same...
    //    printf("CN = %f, FPGAmedian = %f, diff = %E\n", CN, FPGAmedian, CN-FPGAmedian);
    diff->Fill(median - FPGAmedian);
  }

  diff->Draw();
  
  return;
}

void CN_Plot(){

  TFile* foutput = new TFile("foutput.root", "RECREATE");
  foutput->cd();
  
  TF1* fgaus = new TF1("fgaus", "gaus(0)", -10, 10);
  
  //Declaration of leaves types
  float        CNs[2][9][16];
  float        Signals[2][9][1024];

  TTree* t3 = OpenTree(CNs, Signals);

  //--------------------------------------------------------------------------
  
  TH1D* hCN[2][9][16];
  for (long int jj=0; jj<2; jj++) {
    for (long int tt=0; tt<9; tt++) {
      for (long int vv=0; vv<16; vv++) {
	hCN[jj][tt][vv] = new TH1D(Form("hCN_%lu_%lu_%lu", jj, tt, vv), Form("VA%lu; CN (ADC); Entries", vv), 201, -10, 10);
	if (vv%2 == 0 ) {
	  hCN[jj][tt][vv]->SetLineColor(kRed+1);
	  hCN[jj][tt][vv]->SetMarkerColor(kRed+1);
	}
	else {
	  hCN[jj][tt][vv]->SetLineColor(kBlue+1);
	  hCN[jj][tt][vv]->SetMarkerColor(kBlue+1);
	}
	hCN[jj][tt][vv]->SetLineWidth(3);
      }
    }
  }
  
  TH1D* hADC[2][9][8];
  for (long int jj=0; jj<2; jj++) {
    for (long int tt=0; tt<9; tt++) {
      for (long int aa=0; aa<8; aa++) {
	hADC[jj][tt][aa] = new TH1D(Form("hADC_%lu_%lu_%lu", jj, tt, aa), Form("ADC%lu; CN (ADC); Entries", aa), 201, -10, 10);
	hADC[jj][tt][aa]->SetLineColor(kGreen+2);
	hADC[jj][tt][aa]->SetMarkerColor(kGreen+2);
	hADC[jj][tt][aa]->SetLineWidth(3);
      }
    }
  }

  TH1D* hCNdiff[2][9][16][16];
  TH2D* hCNcorr[2][9][16][16];
  for (long int jj=0; jj<2; jj++) {
    for (long int tt=0; tt<9; tt++) {
      for (long int vv=0; vv<16; vv++) {
	for (long int vv2=0; vv2<16; vv2++) {
	  hCNdiff[jj][tt][vv][vv2] = new TH1D(Form("hCNdiff_%lu_%lu_%lu_%lu", jj, tt, vv, vv2), Form("VA%lu-VA%lu; CN (ADC); Entries", vv, vv2), 201, -10, 10);
	  hCNdiff[jj][tt][vv][vv2]->SetLineColor(kOrange);
	  hCNdiff[jj][tt][vv][vv2]->SetMarkerColor(kOrange);
	  hCNdiff[jj][tt][vv][vv2]->SetLineWidth(3);
	  
	  hCNcorr[jj][tt][vv][vv2] = new TH2D(Form("hCNcorr_%lu_%lu_%lu_%lu", jj, tt, vv, vv2), Form("VA%lu_vs_VA%lu; VA%lu CN (ADC); VA%lu CN (ADC)", vv, vv2, vv, vv2), 201, -10, 10, 201, -10, 10);
	  hCNcorr[jj][tt][vv][vv2]->SetMarkerColor(kOrange);
	  hCNcorr[jj][tt][vv][vv2]->SetMarkerSize(1.0);
	}
      }
    }
  }
  
  TH1D* hSENS[2][9];
  TH1D* hSENS_median[2][9];
  for (long int jj=0; jj<2; jj++) {
    for (long int tt=0; tt<9; tt++) {
      hSENS[jj][tt] = new TH1D(Form("hSENS_%lu_%lu", jj, tt), Form("SENS%lu; CN (ADC); Entries", tt), 201, -10, 10);
      hSENS[jj][tt]->SetLineColor(kBlack);
      hSENS[jj][tt]->SetMarkerColor(kBlack);
      hSENS[jj][tt]->SetLineWidth(3);
      
      hSENS_median[jj][tt] = new TH1D(Form("hSENS_median_%lu_%lu", jj, tt), Form("SENS%lu_median; CN (ADC); Entries", tt), 201, -10, 10);
      hSENS_median[jj][tt]->SetLineColor(kViolet+1);
      hSENS_median[jj][tt]->SetMarkerColor(kViolet+1);
      hSENS_median[jj][tt]->SetLineWidth(3);
    }
  }
  
  Long64_t nentries = t3->GetEntries();
  
  Long64_t nbytes = 0;
  for (Long64_t i=0; i<nentries;i++) {
    nbytes += t3->GetEntry(i);
    
    for (long int jj=0; jj<2; jj++) {
      for (long int tt=0; tt<9; tt++) {
	double _SENS = 0.0;
	std::vector<double> _SENS_median;
	for (long int vv = 0; vv < 16; vv++) {
	  _SENS += CNs[jj][tt][vv];
	  _SENS_median.push_back(CNs[jj][tt][vv]);
	  hCN[jj][tt][vv]->Fill(CNs[jj][tt][vv]);
	  if (vv%2==0) {
	    int aa = vv/2;//automatically rounded to even (anyhow here vv is even)
	    hADC[jj][tt][aa]->Fill(0.5*(CNs[jj][tt][vv]+CNs[jj][tt][vv+1]));
	  }
	  for (long int vv2 = 0; vv2 < 16; vv2++) {
	    if (vv<15 && vv2<15) {
	      hCNdiff[jj][tt][vv][vv2]->Fill(CNs[jj][tt][vv]-CNs[jj][tt][vv2]);
	      hCNcorr[jj][tt][vv][vv2]->Fill(CNs[jj][tt][vv], CNs[jj][tt][vv2]);
	    }
	  }
	  if (fabs(CNs[jj][tt][vv]) > 100) {
	    printf("CNs[%lu][%lu][%lu] = %f\n", jj, tt, vv, CNs[jj][tt][vv]);
	  }
	}
	hSENS[jj][tt]->Fill(_SENS/16);
	std::sort(begin(_SENS_median), end(_SENS_median));
	hSENS_median[jj][tt]->Fill(0.5 * (_SENS_median[(_SENS_median.size() / 2) - 1] + _SENS_median[_SENS_median.size() / 2]));
      }
    }
  }
  
  int lef_to_plot = 0;
  int VA_to_plot = 1;
  //  int VA2_to_plot = VA_to_plot+1;//for diff and correlation
  int VA2_to_plot = 13;
  int ADC_to_plot = VA_to_plot/2;//automatically rounded to even
  if (VA_to_plot%2) {
    printf("You required VA %d to be plotted. Also %d will be used.\n", VA_to_plot, VA_to_plot+1);
    printf("The VA could be on the edge and/or the two VAs are not on the same ADC (we'll plot %d).\n", ADC_to_plot);
  }
  int Fig_no = 5;
  
  if (Fig_no == 1) {
    hCN[0][lef_to_plot][VA_to_plot]->Draw();
    fgaus->SetLineColor(kRed+2);
    hCN[0][lef_to_plot][VA_to_plot]->Fit(fgaus);
  }
  else if (Fig_no == 2) {
    hCN[0][lef_to_plot][VA_to_plot]->Draw();
    hCN[0][lef_to_plot][VA_to_plot+1]->Draw("same");
  }
  else if (Fig_no == 3) {
    hADC[0][lef_to_plot][ADC_to_plot]->Draw();
    hCN[0][lef_to_plot][VA_to_plot]->Draw("same");
    hCN[0][lef_to_plot][VA_to_plot+1]->Draw("same");
    hADC[0][lef_to_plot][ADC_to_plot]->Draw("same");
    fgaus->SetLineColor(kGreen+2);
    hADC[0][lef_to_plot][ADC_to_plot]->Fit(fgaus);
    hADC[0][lef_to_plot][ADC_to_plot]->Fit(fgaus);
  }
  else if (Fig_no == 4) {
    hCNdiff[0][lef_to_plot][VA_to_plot][VA2_to_plot]->Draw("same");
  }
  else if (Fig_no == 5) {
    hCNcorr[0][lef_to_plot][VA_to_plot][VA2_to_plot]->Draw("colz");
  }
  else {
    hCN[0][lef_to_plot][VA_to_plot]->Draw();
    //   fgaus->SetLineColor(kRed+2);
    //   hCN[0][lef_to_plot][VA_to_plot]->Fit(fgaus);
    hCN[0][lef_to_plot][VA_to_plot+1]->Draw("same");
    hADC[0][lef_to_plot][ADC_to_plot]->Draw("same");
    //   hSENS[0][lef_to_plot]->Draw("same");
    //   hSENS_median[0][lef_to_plot]->Draw("same");
  }
  
  foutput->Write();
  //  foutput->Close();

  return;
}

void CN_Study() {

  CN_Toy();
  
  return;  
}
