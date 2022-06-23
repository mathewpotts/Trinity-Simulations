#include <iostream>

#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TAxis.h>
#include <TMath.h>
#include <TVectorD.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TString.h>
#include <TProfile.h>
#include <TSystem.h>

using namespace std;

/*
 Global variable for the "Processed_CORSIKA_ROOT_Files" directory
 using the same directory structure that appears in the shared directory
*/
Bool_t cd = gSystem->cd("../"); // must declare type
TString CWD = gSystem->pwd(); // One dir up from this file

/*
 Read in Trinity Shower information from root files and generate a 
 timing_data_histogram_shower%d.root file 
*/
void f(int showerid) // currently sums only one shower tile file at different elevations
{
  int nfiles = 11;
  int nshowers = 9; // there are 9 showers simulated from CORSIKA
  int maxazi = 40;
  int binsperangle = 2;
  
  TString filepath;
  TString internal;
  TString title;
  TString n;
  
  TH1D *arrivaltimes[nfiles * maxazi * binsperangle];
  
  int c = 0;
  int shower = showerid;	
  
  //for(int shower = 0; shower < nshowers; shower++) 
  {
    for(float ele = 0.0; ele <= 1.05; ele += 0.1)
      {
	filepath.Clear();
	filepath.Form("/TrinityShowers/SHOWER%d/22147882_%d_DAT10500%d_ELV%.1f_TELEHEIGHT_2000_TARGET_50.root", shower, shower, shower, ele);
	internal.Clear();
	
	TFile *fo = new TFile(CWD + filepath, "READ");
	
	for(float f = 0.25; f <= 39.751; f += 0.5)
	  {
	    n.Form("h_Arrivaltimes_AzimuthOffset%.2fdeg", f);
	    internal.Form("hist_shr%d_ele%.1f_azi%.2f", shower, ele, f);
	    title.Form("Photon arrival time distribution for shower %d, elevation angle %.1f, azimuth angle %.1f", shower, ele, f);
	    
	    arrivaltimes[c] = (TProfile*)fo->Get(n);
	    arrivaltimes[c]->SetNameTitle(internal, title);
	    c++;
	    
	    cout << c << endl;
	  }
      }
  }
  
  TString filename;
  filename.Form("/processed_root_files/timing_data_histogram_shower%d.root", showerid);
  TFile *fs = new TFile(CWD + filename, "RECREATE");  
  
  for(int i = 0; i < c; i++)
    {
      arrivaltimes[i]->Write();
    }
  
  fs->Close();
}

/*
 Go through all the timing_data_histogram_shower%d.root files to create a
 summed.root file containing all data less than or equal to 1.01 elevation 
 and 39.751 azimuth
 */
void g()
{
  TString filepaths[10];
  TFile *files[10];
  
  for(int i = 0; i < 10; i++)
    {
      filepaths[i].Form("/processed_root_files/timing_data_histogram_shower%d.root", i); 
      files[i] = new TFile(CWD + filepaths[i], "READ");
    }
  
  TH1D *h[880];
  TString title;
  int c = 0;
  
  for(float e = 0.0; e <= 1.01; e += 0.1)
    {
      for(float a = 0.25; a <= 39.751; a += 0.5)
	{
	  title.Form("hist_ele%.1f_azi%.2f", e, a);
	  h[c] = new TH1D(title, title, 1000000, 0, 1e6);
	  c++;
	}
    }
  
  TString internal;
  
  for(int i = 0; i < 10; i++)
    {
      c = 0;
      
      for(float e = 0.0; e <= 1.01; e += 0.1)
	{
	  for(float a = 0.25; a <= 39.751; a += 0.5)
	    {
	      internal.Form("hist_shr%d_ele%.1f_azi%.2f", i, e, a);
	      h[c]->Add((TH1D*) files[i]->Get(internal));
	      cout << internal << endl;
	      c++;
	    }
	}
    }
  
  TString filename;
  filename.Form("/processed_root_files/summed.root");
  
  TFile *fo = new TFile(CWD + filename, "RECREATE");
  
  for(int i = 0; i < 880; i++)
    {
      h[i]->Write();
    }
}

/*
 Average the summed.root file over the 10 showers and create a
 summed_avg.root file 
 */
void h()
{
  TString filename;
  filename.Form("/processed_root_files/summed.root");
  TFile *f = new TFile(CWD + filename, "READ");
  TH1D *h[880];
  
  TString title;
  int c = 0;
  
  for(float e = 0.0; e <= 1.01; e += 0.1)
    {
      for(float a = 0.25; a <= 39.751; a += 0.5)
	{
	  title.Form("hist_ele%.1f_azi%.2f", e, a);
	  //h[c] = new TH1D(title, title, 1000000, 0, 1e6);
	  h[c] = (TH1D*) f->Get(title);
	  
	  for(int bins = 1; bins <= h[c]->GetNbinsX(); bins++)
	    {
	      h[c]->SetBinContent(bins, h[c]->GetBinContent(bins) / 10.);
	    }
	  cout << title << endl;
	  c++;
	}
    }
  
  filename.Form("/processed_root_files/summed_avg.root");
  TFile *fo = new TFile(CWD + filename, "RECREATE");
  
  for(int i = 0; i < 880; i++)
    {
      h[i]->Write();
    }
}

/*
 Read information from the summed_avg.root and summed_timing_data_tprofile.root files and 
 print out to terminal
 
 Probably for debugging purposes
 */
void i()
{
  TString filename;
  filename.Form("/processed_root_files/summed_avg.root");
  TFile *histos = new TFile(CWD + filename, "READ");
  
  filename.Form("/processed_root_files/summed_timing_data_tprofiles.root"); // where are this file generated??
  TFile *profiles = new TFile(CWD + filename, "READ");
  
  //TProfile *p;
  //TH1D *h;
  
  TString h_internal;
  TString p_internal;
  int c = 1;
  
  for(float e = 0.0; e <= 1.01; e += 0.1)
    {
      p_internal.Form("ele_%.1f_deg", e);
      TProfile *p = (TProfile*) profiles->Get(p_internal);
      p->SetErrorOption("S");
      for(float a = 0.25; a <= 39.751; a += 0.5)
	{
	  h_internal.Form("hist_ele%.1f_azi%.2f", e, a);
	  TH1D *h = (TH1D*) histos->Get(h_internal);
	  
	  printf("hist entries: %.0f, hist avg: %.2f, hist sd: %.2f | prof entries: %.0f, prof avg: %.2f, prof sd: %.2f\n", 
		 h->GetEntries(), h->GetMean(), h->GetStdDev(), p->GetBinEntries(c), p->GetBinContent(c), p->GetBinError(c));
	  
	  c++;
	  delete(h);
	}
      delete(p);
    }
}

/*
 Read the summed_avg.root file to extract the Weighted Photon Arrival Timing Spread and create
 a 2D histogram 2dhisto.root

 Produces weird error where the file isn't properly closed:
 Warning in <TFile::Init>: file /storage/hive/project/phy-otte/mpotts32/trinity/trinity-simulations/Processed_CORSIKA_ROOT_Files/processed_root_files/summed_avg.root probably not closed, trying to recover
 */
void j()
{
  TString filename;
  filename.Form("/processed_root_files/summed_avg.root");
  cout << filename << endl;
  TFile *histos = new TFile(CWD + filename, "READ");
  //TH1D *h;
  TString h_internal;
  
  Double_t shower_e = 1e6; // shower energy
  Double_t t_trigger = 30; // trigger pulse width
  
  TH2D *h_w = new TH2D("timings_w", "Weighted Photon Arrival Timing Spread", 80, 0, 40, 10, 0, 1);
  
  for(float e = 0.0; e <= 1.01; e += 0.1)
    {
      for(float a = 0.25; a <= 39.751; a += 0.5)
	{
	  h_internal.Form("hist_ele%.1f_azi%.2f", e, a);
	  TH1D *h = (TH1D*) histos->Get(h_internal);
	  
	  Double_t rms_entries = 0;
	  Double_t w = 0;
	  
	  cout << "shower " << h_internal;
	  printf(": num particles in rms window of %.2f: %.0f\n", h->GetStdDev(), rms_entries);
	  
	  if(h->GetStdDev() < t_trigger)
	    {
	      int lowerbin = h->FindBin(h->GetMean() - h->GetStdDev() * 0.5);
	      int upperbin = h->FindBin(h->GetMean() + h->GetStdDev() * 0.5);
	      for(int i = lowerbin; i <= upperbin; i++)
		{
		  rms_entries += h->GetBinContent(i);
		}
	      w = rms_entries / (h->GetStdDev() * shower_e);
	      printf("rms windows smaller than trigger pulsewidth\n");
	      h_w->Fill(a, e, w);
	      printf("%f\n\n", w);
	    }
	  else
	    {
	      int lowerbin = h->FindBin(h->GetMean() - t_trigger * 0.5);
	      int upperbin = h->FindBin(h->GetMean() + t_trigger * 0.5);
	      for(int i = lowerbin; i <= upperbin; i++)
		{
		  rms_entries += h->GetBinContent(i);
		}
	      w = rms_entries / (t_trigger * shower_e);
	      h_w->Fill(a, e, w);
	      printf("%f\n\n", w);
	    }
	  //printf(": avg: %.2f, sd: %.2f, bin closest to avg: %d, right side: %d, left side: %d\n", h->GetMean(), h->GetStdDev(), h->FindBin(h->GetMean()), h->FindBin(h->GetMean() + h->GetStdDev() * 0.5), h->FindBin(h->GetMean() - h->GetStdDev() * 0.5));
	  delete(h);
	}
    }
  
  filename.Form("/processed_root_files/2dhisto.root");
  TFile *f = new TFile(CWD + filename, "RECREATE");
  h_w->GetXaxis()->SetTitle("Azimuth (degrees)");
  h_w->GetYaxis()->SetTitle("Elevation (degrees)");
  h_w->Write();
  
  f->Close();
  histos->Close();
}

int main(int argc, char **argv)
{
  // Generate timing_data_histogram_shower%d.root files
  // to add if statement eventually...
  //for (int i = 0; i < 10 ; i++)
  //  {
  //    cout << "Processing ShowerID: " << i << endl;
  //    f(i);
  //  }
  
  cout << "g()" << endl;
  g();
  
  cout << "h()" << endl;
  h();
  
  //i();
  
  cout << "j()" << endl;
  j();
  
  return 0;
}
