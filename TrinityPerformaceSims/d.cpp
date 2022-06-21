void d()
{
	TCanvas *c1 = new TCanvas("c1","Test Canvas",1600,750);
	c1->Divide(2,1);
	TH2F *skymapA=new TH2F("skymapA","Skymap Plot A",3601, -180.05, 180.05, 1801, -90.05, 90.05);
	TH2F *skymapB=new TH2F("skymapB","Skymap Plot B",360, -180, 180, 360, -90, 90);
	
	for(int i = 1; i <= 500; i++)
	{
		for(int j = 1; j <= 540; j++)
		{
			skymapA->SetBinContent(i, j, 20000.0);
		}
	}
	
	//~ for(int i = 1; i < 1082; i++){
		//~ skymapA->SetBinContent(1, i, 20000.0);
	//~ }
	
	for(int i = -10; i < 10; i++)
	{
		for(int j = -10; j < 10; j++)
		{
			skymapB->Fill(i, j, 10);
		}
	}
	
	//~ skymapA->Fill(0.,0., -20);
	skymapA->Fill(0., 0., 10);
	skymapA->Fill(0., 0., 10);
	c1->cd(1);
	//~ skymapA->Draw("z aitoff");
	c1->cd(2);
	skymapB->Fill(0.,0.,-1);
	gStyle->SetPalette(56);
	skymapB->Draw("z aitoff");
	cout<<""<<skymapA->GetBinContent(1801,901)<<" "<<(int)ceil(skymapA->GetNbinsX() / 2.0)<<" "<<(int)ceil(skymapA->GetNbinsY() / 2.0)<<endl; //1801, 901 = origin
	Double_t epic = (atan2(sin((228.73 - -50) * (3.141592 / 180.0)), cos((228.73 - -50) * (3.141592 / 180.0)) * sin(40.75639 * (3.141592 / 180.0)) - tan(-15 * (3.141592 / 180.0)) * cos(40.75639 * (3.141592 / 180.0))) * 180 / 3.141592) - 180;
	cout<<epic<<endl; 
	
}

void xd()
{
	Double_t degconv = (3.14159265353 / 180.0), pi = 3.14159265353, LST = 212.33, latitude = 38.52028;
	
	for(int r = -180; r <= 180; r++){
	for(int d = -90; d <= 90; d++){
	Double_t dec = asin(sin(d * degconv) * sin(27.1284 * degconv) + cos(d * degconv) * cos(27.1284 * degconv) * cos((122.9320 - r) * degconv)) * 180 / pi;
	Double_t ra = (atan2((cos(d * degconv) * sin((122.932 - r) * degconv)), (sin(d * degconv) * cos(27.1284 * degconv) - cos(d * degconv) * sin(27.1284 * degconv) * cos((122.9320 - r) * degconv))) * 180 / pi) + 192.8595;
	if(ra > 180)
		ra = ra - 360.;
	if(ra < -180)
		ra = ra + 360.;
	Double_t az = (atan2(sin((LST - ra) * degconv), cos((LST - ra) * degconv) * sin(latitude * degconv) - tan(dec * degconv) * cos(latitude * degconv)) * 180 / pi) - 180;
	Double_t alt = asin(sin(latitude * degconv) * sin(dec * degconv) + cos(latitude * degconv) * cos(dec * degconv) * cos((LST - ra) * degconv)) * 180 / pi;
	if(az > 180.0)
		az = az - 360.0;
	if(az < -180.0)
		az = az + 360.0;
	//~ if(az > 180.0 || az < -180.0 || alt > 90.0 || alt < -90.0)
	cout<<"azi: "<<(int)((az + 180.1) * 10)<<" elv: "<<(int)((alt + 90.1) * 10)<<endl;
	}
	}	
}

void quantiles() {
   // demo for quantiles
   const Int_t nq = 20;
   TH1F *h = new TH1F("h","demo quantiles",100,-3,3);
   h->FillRandom("gaus",5000);
 
   Double_t xq[nq];  // position where to compute the quantiles in [0,1]
   Double_t yq[nq];  // array to contain the quantiles
   for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
   h->GetQuantiles(nq,yq,xq);
 
   //show the original histogram in the top pad
   TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,700,900);
   c1->Divide(1,2);
   c1->cd(1);
   h->Draw();
 
   // show the quantiles in the bottom pad
   c1->cd(2);
   gPad->SetGrid();
   TGraph *gr = new TGraph(nq,xq,yq);
   gr->SetMarkerStyle(21);
   gr->Draw("alp");
}

Double_t xtoaitoffx(Double_t x, Double_t y)
{
	Double_t conv = 3.14159265353/180.;
	Double_t alpha = acos(cos(y * conv) * cos(x * conv / 2.0));
	return (2 * cos(y * conv) * sin(x * conv / 2.0))/(sin(alpha) / alpha) * 180/3.14159265353;
}

Double_t ytoaitoffy(Double_t x, Double_t y)
{
	Double_t conv = 3.14159265353/180.;
	Double_t alpha = acos(cos(y * conv) * cos(x * conv / 2.0));
	
	return sin(y * conv) / (sin(alpha) / alpha) * (1 / conv);
}

void cringe(Double_t x, Double_t y, string c)
{
	if(x > 180.0)
		x -= 360.;
	cout<<-1.0 * xtoaitoffx(x, y)<<" "<<ytoaitoffy(x, y)<<" "<<c<<endl;
}

void crinj(Double_t x, Double_t y, string c)
{
	if(x > 180.0)
		x -= 360.;
	cout<<"(int)("<<-1.0 * x<<" + 181), (int)("<<y<<" + 91)"<<endl;
}

void b()
{
	//INPUT EQ, OUTPUT AITOFF ADJUSTED EQ
	cringe(266.40508920      ,-28.93617470,"center");
	cringe(77.35810683       ,5.69314556,"txs");
	cringe(253.46769620     , 39.76020664,"501");
	cringe(166.11376144     , 38.20889554,"421");
	cringe(197.20311640     , -6.77748316,"v");
	cringe(299.86827636    ,  40.73390122,"cyg");
	cringe(201.36510056   ,   -43.01916045,"cen a");
	cringe(99.73521351   ,    -25.76973734,"dip");
	cringe(50.17781675,-33.73004768  ,"for");
	cringe(133.50301004 ,     43.11656811 ,"ta");
	cout<<endl;
	//INPUT SUPERGAL, OUTPUT ADJUSTED SUPERGAL
	cringe(185.78610785      ,42.31028736,"center");
	cringe(333.73561137      ,-56.21964474,"txs");
	cringe(68.09619712       ,54.30944729 ,"501");
	cringe(71.53130632       ,-10.57060882 ,"421");
	cringe(123.90968227      ,1.46441947,"v");
	cringe(0.44339945        ,61.33797902,"cyg");
	cringe(159.75324669      ,-5.24987415 ,"cen a");
	cringe(226.25147583      ,-79.26256273,"dip");
	cringe(265.52894514      ,-38.73367873,"for");
	cringe(50.03814291       ,-25.15294238,"ta");
	cout<<endl;
	//INPUT GAL, OUTPUT AITOFF ADJUSTED GAL
	cringe(0.,0.,"gal center");
	cringe(195.4053908316515-360.,-19.6360332367710,"txs");
	cringe(63.60003591638,38.85915603678,"mrk 501");
	cringe(179.8316666440783,65.0314774861495,"mrk 421");
	cringe(310.6281841000837-360, 55.8344218001367,"virgo");
	cringe(76.18987441564,05.75538794195,"cyg a");
	cringe(309.51589573409-360,19.41727341133,"cen a");
	cringe(235-360,-14,"dipole");
	cringe(233.8391652073152-360,-57.3349907989605,"Fornax");
	cringe(178,40,"TA hotspot");
	cout<<endl;
	
	crinj(266.40508920      ,-28.93617470,"center");
	crinj(77.35810683       ,5.69314556,"txs");
	crinj(253.46769620     , 39.76020664,"501");
	crinj(166.11376144     , 38.20889554,"421");
	crinj(197.20311640     , -6.77748316,"v");
	crinj(299.86827636    ,  40.73390122,"cyg");
	crinj(201.36510056   ,   -43.01916045,"cen a");
	crinj(99.73521351   ,    -25.76973734,"dip");
	crinj(50.17781675,-33.73004768  ,"for");
	crinj(133.50301004 ,     43.11656811 ,"ta");
	cout<<endl;
	//INPUT SUPERGAL, OUTPUT ADJUSTED SUPERGAL
	crinj(185.78610785      ,42.31028736,"center");
	crinj(333.73561137      ,-56.21964474,"txs");
	crinj(68.09619712       ,54.30944729 ,"501");
	crinj(71.53130632       ,-10.57060882 ,"421");
	crinj(123.90968227      ,1.46441947,"v");
	crinj(0.44339945        ,61.33797902,"cyg");
	crinj(159.75324669      ,-5.24987415 ,"cen a");
	crinj(226.25147583      ,-79.26256273,"dip");
	crinj(265.52894514      ,-38.73367873,"for");
	crinj(50.03814291       ,-25.15294238,"ta");
	cout<<endl;
	//INPUT GAL, OUTPUT AITOFF ADJUSTED GAL
	crinj(0.,0.,"gal center");
	crinj(195.4053908316515-360.,-19.6360332367710,"txs");
	crinj(63.60003591638,38.85915603678,"mrk 501");
	crinj(179.8316666440783,65.0314774861495,"mrk 421");
	crinj(310.6281841000837-360, 55.8344218001367,"virgo");
	crinj(76.18987441564,05.75538794195,"cyg a");
	crinj(309.51589573409-360,19.41727341133,"cen a");
	crinj(235-360,-14,"dipole");
	crinj(233.8391652073152-360,-57.3349907989605,"Fornax");
	crinj(178,40,"TA hotspot");
	
	
	cringe(40.66984680, -0.01329128, "NGC EQ");
	cringe(172.104011, -51.933614, "NGC GAL");
	cringe(304.26270862, -25.83846789, "NGC SUPERGAL");
	crinj(40.66984680, -0.01329128, "NGC");
	
}


void epic(){
   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",1300, 900);
   Double_t x[100], y[100];
   Int_t n = 100;
   for (Int_t i=0;i<n;i++) {
     x[i] = i*0.1;
     y[i] = 10*sin(x[i]+0.2);
   }
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetTitle("Option ACP example");
   gr->GetXaxis()->SetTitle("X title");
   gr->GetYaxis()->SetTitle("Y title");
   gPad->SetLogy(0);
   gPad->SetLogx(0);
   gr->Draw("AC");
}


void epic2()
{
	Double_t epic[21];
	epic[0] = 71717.2;
	epic[1] = 239867;
	epic[2] = 648930;
	epic[3] = 1.6161e+06;
	epic[4] = 3.84329e+06;
	epic[5] = 8.93447e+06;
	epic[6] = 1.87321e+07;
	epic[7] = 3.80246e+07;
	epic[8] = 7.23577e+07;
	epic[9] = 1.25462e+08;
	epic[10] = 1.97201e+08;
	epic[11] = 2.82618e+08;
	epic[12] = 3.54912e+08;
	epic[13] = 4.41724e+08;
	epic[14] = 5.20205e+08;
	epic[15] = 5.72186e+08;
	epic[16] = 6.26571e+08;
	epic[17] = 6.40303e+08;
	epic[18] = 6.58683e+08;
	epic[19] = 6.12582e+08;
	epic[20] = 6.26538e+08;
	
	Double_t integral = 0;
	for(int i = 0; i < 20; i++)
	{
		integral += (epic[i] + epic[i + 1]) / 2.0 * (pow(10, 6 + ((double)(i + 1)) * 0.2) - pow(10, 6 + ((double)(i)) * 0.2));
	}
	cout<<integral<<endl;
	cout<<integral * 10*365*24*3600*0.20<<endl;
	
	Double_t e[21];
	e[0] = 5.68602e-07 ;
	e[1] = 2.6944e-07 ;
	e[2] = 1.57846e-07 ;
	e[3] = 1.00453e-07 ;
	e[4] = 6.69466e-08 ;
	e[5] = 4.56418e-08 ;
	e[6] = 3.4502e-08 ;
	e[7] = 2.69382e-08 ;
	e[8] = 2.24361e-08 ;
	e[9] = 2.05079e-08 ;
	e[10] = 2.06787e-08 ;
	e[11] = 2.28682e-08 ;
	e[12] = 2.8861e-08 ;
	e[13] = 3.6752e-08 ;
	e[14] = 4.94604e-08 ;
	e[15] = 7.12681e-08 ;
	e[16] = 1.03148e-07 ;
	e[17] = 1.59973e-07 ;
	e[18] = 2.46465e-07 ;
	e[19] = 4.20018e-07 ;
	e[20] = 6.50856e-07 ;
	
	integral = 0;
		
	for(int i = 0; i < 20; i++)
	{
		integral += (e[i] + e[i + 1]) / 2.0 * (pow(10, 6 + ((double)(i + 1)) * 0.2) - pow(10, 6 + ((double)(i)) * 0.2));
	}
	//~ cout<<integral<<endl;
	//~ cout<<1.0 / integral<<endl; 
}

void wow()
{
	TFile *file2 =  TFile::Open("eq.root");
	TH2F *a = (TH2F*)file2->Get("skymapProjEq");
	a->Draw("z aitoff");
	
	//~ TFile *file22 =  TFile::Open("gal.root");
	//~ TH2F *a2 = (TH2F*)file22->Get("skymapFullProjection");
	//~ a2->Draw("z aitoff");
	
	const Int_t Number = 9;
	Double_t Red[Number]    = { 242./255., 234./255., 237./255., 230./255., 212./255., 156./255., 99./255., 45./255., 0./255.};
	Double_t Green[Number]  = { 243./255., 238./255., 238./255., 168./255., 101./255.,  45./255.,  0./255.,  0./255., 0./255.};
	Double_t Blue[Number]   = { 230./255.,  95./255.,  11./255.,   8./255.,   9./255.,   3./255.,  1./255.,  1./255., 0./255.};
	Double_t Length[Number] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
	Int_t nb=99;
	TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
	//~ TColor::InvertPalette();
	a->SetContour(nb);
}

void dekoms()
{
	TCanvas *c1 = new TCanvas("c1","Examples of Gaxis",10,10,900,500);
   c1->Range(-6,-0.1,6,0.1);
   TGaxis *axis1 = new TGaxis(-5.5,0.,5.5,0.,0.0,100,510,"");
   axis1->SetName("axis1");
   axis1->SetTitle("Axis Title");
   axis1->SetTitleSize(0.05);
   axis1->SetTitleColor(kBlue);
   axis1->SetTitleFont(42);
   axis1->ChangeLabel(1,-1,-1,-1,2);
   axis1->ChangeLabel(3,-1,0.);
   axis1->ChangeLabel(5,30.,-1,0);
   axis1->ChangeLabel(6,-1,-1,-1,3,-1,"6th label");
   axis1->ChangeLabel(-2,-1,-1,-1,-1,-1,"0");
   axis1->Draw();
}

void Piss()
{
	Double_t latitude = 38.52028; //lat of frisco peak, utah
	Double_t tStep = 0.25;
	Double_t LST = 0;
	
	TH2F *teleFOV = new TH2F("telefov", "Telescope FoV", 3601, -180.05, 180.05, 1801, -90.05, 90.05);
	
	int yBinMax = 900;
	int yBinMin = 662;
	
	for(int yBins = yBinMin; yBins <= yBinMax; yBins++)
	{
		for(int xBins = 1; xBins <= teleFOV->GetNbinsX(); xBins++)
		{
			teleFOV->SetBinContent(xBins, yBins, tStep * 4.0);
		}
	}
	
	teleFOV->Draw("z");
}

void hammer(Double_t lo, Double_t la)
{
	Double_t conv = 3.14159265353/180.;
	Double_t z  = sqrt(1+cos(la*conv)*cos(lo*conv/2));
    Double_t x  = 180*cos(la*conv)*sin(lo*conv/2)/z;
    Double_t yy  = 90*sin(la*conv)/z;  
      
    cout<<x * -1.<<", "<<yy<<endl;
}
