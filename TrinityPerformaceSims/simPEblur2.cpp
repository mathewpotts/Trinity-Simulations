Int_t iConfig = 2;
double DetectorAltitude[] = {0, 1, 2, 3};
double pi = 3.14159265359;
double lincorr[] = { 0.0509251, 0.0522854, 0.0595455, 0.0642221};
double scalefirst[] = { 1.00001, 1.03346, 1.53535, 2.36961};
double eleScaling[] = { 0.00163848, 0.00191408, 0.0071185, 0.0513182};
double absorptionlength[] = { 16.7049, 16.6106, 17.9808, 19.0274};
double parPEF[] = { 0.000332267, 0.580756, 4.25751e-05, 1.9491, 0.0427249};
double REarth = 6371; //km
Double_t c = 299792.458; //km/s
Double_t DecayTime = 0.290e-12; //s
Double_t Mtau = 1.7768; //GeV

Double_t myPEfunction2(Double_t azi, Double_t elv, Double_t l)
{
   if(azi>0.69813170)
     return 0;

   //Calculate azimuth angle in frame of master pe distribution (50km, 0ele,
   //0altitude
   Double_t dTelAngle = atan(DetectorAltitude[iConfig]*1e-3/l);
   Double_t dAngle = sqrt(azi*azi + (elv-dTelAngle)*(elv-dTelAngle))*57.295780; //in deg
   //calculate how many PEs / per m2 per GeV
   Double_t f = 0;
   if(dAngle<1.3)
     f = parPEF[0]*exp(-1.1/parPEF[1])+parPEF[2]*exp(-1.1/(parPEF[3]+1.1*parPEF[4]));
   else
     f = parPEF[0]*exp(-dAngle/parPEF[1])+parPEF[2]*exp(-dAngle/(parPEF[3]+dAngle*parPEF[4]));

//other parameters to get PE intensity for different distance, azimuth and
//elevation
   //scale PE distribution to first PE distribution at 50km distance
   f*=   scalefirst[iConfig];
   //Get elevation dependence
   f*= (2-exp(-elv/eleScaling[iConfig]));
   //Get Distance dependence
   f*=  exp(-(l-55)/
           (absorptionlength[iConfig]+(l-55)*lincorr[iConfig])); //55km is the distance for which the normalized PE distribution is extracted

   return f;
}

Double_t GH(Double_t *x, Double_t *p)
{
	Double_t xx = x[0];
	
	Double_t X_0 = p[0];
	Double_t X_max = p[1];
	Double_t lmbda = p[2];
	Double_t N_0 = p[3];
	
	return N_0 * pow(((xx - X_0) / (X_max - X_0)), ((X_max - X_0) / lmbda)) * exp((X_max - xx) / lmbda);
}

Double_t Gauss2D(Double_t *x, Double_t *p)
{
	Double_t xx = x[0];
	Double_t yy = x[1];
	
	Double_t sigma = p[0];
	
	return 1. / (2. * pi * sigma * sigma) * exp(-1. * (xx * xx + yy * yy) / (2. * sigma * sigma));
}

Double_t DistanceFovBelow(Double_t fovbelow, Double_t hdist)
{
	Double_t slope = tan(fovbelow * pi / 180.);
	Double_t x = (sqrt(slope * slope * REarth * REarth + 2. * hdist * slope * REarth - hdist * hdist * slope * slope) + slope * slope * hdist - slope * REarth) / (slope * slope + 1.);
	Double_t y = sqrt(REarth * REarth - x * x);
	
	return sqrt((hdist - x) * (hdist - x) + (REarth - y) * (REarth - y));
}

Double_t TauDistanceBehindHorizonToEarth(Double_t hdist, Double_t azi, Double_t elv, Double_t l)
{
	Double_t x = hdist - l;
	Double_t y = 0;
	Double_t z = REarth;
	
	Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t k = sin(elv * pi / 180.);
	
	Double_t a = 1; //unit vector
	Double_t b = -2. * (i * -1.* x + j * -1. * y + k * -1. * z);
	Double_t c = x * x + y * y + z * z - REarth * REarth;
	
	if(b * b - 4. * a * c <= 0) //tau does not intersect earth
		return -1;
	
	Double_t t = (-1. * b + sqrt(b * b - 4. * a * c)) / (2 * a);
	
	Double_t x1 = x + t * i; 
	Double_t y1 = y + t * j;
	Double_t z1 = z + t * k;
	
	return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1));
}

Double_t MaxdInCamBehindH(Double_t hdist, Double_t azi, Double_t elv, Double_t l, Double_t vfovaboveh, Double_t s)
{
	Double_t a = hdist; //cone
	Double_t b = 0;
	Double_t c = REarth;
	
	Double_t x = hdist - l; //AS
	Double_t y = 0;
	Double_t z = REarth;
	
	Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t k = sin(elv * pi / 180.);
	
	Double_t m = tan(vfovaboveh * pi / 180.);
	
	Double_t a1 = j * j * m * m - k * k + i * i * m * m;
	Double_t b1 = (2. * a * i * m * m + 2. * b * j * m * m - 2. * c * k - 2. * j * m * m * y + 2. * k * z - 2 * i * m * m * x) * -1.;
	Double_t c1 = a * a * m * m - 2. * a * m * m * x + b * b * m * m - 2. * b * m * m * y - c * c + 2. * c * z + m * m * x * x + m * m * y * y - z * z;
	
	
	if(b1 * b1 - 4. * a1 * c1 <= 0) //tau trajectory does not intersect with veritcal fov
		return -1;
	
	Double_t t = (-1. * b1 - sqrt(b1 * b1 - 4. * a1 * c1)) / (2. * a1);
	
	Double_t x1 = x + t * i; 
	Double_t y1 = y + t * j;
	Double_t z1 = z + t * k;
	
	Double_t d = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + (z - z1) * (z - z1)) - s;
	
	if(d <= 0) //length of air shower cannot fit within d
		return -1;
	
	return d;
}

Double_t TauTrajectoryHorizonDist(Double_t azi, Double_t elv, Double_t l, Double_t H)
{
	Double_t R = REarth;
	
	Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t k = sin(elv * pi / 180.);
	
	Double_t denom = H*H*j*j+H*H*k*k+k*k*R*R;
	Double_t radicand = -H*H*H*H*j*j*k*k*R*R+H*H*H*H*(-k*k*k*k)*R*R+2.*H*H*j*j*k*k*l*l*R*R+2.*H*H*k*k*k*k*l*l*R*R-j*j*k*k*l*l*l*l*R*R-k*k*k*k*l*l*l*l*R*R+4.*k*k*k*k*l*l*R*R*R*R;
	
	if(denom == 0 || radicand <= 0)
		return -1;

	return (-1.*(H*H*j*j*i*l*l)/(denom)-(H*H*k*k*i*l*l)/(denom)+(2.*H*H*k*k*i*R*R)/(denom)-(H*j*j*k*l*l*R)/(denom)-
					(H*k*k*k*l*l*R)/(denom)+(2.*H*k*k*k*R*R*R)/(denom)-H*H*i+(H*H*H*H*j*j*i)/(denom)+(H*H*H*H*k*k*i)/(denom)-(H*i*sqrt(radicand))/(denom)-
					(k*R*sqrt(radicand))/(denom)+(H*H*H*j*j*k*R)/(denom)+(H*H*H*k*k*k*R)/(denom)+i*l*l)/(2.*k*R);
}

Double_t * TauInFrontOfH(Double_t azi, Double_t elv, Double_t l, Double_t H, Double_t dist, Double_t vfovaboveh, Double_t s, Double_t hfov)
{
	static Double_t ret[8] = {0, 0, 0, 0, 0, 0, -1, 0};
	
	Double_t x_as = dist;
	Double_t y_as = 0;
	Double_t z_as = REarth;
	
	Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
	Double_t k = sin(elv * pi / 180.);
	
	Double_t a = 1; //unit vector
	Double_t b = -2. * (i * -1.* x_as + j * -1. * y_as + k * -1. * z_as);
	Double_t c = x_as * x_as + y_as * y_as + z_as * z_as - REarth * REarth;
	
	if(b * b - 4. * a * c <= 0) //tau does not intersect earth
	{
		cout << "Tau trajectory does not intersect Earth." << endl;
		return ret;
	}
	
	Double_t t = (-1. * b + sqrt(b * b - 4. * a * c)) / (2 * a);
	
	Double_t x_earth = x_as + t * i; 
	Double_t y_earth = y_as + t * j;
	Double_t z_earth = z_as + t * k;
	
	Double_t x_det = H; //detector pos
	Double_t y_det = 0;
	Double_t z_det = REarth;
	
	if(dist <= H) 
	{	
		Double_t m = tan(vfovaboveh * pi / 180.);
		
		Double_t a1 = j * j * m * m - k * k + i * i * m * m;
		Double_t b1 = (2. * x_det * i * m * m + 2. * y_det * j * m * m - 2. * z_det * k - 2. * j * m * m * y_as + 2. * k * z_as - 2 * i * m * m * x_as) * -1.;
		Double_t c1 = x_det * x_det * m * m - 2. * x_det * m * m * x_as + y_det * y_det * m * m - 2. * y_det * m * m * y_as - z_det * z_det + 2. * z_det * z_as + m * m * x_as * x_as + m * m * y_as * y_as - z_as * z_as;
		
		if(b1 * b1 - 4. * a1 * c1 <= 0)//tau trajectory does not intersect with veritcal fov
		{
			cout << "Tau trajectory does not intersect vertical FoV." << endl;
			return ret;
		}
	
		t = (-1. * b1 - sqrt(b1 * b1 - 4. * a1 * c1)) / (2. * a1);
		
		Double_t x_vfov = x_as + t * i; 
		Double_t y_vfov = y_as + t * j;
		Double_t z_vfov = z_as + t * k;
		
		Double_t d = sqrt((x_vfov - x_earth) * (x_vfov - x_earth) + (y_vfov - y_earth) * (y_vfov - y_earth) + (z_vfov - z_earth) * (z_vfov - z_earth)) - s;
		
		if(d <= 0)
		{
			cout << "Air shower cannot be contained within camera," << endl;
			return ret;
		}
		
		ret[0] = x_vfov - i * s;
		ret[1] = y_vfov - j * s;
		ret[2] = z_vfov - k * s;
		ret[3] = x_vfov;
		ret[4] = y_vfov;
		ret[5] = z_vfov;
		ret[6] = 1;
		ret[7] = d;
		
		cout << "x_earth: " << x_earth << endl;
		cout << "y_earth: " << y_earth << endl;
		cout << "z_earth: " << z_earth - REarth << endl;
		cout << "tau earth to det: " << l << endl;
		
		return ret;
	}
	
	Double_t m = -1. * tan(hfov * 0.5 * pi / 180.);
	
	Double_t denom = -m * i + j;
	
	if(denom == 0)
	{
		cout << "Air shower does not intersect horizontal FoV" << endl;
		return ret;
	}
	
	t = -(-m * x_as + y_as + m * x_det - y_det) / denom;
	
	Double_t x_hfov = x_as + t * i;
	Double_t y_hfov = y_as + t * j;
	Double_t z_hfov = z_as + t * k;
	
	Double_t d = sqrt((x_hfov - x_earth) * (x_hfov - x_earth) + (y_hfov - y_earth) * (y_hfov - y_earth) + (z_hfov - z_earth) * (z_hfov - z_earth)) - s;
	
	if(d <= 0)
	{
		cout << "Air shower cannot be contained within camera." << endl;
		return ret;
	}
	
	ret[0] = x_hfov - i * s;
	ret[1] = y_hfov - j * s;
	ret[2] = z_hfov - k * s;
	ret[3] = x_hfov;
	ret[4] = y_hfov;
	ret[5] = z_hfov;
	ret[6] = 2;
	ret[7] = d;
	
	return ret;
}

void SimPEinCam(Double_t azi, Double_t elv, Double_t l, Double_t E)
{
	Double_t E_c = 0.088; //crit energy, GeV
	Double_t X_0 = 33.662; //rad length, g/cm^2
	Double_t lmbda = 70; //decay length, g/cm^2
	Double_t X_n = log(E / E_c) / log(2); //# of rad lengths
	Double_t X_max = X_n * X_0; //Total rad length, g/cm^2
	Double_t air_p = 1.1073e-3; //air density, g/cm^3
	Double_t s = X_max / air_p * 1e-5; //shower length, km
	Double_t N_0 = 1; //max particles at X_max
		
	Double_t eff_mirror_size = 10; //m^2
	iConfig = 2; //km detector altitude
	
	Double_t totalPEs = myPEfunction2((azi * pi / 180.), (elv * pi / 180.), l) * eff_mirror_size * E;
	
	TString ghtitle;
	ghtitle.Form("PE Distribution for %1.0e GeV particle", E);
	
	TCanvas *ghc = new TCanvas(ghtitle,ghtitle,1000,800);
	ghc->cd(1);
	ghc->SetLeftMargin(0.12);
	
	TF1 *gh = new TF1("gh", GH, X_0, X_max, 4);
	gh->SetParameters(X_0, X_max, lmbda, N_0);
	Double_t norm = 1.0 / gh->Integral(X_0, X_max);
	N_0 = norm * totalPEs;
	
	gh->SetParameters(X_0, X_max, lmbda, N_0);
	gh->SetTitle(ghtitle);
	gh->GetXaxis()->SetTitle("Interaction depth (g/cm^{2})");
	gh->GetYaxis()->SetTitle("Number of PEs");
	ghc->SetLogy(0);
	ghc->SetLogx(0);
	gh->Draw();
	
	cout << "Total PEs produced: " << totalPEs << endl;
	cout << "Max number of PEs produced: " << N_0 << endl;
	cout << "Initial particle energy: " << E << " GeV" << endl;
	cout << "Critical energy: " << E_c << " Gev" << endl;
	cout << "Radiation length: " << X_0 << " g/cm^2" << endl;
	cout << "Decay length: " << lmbda << " g/cm^2" << endl;
	cout << "Total number of radiation lengths: " << X_n << endl;
	cout << "Total radiation length: " << X_max << " g/cm^2" << endl;
	cout << "Shower length: " << s << " km" << endl;
	cout << "Air density: " << air_p << " g/cm^3" << endl;
	
	cout << "\nShower azimuth: " << azi << " deg" << endl;
	cout << "Shower elevation: " << elv << " deg" << endl;
	cout << "Distance from particle emergence to detector: " << l << " km" << endl;
	
	Double_t hfov = 60; //deg
	Double_t vfov = 5; //deg
	Double_t pxwidth = 0.3; //deg
	Double_t pxheight = pxwidth; //deg
	Double_t vfovbelow = 3.0; //deg
	Double_t vfovabove = 2.0; //deg
	
	Double_t nxpx = ceil(hfov / pxwidth);
	Double_t nypx = ceil(vfov / pxheight);
	
	Double_t proj_alt, proj_hwidth, proj_length, ang_proj_alt, ang_proj_hwidth, ang_proj_length, length_to_tip;
	Double_t x1, y1, z1, x2, y2, z2;
	
	Double_t horizondist = sqrt((DetectorAltitude[iConfig] + REarth) * (DetectorAltitude[iConfig] + REarth) - REarth * REarth); //km
	Double_t vfovbelowdist = DistanceFovBelow(vfovbelow, horizondist); //km
	Double_t DecayLength = E * c * DecayTime / Mtau; //km
	
	bool infrontofH;
	
	cout << "Horizon distance: " << horizondist << " km" << endl;
	cout << "Max vertical FoV below distance : " << vfovbelowdist << " km" << endl;
	cout << "Tau decay length: " << DecayLength << " km" << endl;
	
	if(l >= horizondist) //tau emerges behind horizon
	{
		Double_t dprime = TauDistanceBehindHorizonToEarth(horizondist, azi, elv, l); //km, the distance the tau travels before entering FoV (when it crosses the horizon)
		
		if(dprime < 0) //tau does not intesect earth, abort and exit
		{
			cout << "Tau does not intersect Earth, no air shower can be created!" << endl;
			return;
		}
		
		if(exp(-1. * dprime / DecayLength) < 0.1) //90% of taus will decay before entering fov
		{
			cout << "90% of taus will decay before emerging into field of view, no air shower created!" << endl;
			return;
		}
		
		Double_t d = MaxdInCamBehindH(horizondist, azi, elv, l, vfovabove, s); //km, the distance tau travels from point of first emergence in the fov to the point of air shower initiation
		
		if(d < 0) //air shower to too large to be contained fully within the camera
		{
			cout << "Air shower cannot be contained within camera." << endl;
			return;
		}
		
		cout << "Percentage of taus that will decay and produce a fully-contained air shower: " << (1 - exp(-1. * (dprime + d) / DecayLength)) * 100. << "%" << endl;
		
		x1 = d * cos(azi * pi / 180.) * cos(elv * pi / 180.); //getting the coordinates of the shower base
		y1 = d * sin(azi * pi / 180.) * cos(elv * pi / 180.);
		z1 = d * sin(elv * pi / 180.);
		
		x2 = (d + s) * cos(azi * pi / 180.) * cos(elv * pi / 180.); //getting the coordinates of the shower tip
		y2 = (d + s) * sin(azi * pi / 180.) * cos(elv * pi / 180.);
		z2 = (d + s) * sin(elv * pi / 180.);
		
		length_to_tip = sqrt((-l + x2) * (-l + x2) + (y2) * (y2) + (z2) * (z2));
		
		proj_alt = s * sin(elv * pi / 180.);
		proj_hwidth = sin(azi  * pi / 180.) * s * cos(elv  * pi / 180.);
		proj_length = sqrt(proj_alt * proj_alt + proj_hwidth * proj_hwidth);
		
		ang_proj_alt = 180. / pi * (2 * atan2(proj_alt, (2 * length_to_tip)));
		ang_proj_hwidth = 180. / pi * (2 * atan2(proj_hwidth, (2 * length_to_tip)));
		ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
		
		infrontofH = false;
	}
	else if(l < horizondist && l >= vfovbelowdist) //tau emerges infront of horizon and in view of detector
	{
		Double_t x = TauTrajectoryHorizonDist(azi, elv, l, horizondist);
		Double_t *as_coords = TauInFrontOfH(azi, elv, l, horizondist, x, vfovabove, s, hfov);
		
		if(as_coords[6] < 0)
		{
			cout << "Air shower could not be constructed." << endl;
			return;
		}
		
		if(as_coords[6] == 1)
			cout << "Air shower reached top of FoV" << endl;
		else
			cout << "Air shower reached max horizontal FoV" << endl;
		
		cout << "Percentage of taus that will decay and produce a fully-contained air shower: " << (1 - exp(-1. * as_coords[7] / DecayLength)) * 100. << "%" << endl;
		
		x1 = as_coords[0]; //base
		y1 = as_coords[1];
		z1 = as_coords[2] - REarth;
		
		x2 = as_coords[3]; //tip
		y2 = as_coords[4];
		z2 = as_coords[5] - REarth;
		
		cout << "x1: " << x1 << endl;
		cout << "y1: " << y1 << endl;
		cout << "z1: " << z1 << endl;
		
		cout << "x2: " << x2 << endl;
		cout << "y2: " << y2 << endl;
		cout << "z2: " << z2 << endl;
		
		length_to_tip = sqrt((x2 - horizondist) * (x2 - horizondist) + (y2) * (y2) + (z2) * (z2));
		cout << "length to top: " << length_to_tip << endl;
		
		proj_alt = s * sin(elv * pi / 180.);
		proj_hwidth = sin(azi  * pi / 180.) * s * cos(elv  * pi / 180.);
		proj_length = sqrt(proj_alt * proj_alt + proj_hwidth * proj_hwidth);
		
		ang_proj_alt = 180. / pi * (2 * atan2(proj_alt, (2 * length_to_tip)));
		ang_proj_hwidth = 180. / pi * (2 * atan2(proj_hwidth, (2 * length_to_tip)));
		ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
		
		infrontofH = true;
	}
	else //tau emerges too close to the detector
	{
		infrontofH = true;
	}
	
	cout << "Altitude: " << proj_alt * 1000. << " m" << endl;
	cout << "Projected width: " << proj_hwidth * 1000. << " m" << endl;
	cout << "Projected length: " << proj_length * 1000. << " m" << endl;
	cout << "Angular altitude: " << ang_proj_alt << " deg" << endl;
	cout << "Projected angular width: " << ang_proj_hwidth << " deg" << endl;
	cout << "Projected angular length: " << ang_proj_length << " deg" << endl;
	
	Double_t distobase;
	Double_t showerbasex, showerbasey;
	
	if(!infrontofH)
	{
		distobase = sqrt((-l + x1) * (-l + x1) + (y1) * (y1) + (z1) * (z1));
		showerbasex = 180. / pi * (2 * atan2(y1, (2 * distobase))) + hfov * 0.5;
		showerbasey = 180. / pi * (2 * atan2(z1, (2 * distobase))) + vfovbelow;
	}
	else
	{
		distobase = sqrt((x1 - horizondist) * (x1 - horizondist) + (y1) * (y1) + (z1) * (z1));
		showerbasex = 180. / pi * (2 * atan2(y1, (2 * distobase)));
		showerbasey = 180. / pi * (2 * atan2(z1, (2 * distobase)));
	}
	
	cout << "showerbasex: " << showerbasex << endl;
	cout << "showerbasey: " << showerbasey << endl;
	cout << "showertipy: " << showerbasey + ang_proj_alt << endl;
	cout << "ang_proj_alt: " << ang_proj_alt << endl;
	cout << "dist to base: " << distobase << endl;
	cout << "alt alt: " << 180. / pi * (2 * atan2(z2, (2 * length_to_tip))) - 180. / pi * (2 * atan2(z1, (2 * distobase))) << endl;
	//~ return;
	
	Double_t showertipx = showerbasex + ang_proj_hwidth;
	Double_t showerslope = ang_proj_alt / ang_proj_hwidth;
	Double_t showerintercepty = showerbasey - showerslope * showerbasex;
	cout << "Shower formula: y = " << showerslope << "x + " << showerintercepty << endl;
	
	Double_t triggerthresh = 0.99999;
	Double_t oversample = 10;
	
	TString camtitle;
	camtitle.Form("Projected Air Shower Image In Camera (%.0f deg azimuth, %.0f deg elevation, %.0f km distance)", azi, elv, l);
	
	TCanvas *pec = new TCanvas(camtitle,camtitle,1700,800);
	pec->cd(1);
	pec->SetRightMargin(0.12);
	
	TH2F *camera = new TH2F("camera", camtitle, nxpx, 0, nxpx, nypx, 0, nypx);
	camera->GetXaxis()->SetTitle("Pixel x");
	camera->GetYaxis()->SetTitle("Pixel y");
	camera->GetZaxis()->SetTitle("PE count");
	
	TH2F *cameraOS = new TH2F("camOS", "camOS", nxpx * oversample, 0, nxpx, nypx * oversample, 0, nypx);
	
	Double_t xint_a, xint_b, x, y, totalcamPE = 0;
	Double_t xrange[2];
	Double_t int_range[2];
	
	for(int xbin = 1 ; xbin <= cameraOS->GetNbinsX(); xbin++) 
	{
		for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++) 
		{
			xint_a = xint_b = -1;
			
			//bottom
			y = (ybin - 1) * pxheight / oversample;
			x = (y - showerintercepty) / showerslope;
			if(x >= (xbin - 1) * pxwidth / oversample && x <= xbin * pxwidth / oversample && ((x >= showerbasex && x <= showertipx) || (x <= showerbasex && x >= showertipx)))
				xint_a = x;
			
			//top
			y = ybin * pxheight / oversample;
			x = (y - showerintercepty) / showerslope;
			if(x >= (xbin - 1) * pxwidth / oversample && x <= xbin * pxwidth / oversample && ((x >= showerbasex && x <= showertipx) || (x <= showerbasex && x >= showertipx)))
			{
				if(xint_a <= 0)
					xint_a = x;
				else
					xint_b = x;
			}
			
			//left
			x = (xbin - 1) * pxwidth / oversample;
			y = showerslope * x + showerintercepty;
			if(y >= (ybin - 1) * pxheight / oversample && y <= ybin * pxheight / oversample && ((x >= showerbasex && x <= showertipx) || (x <= showerbasex && x >= showertipx)))
			{
				if(xint_a <= 0)
					xint_a = x;
				else
					xint_b = x;
			}
					
			//right
			x = xbin * pxwidth / oversample;
			y = showerslope * x + showerintercepty;
			if(y >= (ybin - 1) * pxheight / oversample && y <= ybin * pxheight / oversample && ((x >= showerbasex && x <= showertipx) || (x <= showerbasex && x >= showertipx)))
			{
				if(xint_a <= 0)
					xint_a = x;
				else
					xint_b = x;
			}
			
			//~ cout << xint_a << ", " << xint_b << " | " << xbin << ", " << ybin << endl;
			
			//single crossing; only one variable fills
			if(xint_a > xint_b && xint_b < 0)
			{
				cout << "Single crossing at x = " << xint_a << " for bin (" << xbin << ", " << ybin << ")" << endl;
				if(abs(showerbasex - xint_a) < abs(showertipx - xint_a)) //base crossing
				{
					xrange[0] = showerbasex;
					xrange[1] = xint_a;
					cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
					
					int_range[0] = X_0;
					int_range[1] = sqrt((xrange[1] - showerbasex)*(xrange[1] - showerbasex) + 
						((xrange[1] * showerslope + showerintercepty) - showerbasey) * ((xrange[1] * showerslope + showerintercepty) - showerbasey)) / ang_proj_length * (X_max - X_0) + X_0;
					
					Double_t PEs = gh->Integral(int_range[0], int_range[1]);
					
					cameraOS->SetBinContent(xbin, ybin, PEs);
					camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), PEs);
				}
				else
				{
					xrange[0] = xint_a;
					xrange[1] = showertipx;
					cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
					
					int_range[1] = X_max;
					int_range[0] = sqrt((xrange[0] - showerbasex) * (xrange[0] - showerbasex) + 
						((xrange[0] * showerslope + showerintercepty) - showerbasey) * ((xrange[0] * showerslope + showerintercepty) - showerbasey)) / ang_proj_length * (X_max - X_0) + X_0;
					
					Double_t PEs = gh->Integral(int_range[0], int_range[1]);
					
					cameraOS->SetBinContent(xbin, ybin, PEs);
					camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), PEs);
				}
			}
			else if(xint_a > xint_b)
			{
				cout << "Double crossing at x1 = " << xint_b << " x2 = " << xint_a << " for bin (" << xbin << ", " << ybin << ")" << endl;
				xrange[0] = xint_b;
				xrange[1] = xint_a;
				cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
				
				int_range[0] = sqrt((xrange[0] - showerbasex) * (xrange[0] - showerbasex) + 
						((xrange[0] * showerslope + showerintercepty) - showerbasey) * ((xrange[0] * showerslope + showerintercepty) - showerbasey)) / ang_proj_length * (X_max - X_0) + X_0;
				int_range[1] = sqrt((xrange[1] - showerbasex)*(xrange[1] - showerbasex) + 
						((xrange[1] * showerslope + showerintercepty) - showerbasey) * ((xrange[1] * showerslope + showerintercepty) - showerbasey)) / ang_proj_length * (X_max - X_0) + X_0;
				
				Double_t PEs = gh->Integral(int_range[0], int_range[1]);
					
				cameraOS->SetBinContent(xbin, ybin, PEs);
				camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), PEs);
			}
			else if(xint_a < xint_b)
			{
				cout << "Double crossing at x1 = " << xint_a << " x2 = " << xint_b << " for bin (" << xbin << ", " << ybin << ")" << endl;
				xrange[0] = xint_a;
				xrange[1] = xint_b;
				cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
				
				int_range[0] = sqrt((xrange[0] - showerbasex) * (xrange[0] - showerbasex) + 
						((xrange[0] * showerslope + showerintercepty) - showerbasey) * ((xrange[0] * showerslope + showerintercepty) - showerbasey)) / ang_proj_length * (X_max - X_0) + X_0;
				int_range[1] = sqrt((xrange[1] - showerbasex)*(xrange[1] - showerbasex) + 
						((xrange[1] * showerslope + showerintercepty) - showerbasey) * ((xrange[1] * showerslope + showerintercepty) - showerbasey)) / ang_proj_length * (X_max - X_0) + X_0;
				
				Double_t PEs = gh->Integral(int_range[0], int_range[1]);
					
				cameraOS->SetBinContent(xbin, ybin, PEs);
				camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), PEs);
			}
		}
	}
	
	for(int i = 1 ; i <= camera->GetNbinsX(); i++)
	{
		for(int j = 1 ; j <= camera->GetNbinsY(); j++)
		{
			if(camera->GetBinContent(i, j) > triggerthresh)
				totalcamPE += camera->GetBinContent(i, j);
			else if(camera->GetBinContent(i, j) > 0)
				camera->SetBinContent(i, j, 0);
		}
	}
	
	cout.precision(17);
	cout << "Total PEs produced: " << (int) totalPEs << endl;
	cout << "Total captured PEs: " << (int) totalcamPE << endl;
	
	camera->SetContour(99);
	camera->SetStats(0);
	camera->Draw("COLZ");
	
	Double_t blursigma = 1 * oversample;
	Double_t kernrad = 6 * blursigma; //sigmas
	
	TString camtitleb;
	camtitleb.Form("Projected Air Shower Image In Camera With %1.1f deg Gaussian Blur (%.0f deg azimuth, %.0f deg elevation, %.0f km distance)", pxwidth * blursigma / oversample, azi, elv, l);
	
	TCanvas *pecBlur = new TCanvas(camtitleb,camtitleb,1200,1200);
	pecBlur->cd(1);
	pecBlur->SetRightMargin(0.12);
	
	TH2F *cameraBlur = new TH2F(camtitleb, camtitleb, nxpx, 0, nxpx, nypx, 0, nypx);
	cameraBlur->GetXaxis()->SetTitle("Pixel x");
	cameraBlur->GetYaxis()->SetTitle("Pixel y");
	cameraBlur->GetZaxis()->SetTitle("PE count");
	
	Double_t GMatrix[(int)((1 + 2 * kernrad))][(int)((1 + 2 * kernrad))];
	
	TF2 *G2D = new TF2("g2d", Gauss2D, -1. * kernrad, kernrad, -1. * kernrad, kernrad , 1);
	G2D->SetParameter(0, blursigma);
	Double_t gweight = 0;
	
	for(int i = 0; i < (int)(1 + 2 * kernrad); i++)
	{
		for(int j = 0; j < (int)(1 + 2 * kernrad); j++)
		{
			GMatrix[i][j] = G2D->Eval(j - kernrad, kernrad - i);
			gweight += G2D->Eval(j - kernrad, kernrad - i);
		}
	}
	
	for(int i = 0; i < (int)(1 + 2 * kernrad); i++)
	{
		for(int j = 0; j < (int)(1 + 2 * kernrad); j++)
		{
			GMatrix[i][j] = GMatrix[i][j] / gweight;
		}
	}
	
	for(int xbin = 1; xbin <= cameraOS->GetNbinsX(); xbin++)
	{
		for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++)
		{
			if(cameraOS->GetBinContent(xbin, ybin) > 0)
			{
				for(int x = xbin - kernrad; x <= xbin + kernrad; x++)
				{
					for(int y = ybin - kernrad; y <= ybin + kernrad; y++)
					{
						if((x >= 1 && x <= cameraOS->GetNbinsX()) && (y >= 1 && y <= cameraOS->GetNbinsY()))
						{	
							cameraBlur->Fill(cameraOS->GetXaxis()->GetBinCenter(x), cameraOS->GetYaxis()->GetBinCenter(y), 
							GMatrix[(int)(y - ybin + kernrad)][(int)(x - xbin + kernrad)] * cameraOS->GetBinContent(xbin, ybin));
						}
					}
				}
			}
		}
	}
	
	totalcamPE = 0;
	
	for(int i = 1 ; i <= cameraBlur->GetNbinsX(); i++)
	{
		for(int j = 1 ; j <= cameraBlur->GetNbinsY(); j++)
		{
			if(cameraBlur->GetBinContent(i, j) > triggerthresh)
				totalcamPE += cameraBlur->GetBinContent(i, j);
			else if(cameraBlur->GetBinContent(i, j) > 0)
				cameraBlur->SetBinContent(i, j, 0);
		}
	}
	
	cout << "Total captured PEs (blur): " << (int) totalcamPE << endl;
	
	int minbinx = cameraBlur->GetNbinsX() + 20, maxbinx = -10, minbiny = cameraBlur->GetNbinsY() + 20, maxbiny = -10;
	
	for(int i = 1; i <= cameraBlur->GetNbinsX(); i++)
	{
		for(int j = 1; j <= cameraBlur->GetNbinsY(); j++)
		{
			if(cameraBlur->GetBinContent(i, j))
			{
				if(i < minbinx)
					minbinx = i;
				if(i > maxbinx)
					maxbinx = i;
				if(j < minbiny)
					minbiny = j;
				if(j > maxbiny)
					maxbiny = j;
			}
		}
	}
	
	if(maxbinx - minbinx > maxbiny - minbiny)
	{
		int range = maxbinx - minbinx;
		
		if(range > cameraBlur->GetNbinsY())
		{
			TH2F *cameraBlurBigger = new TH2F("bigger", camtitleb, nxpx, 0, nxpx, nxpx, 0, nxpx);
			cameraBlurBigger->GetXaxis()->SetTitle("Pixel x");
			cameraBlurBigger->GetYaxis()->SetTitle("Pixel y");
			cameraBlurBigger->GetZaxis()->SetTitle("PE count");
			
			for(int i = 1; i < cameraBlur->GetNbinsX(); i++)
			{
				for(int j = 1; j < cameraBlur->GetNbinsY(); j++)
				{
					cameraBlurBigger->SetBinContent(i, j, cameraBlur->GetBinContent(i, j));
				}
			}
			
			cameraBlurBigger->GetXaxis()->SetRangeUser(minbinx, maxbinx);
			cameraBlurBigger->GetYaxis()->SetRangeUser((maxbiny + minbiny) / 2. - range / 2., (maxbiny + minbiny) / 2. + range / 2.);
			
			cameraBlurBigger->SetContour(99);
			cameraBlurBigger->SetStats(0);
			cameraBlurBigger->Draw("COLZ");
		}
		else
		{
			cameraBlur->GetXaxis()->SetRangeUser(minbinx, maxbinx);
			cameraBlur->GetYaxis()->SetRangeUser((maxbiny + minbiny) / 2. - range / 2., (maxbiny + minbiny) / 2. + range / 2.);
			cameraBlur->SetContour(99);
			cameraBlur->SetStats(0);
			cameraBlur->Draw("COLZ");
		}
	}
	else
	{
		int range = maxbiny - minbiny;
		
		cameraBlur->GetYaxis()->SetRangeUser(minbiny, maxbiny);
		cameraBlur->GetXaxis()->SetRangeUser((maxbinx + minbinx) / 2. - range / 2., (maxbinx + minbinx) / 2. + range / 2.);
		
		cameraBlur->SetContour(99);
		cameraBlur->SetStats(0);
		cameraBlur->Draw("COLZ");
	}
}
