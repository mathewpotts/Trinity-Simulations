#include <iostream>
#include <string>

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
#include <TTree.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <Getline.h>

using namespace std;

Int_t iConfig = 2;
double DetectorAltitude[] = {0, 1, 2, 3};
double pi = 3.14159265359;
double REarth = 6371; //km
Double_t c = 299792.458; //km/s
Double_t DecayTime = 0.290e-12; //s
Double_t Mtau = 1.7768; //GeV

// Obtained from 3e4 GeV gamma rays
double lincorr[] = { 0.0509251, 0.0522854, 0.0595455, 0.0642221};
double scalefirst[] = { 1.00001, 1.03346, 1.53535, 2.36961};
double eleScaling[] = { 0.00163848, 0.00191408, 0.0071185, 0.0513182};
double absorptionlength[] = { 16.7049, 16.6106, 17.9808, 19.0274};

// Parameterization of PE distribution at 50km,0ele, and 0 altitude
double parPEF[] = { 0.000332267, 0.580756, 4.25751e-05, 1.9491, 0.0427249};

/////////////////////////////
// Functions:

// Function is from Nepomuk's TrinityPerformaceCalculation Code
// Slightly modifed to have 3 inputs rather than 2
Double_t myPEfunction2(Double_t azi, Double_t elv, Double_t l)
{
  // azi is azimuth angle in rad
  // elv is elevation angle in rad
  // l is distance to where tau comes out in km
  
  if(azi>0.69813170) // if > 40 degrees then return 0
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

// Just a place holder now, modeled after Nepomuk's PE function. 
Double_t myPhotonArrivalFunction(Double_t azi, Double_t elv, Double_t l)
{
  // azi is azimuth angle in rad
  // elv is elevation angle in rad
  // l is distance to where tau comes out in km

  if (azi>0.69813170) // if > 40 degrees then return 0
    return 0;
  
  // Calculate azimuth angle in frame of master photon arrival distribution
  Double_t dTelAngle = atan(DetectorAltitude[iConfig]*1e-3/l);
  Double_t dAngle = sqrt(azi*azi + (elv-dTelAngle)*(elv-dTelAngle))*57.295780; //in deg
  
  // Calculate the timing spread of photons
  Double_t f = 0;
  if (dAngle < 1.3)
    f = ;
  else
    f = ;
  
  return f;
}

// Function is from Nepomuk's TrinityPerformaceCalculation Code
Double_t DistanceThroughEarth(Double_t azimuth, Double_t elevation, Double_t y)
{
  elevation = elevation * TMath::DegToRad(); //elevation angle in rad 
                                             //(determines path through Earth)
  azimuth = azimuth * TMath::DegToRad();  //azimuth angle in rad
  
  Double_t l = y; //Distance from detector to where the tau comes out 
                  //detector is always at z=0
  
  Double_t v = sqrt((REarth+DetectorAltitude[iConfig])*(REarth+DetectorAltitude[iConfig])-REarth*REarth);
  
  //shortest distance d between tau trajectory and detector
  Double_t nproj = y*sqrt( 1 + tan(azimuth)*tan(azimuth) ); //projection of trajectory to x-y plane
  Double_t denomsquared= y*tan(azimuth)*y*tan(azimuth) + y*y + nproj*nproj*tan(elevation)*tan(elevation);
  
  //normalized trajectory vector of tau
  Double_t dNormalize = y/sqrt(denomsquared);
  //Double_t dNx = dNormalize * tan(azimuth);
  Double_t dNy = -dNormalize;
  Double_t dNz = dNormalize * sqrt( 1 + tan(azimuth)*tan(azimuth) ) * tan(elevation);
  
  Double_t p = 2 * ( REarth*dNz - (v-l)*dNy );
  Double_t q = (v-l)*(v-l);
  
  if(q-p*p/4>=0) //trajectory does not intersect with Earth
    return 0;
  
  Double_t i1 = p/2. - sqrt(p*p/4.-q);
  Double_t i2 = p/2. + sqrt(p*p/4.-q);
  
  return fabs(i2-i1);
}

// Gaisser-Hillas Function
Double_t GH(Double_t *x, Double_t *p) 
{
  Double_t xx = x[0];
  
  Double_t X_0 = p[0];
  Double_t X_max = p[1];
  Double_t lmbda = p[2];
  Double_t N_0 = p[3];
  
  return N_0 * pow(((xx - X_0) / (X_max - X_0)), ((X_max - X_0) / lmbda)) * exp((X_max - xx) / lmbda);
}

//2D Gaussian Function
Double_t Gauss2D(Double_t *x, Double_t *p) 
{
  Double_t xx = x[0];
  Double_t yy = x[1];
  
  Double_t sigma = p[0];
  
  return 1. / (2. * pi * sigma * sigma) * exp(-1. * (xx * xx + yy * yy) / (2. * sigma * sigma));
}

// Returns the distance in km from the detector to the lowest point 
// of the detector's lower fov
Double_t DistanceFovBelow(Double_t fovbelow, Double_t hdist) 
{
  Double_t slope = tan(fovbelow * pi / 180.);
  Double_t x = (sqrt(slope * slope * REarth * REarth + 2. * hdist * slope * REarth - hdist * hdist * slope * slope) + slope * slope * hdist - slope * REarth) / (slope * slope + 1.);
  Double_t y = sqrt(REarth * REarth - x * x);
  
  return sqrt((hdist - x) * (hdist - x) + (REarth - y) * (REarth - y));
}

// Returns the coordinates of tau emergence wrt to the center of the earth
Double_t *GetTauEarthCoords(Double_t azi, Double_t elv, Double_t l, Double_t hdist) 
{
  static Double_t ret[4] = {0, 0, 0, -1}; //4th element is an error code
  
  Double_t x_as = hdist - l; //coords of tau trajectory intersecting detector-horizon plane
  Double_t y_as = 0;
  Double_t z_as = REarth;
  
  Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.); //unit vector
  Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
  Double_t k = sin(elv * pi / 180.);
  
  //calculating the intersection of a line and a sphere
  Double_t a = 1; //unit vector has length 1 
  Double_t b = -2. * (i * -1.* x_as + j * -1. * y_as + k * -1. * z_as);
  Double_t c = x_as * x_as + y_as * y_as + z_as * z_as - REarth * REarth;
  
  if(b * b - 4. * a * c <= 0) //tau does not intersect earth
    return ret;
  
  Double_t t = (-1. * b + sqrt(b * b - 4. * a * c)) / (2 * a);
  
  ret[0] = x_as + t * i;
  ret[1] = y_as + t * j;
  ret[2] = z_as + t * k;
  ret[3] = 1;
  
  return ret;
}

// Returns the coordinates of the intersecting point of the tau 
// trajectory and the upper or lower fov
Double_t *TauTrajectoryFoV(Double_t *coords, Double_t azi, Double_t elv, Double_t fov, Double_t hdist, bool lower) 
{
  static Double_t ret[4] = {0, 0, 0, -1}; //4th element is an error code
  
  Double_t x_e = coords[0]; //tau emergence coords
  Double_t y_e = coords[1];
  Double_t z_e = coords[2];
  
  Double_t x_det = hdist; //detector coords
  Double_t y_det = 0;
  Double_t z_det = REarth;
  
  Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.); //unit vector
  Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
  Double_t k = sin(elv * pi / 180.);
  
  Double_t m = tan(fov * pi / 180.); //'slope' of the cone
  
  //calculating the intersection between a line and a cone
  Double_t a = j * j * m * m - k * k + i * i * m * m;
  Double_t b = (2. * x_det * i * m * m + 2. * y_det * j * m * m - 2. * z_det * k - 2. * j * m * m * y_e + 2. * k * z_e - 2 * i * m * m * x_e) * -1.;
  Double_t c = x_det * x_det * m * m - 2. * x_det * m * m * x_e + y_det * y_det * m * m - 2. * y_det * m * m * y_e - z_det * z_det + 2. * z_det * z_e + m * m * x_e * x_e + m * m * y_e * y_e - z_e * z_e;
  
  Double_t t1 = (-1. * b + sqrt(b * b - 4. * a * c)) / (2. * a);
  Double_t t2 = (-1. * b - sqrt(b * b - 4. * a * c)) / (2. * a);
  
  Double_t x1 = x_e + t1 * i;
  Double_t y1 = y_e + t1 * j;
  Double_t z1 = z_e + t1 * k;
  
  Double_t x2 = x_e + t2 * i;
  Double_t y2 = y_e + t2 * j;
  Double_t z2 = z_e + t2 * k;
  
  //~ printf("1) %f, %f, %f\n", x1, y1, z1 - z_det);
  //~ printf("2) %f, %f, %f\n", x2, y2, z2 - z_det);
  
  if(lower) //if we are dealing with the lower fov cone
    {
      if(x1 >= x_e && x1 < x_det && z1 < z_det) //check if the coordinate is in front of where the tau emerged, in front and below the detector
	{
	  ret[0] = x1;
	  ret[1] = y1;
	  ret[2] = z1;
	  ret[3] = 1;
	  
	  return ret;
	}
      if(x2 >= x_e && x2 < x_det && z2 < z_det)
	{
	  ret[0] = x2;
	  ret[1] = y2;
	  ret[2] = z2;
	  ret[3] = 1;
	  printf("second one lower\n");
	  return ret;
	}
    }
  else //we deal withthe upper fov cone
    {
      if(x1 >= x_e && x1 < x_det && z1 > z_det) //check if the coordinate is in front of where the tau emerged, in front and above the detector
	{
	  ret[0] = x1;
	  ret[1] = y1;
	  ret[2] = z1;
	  ret[3] = 1;
	  
	  return ret;
	}
      if(x2 >= x_e && x2 < x_det && z2 > z_det)
	{
	  ret[0] = x2;
	  ret[1] = y2;
	  ret[2] = z2;
	  ret[3] = 1;
	  //~ printf("second one upper\n");
	  return ret;
	}
    }
  
  return ret;
}

// Check if the tau emerges below the lower fov of the detector
bool TauEmergeBelowVFoV(Double_t *coords, Double_t vfovbelow, Double_t hdist) 
{
  Double_t x_e = coords[0]; //tau emergence coords
  Double_t y_e = coords[1];
  Double_t z_e = coords[2];
  
  Double_t x_det = hdist; //detector coords
  Double_t y_det = 0;
  Double_t z_det = REarth;
  
  Double_t m = tan(vfovbelow * pi / 180.); //'slope' of the cone
  
  //evaluating the cone formula
  Double_t LHS = m * m * (x_e - x_det) * (x_e - x_det) + m * m * (y_e - y_det) * (y_e - y_det);
  Double_t RHS = (z_e - z_det) * (z_e - z_det);
  
  if(LHS > RHS) //if the emergence point falls outside lower fov cone
    return false;
  
  return true;
}

// Checks to see if the air shower is observable
bool AirShowerAppearsInVFoV(Double_t *coords, Double_t dprime, Double_t d, Double_t azi, Double_t elv, Double_t l, Double_t vfovabove, Double_t vfovbelow, Double_t hdist, Double_t s, short tau_config)
{
  Double_t x_e = coords[0]; //tau emergence coords
  Double_t y_e = coords[1];
  Double_t z_e = coords[2];
  
  Double_t x_det = hdist; //detector coords
  Double_t y_det = 0;
  Double_t z_det = REarth;
  
  Double_t i = cos(azi * pi / 180.) * cos(elv * pi / 180.); //unit vector
  Double_t j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
  Double_t k = sin(elv * pi / 180.);
  
  Double_t x_base = (d + dprime) * i + x_e;  //air shower base coords
  Double_t y_base = (d + dprime) * j + y_e;
  Double_t z_base = (d + dprime) * k + z_e;
  
  Double_t m_u = tan(vfovabove * pi / 180.); //'slope' of the upper fov cone
  Double_t m_l = tan(vfovbelow * pi / 180.); //'slope' of the lower fov cone
  
  //evaluating the upper and lower cone equations
  Double_t LHS_upper = m_u * m_u * (x_base - x_det) * (x_base - x_det) + m_u * m_u * (y_base - y_det) * (y_base - y_det);
  Double_t RHS_upper = (z_base - z_det) * (z_base - z_det);
  Double_t LHS_lower = m_l * m_l * (x_base - x_det) * (x_base - x_det) + m_l * m_l * (y_base - y_det) * (y_base - y_det);
  Double_t RHS_lower = (z_base - z_det) * (z_base - z_det);
  
  //~ printf("(%f, %f, %f)\n", x_base, y_base, z_base - REarth);
  //~ cout << sqrt((x_det - x_base) * (x_det - x_base) + (y_det - y_base) * (y_det - y_base) + (z_det - z_base) * (z_det - z_base)) << endl;
  
  if(LHS_upper >= RHS_upper && z_base >= z_det && x_base <= x_det) //base is within upper fov, above horizon, and in front of detector
    return true;
  if(LHS_lower >= RHS_lower && z_base < z_det && x_base <= x_det) //base is within lower fov, below horizon, and in front of detector
    {
      if(tau_config == 0) //base is beyond the horizon
	{
	  Double_t x_plane = hdist - l; //setting up the coordinates of the tau trajectory and detector-horizon plane intersection
	  Double_t y_plane = 0;
	  Double_t z_plane = REarth;
	  
	  //calculating the distance between the base of the shower and the plane intersection point
	  Double_t dist = sqrt((x_plane - x_base) * (x_plane - x_base) + (y_plane - y_base) * (y_plane - y_base) + (z_plane - z_base) * (z_plane - z_base));
	  if(s > dist) //if the length of the shower is larger than the distance between the base of the shower and the horizon-detector plane
	    return true;
	  return false;
	}
      return true;
    }
  if(LHS_lower < RHS_lower && x_base < x_det && tau_config == 2) //tau emerges close to and below the lower fov of detector
    {
      Double_t *as_c = TauTrajectoryFoV(coords, azi, elv, vfovbelow, hdist, true); //get the coordinates of the lower fov tau intersection
      
      //calculat the distance between the base of the shower and the lower fov intersection
      Double_t dist = sqrt((as_c[0] - x_base) * (as_c[0] - x_base) + (as_c[1] - y_base) * (as_c[1] - y_base) + (as_c[2] - z_base) * (as_c[2] - z_base));
      
      if(s > dist) //if the distance from the base to the intersection of the lower fov is less that the air shower length
	{
	  //~ printf("dist: %f, s: %f\n", dist, s);
	  return true;
	}
      return false;
    }
  if(LHS_lower < RHS_lower && x_base < x_det && tau_config == 0) //extreme edge case for beyond the horizon emergence where base of shower is not in lower fov
    {
      Double_t x_plane = hdist - l; //setting up the coordinates of the tau trajectory and detector-horizon plane intersection
      Double_t y_plane = 0;
      Double_t z_plane = REarth;
      
      //calculating the distance between the base of the shower and the plane intersection point
      Double_t dist = sqrt((x_plane - x_base) * (x_plane - x_base) + (y_plane - y_base) * (y_plane - y_base) + (z_plane - z_base) * (z_plane - z_base));
      if(s > dist) //if the length of the shower is larger than the distance between the base of the shower and the horizon-detector plane
	return true;
      return false;
    }
  
  return false;
}

bool *ASCutoff(Double_t *base, Double_t *tip, Double_t vfovbelow, Double_t vfovabove, Double_t hdist)
{
  static bool ret[2] = {true, true}; //base, tip
  
  Double_t x_det = hdist; //detector coords
  Double_t y_det = 0;
  Double_t z_det = REarth;
  
  Double_t m_u = tan(vfovabove * pi / 180.); //'slope' of the upper fov cone
  Double_t m_l = tan(vfovbelow * pi / 180.); //'slope' of the lower fov cone
  
  //evaluating the upper and lower fov cone equations for the base
  Double_t LHS_upper = m_u * m_u * (base[0] - x_det) * (base[0] - x_det) + m_u * m_u * (base[1] - y_det) * (base[1] - y_det); //base first
  Double_t RHS_upper = (base[2] - z_det) * (base[2] - z_det);
  Double_t LHS_lower = m_l * m_l * (base[0] - x_det) * (base[0] - x_det) + m_l * m_l * (base[1] - y_det) * (base[1] - y_det);
  Double_t RHS_lower = (base[2] - z_det) * (base[2] - z_det);
  
  if(LHS_upper >= RHS_upper && base[2] >= z_det && base[0] <= x_det) //base is within upper fov
    ret[0] = false;
  else if(LHS_lower >= RHS_lower && base[0] >= 0 && base[2] < z_det) //base is within lower fov
    ret[0] = false;
  
  //evaluating the upper and lower fov cone equations for the tip
  LHS_upper = m_u * m_u * (tip[0] - x_det) * (tip[0] - x_det) + m_u * m_u * (tip[1] - y_det) * (tip[1] - y_det); //now tip
  RHS_upper = (tip[2] - z_det) * (tip[2] - z_det);
  LHS_lower = m_l * m_l * (tip[0] - x_det) * (tip[0] - x_det) + m_l * m_l * (tip[1] - y_det) * (tip[1] - y_det);
  RHS_lower = (tip[2] - z_det) * (tip[2] - z_det);
  
  if(LHS_upper >= RHS_upper && tip[2] >= z_det && tip[0] <= x_det) //tip is within upper fov
    ret[1] = false;
  else if(LHS_lower >= RHS_lower && tip[0] >= 0 && tip[2] < z_det) //base is within lower fov
    ret[1] = false;
  
  return ret;
}

// We expect to have showers with energies greater than 1 PeV
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
  
  Double_t totalPEs = myPEfunction2(azi * TMath::DegToRad(), elv * TMath::DegToRad(), l) * eff_mirror_size * E; //calculating the average total PEs expected from the shower
  
  // Checks to exit
  if (totalPEs == 0) //exit if shower azimuth angle >40 degrees
    {
      cout << "Azimuth angle too large." << endl;
      return;
    }
  
  TString ghtitle;
  ghtitle.Form("PE Distribution for %.0e GeV particle", E);
  
  TCanvas *ghc = new TCanvas(ghtitle,ghtitle,1000,800);
  ghc->cd(1);
  ghc->SetLeftMargin(0.12);
  
  TF1 *gh = new TF1("gh", GH, X_0, X_max, 4); //create the gaisser hillas function and initiate parameters and domain
  gh->SetParameters(X_0, X_max, lmbda, N_0);
  Double_t norm = 1.0 / gh->Integral(X_0, X_max); //normalizing the gaisser hillas function s.t. integrating from X_0 -> X_max yields totalPEs
  N_0 = norm * totalPEs;
  
  gh->SetParameters(X_0, X_max, lmbda, N_0);
  gh->SetTitle(ghtitle);
  gh->GetXaxis()->SetTitle("Interaction depth (g/cm^{2})");
  gh->GetYaxis()->SetTitle("Number of PEs");
  ghc->SetLogy(0);
  ghc->SetLogx(0);
  gh->Draw(); //drawing a graph of the gaisser hillas function
  
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
  cout << "Distance from detector to tau trajectory intersection with detector-horizon plane: " << l << " km" << endl;
  
  Double_t hfov = 60; //detector horizontal fov, deg
  Double_t vfov = 5; //detector vertical fov, deg
  Double_t pxwidth = 0.3; //horizontal angular resolution of each pixel, deg
  Double_t pxheight = pxwidth; //ditto, vertical, deg
  Double_t vfovbelow = 3.0; //vertical fov below the horizon, deg
  Double_t vfovabove = vfov - vfovbelow; //vertical fov above the horizon, deg
  
  Double_t nxpx = ceil(hfov / pxwidth); //number of pixels, x
  Double_t nypx = ceil(vfov / pxheight); //number of pixels, y
  
  Double_t horizondist = sqrt((DetectorAltitude[iConfig] + REarth) * (DetectorAltitude[iConfig] + REarth) - REarth * REarth); //distance from detector to horizon, km
  Double_t vfovbelowdist = DistanceFovBelow(vfovbelow, horizondist); //distance from detector to the lowest point of the lower fov, km
  Double_t DecayLength = E * c * DecayTime / Mtau; //tau decay length, km
  
  short tau_config = 0; //0 = beyond the horizon, 1 = in front of horizon but beyond detector fov, 2 = in front of horizon and in front of detector fov
  
  cout << "Horizon distance: " << horizondist << " km" << endl;
  cout << "Max vertical FoV below distance : " << vfovbelowdist << " km" << endl;
  cout << "Tau decay length: " << DecayLength << " km" << endl;
  
  /* For all coordinates, the center of the earth is the origin, and all coordinates are measured wrt the origin. The horizon point exists directly above the center of the earth,
   * (0, 0, REarth). The detector exists on the positive x-axis, (horizondist, 0, REarth). The intersection of the tau trajectory and the horizon-detector plane exists on the
   * same axis as the detector, (horizondist - l, 0, REarth). */
  
  Double_t *tau_earth_coords = GetTauEarthCoords(azi, elv, l, horizondist); //gets the emergent coordinates of tau from earth, km
  
  printf("tau earth coords: %.3f, %.3f, %.3f, %1.0f | %.3f\n", horizondist - tau_earth_coords[0], tau_earth_coords[1], tau_earth_coords[2] - REarth, tau_earth_coords[3], DistanceThroughEarth(azi, elv, l));
  if(tau_earth_coords[3] < 0) //check if the tau even intersects with the earth, if not, terminate
    {
      cout << "Tau trajectory does not intersect with Earth. No air shower will be created." << endl;
      return;
    }
  
  if(horizondist - tau_earth_coords[0] < vfovbelowdist && TauEmergeBelowVFoV(tau_earth_coords, vfovbelow, horizondist)) //check to see if the tau emerges in front of and out of view of the lower fov of the detector
    tau_config = 2;
  else if(tau_earth_coords[0] >= 0) //check if tau emerges in front of the horizon
    tau_config = 1;
  
  Double_t dprime, d, p_decay, p_decay_max = 0.9999999999999; //exclude P = 1 as that would yield infinity due to ln
  
  
  //calculating d' and d
  if(tau_config == 0) //beyond the horizon
    {
      dprime = sqrt((horizondist - l - tau_earth_coords[0]) * (horizondist - l - tau_earth_coords[0]) + tau_earth_coords[1] * tau_earth_coords[1] + 
		    (REarth - tau_earth_coords[2]) * (REarth - tau_earth_coords[2])); //distance between tau emergence and intersection point between tau trajectory and detector-horizon plane
      
      cout << "Tau emerges beyond the horizon." << endl;
      
      //~ cout << AirShowerAppearsInVFoV(tau_earth_coords, dprime, - 0.5 * dprime, azi, elv, l, vfovabove, vfovbelow, horizondist, s, tau_config) << endl; 
      //~ return;
    }
  else if (tau_config == 1) //in front of horizon and in view, tau immediately becomes within view of the detector (d' = 0)
    {
      dprime = 0;
      cout << "Tau emerges in front of horizon and within lower detector FoV." << endl;
    }
  else //tau emerges close to the detector and is initially out of view of the lower fov
    {
      Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovbelow, horizondist, true); //get the coordinates of the intersection point between the tau trajectory and the lower fov
      
      if(as_c[3] < 0) //trajectory does not intersect the lower fov in front of the detector
	{
	  cout << "Tau trajectory intersects lower FoV beyond detector. No air shower will be visible." << endl;
	  return;
	}
      
      dprime = sqrt((as_c[0] - tau_earth_coords[0]) * (as_c[0] - tau_earth_coords[0]) + (as_c[1] - tau_earth_coords[1]) * (as_c[1] - tau_earth_coords[1]) + 
		    (as_c[2] - tau_earth_coords[2]) * (as_c[2] - tau_earth_coords[2])); //distance between tau emergence and the point at which the tau enters the lower fov
      
      cout << "Tau emerges in front of horizon and out of lower detector FoV." << endl;
      
      //~ cout << AirShowerAppearsInVFoV(tau_earth_coords, dprime, 0.5 * dprime, azi, elv, l, vfovabove, vfovbelow, horizondist, s, tau_config) << endl; 
      //~ return;
    }
  
  Double_t D; //sum of d' and d
  
  for(int l = 0;;) //here we begin rolling for tau decay probabilities
    {
      p_decay = rand() / double(RAND_MAX) * p_decay_max; //get random number (0, 1)
      //~ p_decay = 1.71444 / 100.; //DEBUGGING PROBABILITY
      D = -1. * log(1. - p_decay) * DecayLength; //from the probability, calculate the tau's decay length
      d = D - dprime; //get d
      l++; //increment a counter for the number of rolls
      
      if(AirShowerAppearsInVFoV(tau_earth_coords, dprime, d, azi, elv, l, vfovabove, vfovbelow, horizondist, s, tau_config)) //check to see if the decay probability will yeild an observable shower
	{
	  cout << "Success probability: " << p_decay * 100. << "%" << endl;
	  cout << "Number of trials: " << l << endl;
	  break; //if the shower is observable, break the loop
	}
      
      cout << "Trial probability: " << p_decay * 100. << "%" << endl;
    }
  //~ printf("d: %f, dprime: %f, p_decay: %f\n", d, dprime, p_decay);
  //~ cout << tau_config << endl;
  
  cout << "Distance tau travels before entering FoV: " << dprime << " km" << endl;
  cout << "Distance tau travels before decaying: " << D << " km" << endl;
  if(d < 0)
    cout << "Tau decays before entering FoV." << endl;
  else
    cout << "Distance tau travels within FoV: " << d << " km" << endl;
  cout << "Tau decay probability: " << p_decay * 100. << "%" << endl;
  
  Double_t ang_proj_alt, ang_proj_hwidth, ang_proj_length; //angular altitude, width, and length of the shower within the detector
  Double_t image_start_x, image_start_y, image_end_x, image_end_y; //start and end coordinates for the shower within the camera
  Double_t i, j, k; //unit vector for the air shower
  Double_t shower_percent; //the amount of the shower starting from the base to the topmost observable point of the shower
  Double_t *as_base = (Double_t*) calloc(3, sizeof(Double_t)); //shower base coordinates, km
  Double_t *as_tip = (Double_t*) calloc(3, sizeof(Double_t)); //shower tip coordinates, km
  Double_t *as_actual_ang = (Double_t*) calloc(3, sizeof(Double_t)); //used for obtaining the PEs incident on each pixel
  
  i = cos(azi * pi / 180.) * cos(elv * pi / 180.); //getting the unit vector pointing in the direction of the air shower
  j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
  k = sin(elv * pi / 180.);
  
  as_base[0] = tau_earth_coords[0] + i * D; //getting the coordinates for the base of the air shower
  as_base[1] = tau_earth_coords[1] + j * D;
  as_base[2] = tau_earth_coords[2] + k * D;
  
  as_tip[0] = tau_earth_coords[0] + i * (D + s); //getting the coordinates for the tip of the air shower
  as_tip[1] = tau_earth_coords[1] + j * (D + s);
  as_tip[2] = tau_earth_coords[2] + k * (D + s);
  
  as_actual_ang[0] = 180. / pi * atan2(as_base[1], horizondist - as_base[0]) + hfov * 0.5; //the base x coordinate of the shower in the camera
  as_actual_ang[1] = 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]) + vfovbelow; //the base y coordinate of the shower in the camera
  
  printf("base: %f, %f, %f\n", horizondist - as_base[0], as_base[1], as_base[2] - REarth);
  printf("tip: %f, %f, %f\n", horizondist - as_tip[0], as_tip[1], as_tip[2] - REarth);
  printf("length to base: %f\n", sqrt((as_base[0] - horizondist) * (as_base[0] - horizondist) + (as_base[1]) * (as_base[1]) + (as_base[2] - REarth) * (as_base[2] - REarth)));
  
  bool *as_cutoff = ASCutoff(as_base, as_tip, vfovbelow, vfovabove, horizondist); //boolean array that denotes whether or not the base and/or tip of the shower is truncated
  // as_cutoff[0] = base, as_cutoff[1] = tip
  
  printf("base bool: %d, tip bool: %d\n", as_cutoff[0], as_cutoff[1]);
  
  //here we calculate and setup the size of the air shower within the camera
  if(!as_cutoff[0] && !as_cutoff[1]) //both base and tip are in fov, fully contained shower
    {
      //calculating the angular dimentions of the shower within the camera
      ang_proj_alt = fabs(180. / pi * atan2(as_tip[2] - REarth, horizondist - as_tip[0]) - 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]));
      ang_proj_hwidth = fabs(180. / pi * atan2(as_tip[1], horizondist - as_tip[0]) - 180. / pi * atan2(as_base[1], horizondist - as_base[0]));
      ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
      
      //setting the start and end coordinates of the shower within the camera; the middle of the camera is treated as the origin
      image_start_x = 180. / pi * atan2(as_base[1], horizondist - as_base[0]) + hfov * 0.5;
      image_start_y = 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]) + vfovbelow;
      image_end_x = image_start_x + ang_proj_hwidth;
      image_end_y = image_start_y + ang_proj_alt;
      
      //angular distance between the start and end of the shower in the camera
      as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
      //percentage of the shower from the base to the topmost visible point of the shower
      shower_percent = 1;
      
      cout << "Shower fully contained within vertical FoV." << endl;
    }
  else if(as_cutoff[0] && !as_cutoff[1]) //only base is cut off
    {
      if(as_base[0] < 0) //the base is truncated by the horizon
	{
	  //calculating the angular dimentions of the shower within the camera
	  ang_proj_alt = 180. / pi * atan2(as_tip[2] - REarth, horizondist - as_tip[0]);
	  ang_proj_hwidth = fabs(180. / pi * atan2(as_tip[1], horizondist - as_tip[0]));
	  
	  //setting the start and end coordinates of the shower within the camera; the middle of the camera is treated as the origin
	  image_start_x = hfov * 0.5;
	  image_start_y = vfovbelow;
	  image_end_x = image_start_x + ang_proj_hwidth;
	  image_end_y = image_start_y + ang_proj_alt;
	  
	  //angular distance between the start and end of the shower in the camera
	  as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	  //percentage of the shower from the base to the topmost visible point of the shower
	  shower_percent = 1;
	  
	  cout << "Shower base truncated by horizon." << endl;
	}
      else //the base is truncated by the lower vertical fov 
	{
	  Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovbelow, horizondist, true); //get the coordinates for the intersection point between the tau trajectory and the lower fov
	  
	  printf("lower fov : %f, %f, %f\n", horizondist - as_c[0], as_c[1], as_c[2] - REarth);
	  
	  //calculating the angular dimentions of the shower within the camera
	  ang_proj_alt = fabs(180. / pi * atan2(as_tip[2] - REarth, horizondist - as_tip[0]) - 180. / pi * atan2(as_c[2] - REarth, horizondist - as_c[0]));
	  ang_proj_hwidth = fabs(180. / pi * atan2(as_tip[1], horizondist - as_tip[0]) - 180. / pi * atan2(as_c[1], horizondist - as_c[0]));
	  
	  //setting the start and end coordinates of the shower within the camera; the middle of the camera is treated as the origin
	  image_start_x = 180. / pi * atan2(as_c[1], horizondist - as_c[0]) + hfov * 0.5;
	  image_start_y = 0;
	  image_end_x = image_start_x + ang_proj_hwidth;
	  image_end_y = image_start_y + ang_proj_alt;
	  
	  //angular distance between the start and end of the shower in the camera
	  as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	  //percentage of the shower from the base to the topmost visible point of the shower
	  shower_percent = 1;
	  
	  cout << "Shower base truncated by lower FoV." << endl;
	}
      
      ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth); //calculating the length of the shower within the camera
    }
  else if(!as_cutoff[0] && as_cutoff[1]) //only tip is cut off
    {
      Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovabove, horizondist, false); //get the coordinates for the intersection point between the tau trajectory and the upper fov
      
      //calculating the angular dimentions of the shower within the camera
      ang_proj_alt = fabs(180. / pi * atan2(as_c[2] - REarth, horizondist - as_c[0]) - 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]));
      ang_proj_hwidth = fabs(180. / pi * atan2(as_c[1], horizondist - as_c[0]) - 180. / pi * atan2(as_base[1], horizondist - as_base[0]));
      ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
      
      //setting the start and end coordinates of the shower within the camera; the middle of the camera is treated as the origin
      image_start_x = 180. / pi * atan2(as_base[1], horizondist - as_base[0]) + hfov * 0.5;
      image_start_y = 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]) + vfovbelow;
      image_end_x = image_start_x + ang_proj_hwidth;
      image_end_y = image_start_y + ang_proj_alt;
      
      //angular distance between the start and end of the shower in the camera
      as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
      //percentage of the shower from the base to the topmost visible point of the shower
      shower_percent = sqrt((as_base[0] - as_c[0]) * (as_base[0] - as_c[0]) + (as_base[1] - as_c[1]) * (as_base[1] - as_c[1]) + (as_base[2] - as_c[2]) * (as_base[2] - as_c[2])) / s;
      //~ printf("percent of shower in fov: %f\n", shower_percent * 100.);
      //~ printf("pes: %f\n", gh->Integral(X_0, shower_percent * (X_max - X_0) + X_0));
      cout << "Shower tip truncated by upper FoV." << endl;
    }
  else //both base and tip are cut off
    {
      if(as_base[0] < 0) //the base is truncated by the horizon, the tip is truncated by the upper fov
	{
	  Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovabove, horizondist, false); //get the coordinates for the intersection point between the tau trajectory and the upper fov
	  
	  //calculating the angular dimentions of the shower within the camera
	  ang_proj_alt = 180. / pi * atan2(as_c[2] - REarth, horizondist - as_c[0]);	
	  ang_proj_hwidth = fabs(180. / pi * atan2(as_c[1], horizondist - as_c[0]));
	  
	  //setting the start and end coordinates of the shower within the camera; the middle of the camera is treated as the origin
	  image_start_x = hfov * 0.5;
	  image_start_y = vfovbelow;
	  image_end_x = image_start_x + ang_proj_hwidth;
	  image_end_y = image_start_y + ang_proj_alt;
	  
	  //angular distance between the start and end of the shower in the camera
	  as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	  //percentage of the shower from the base to the topmost visible point of the shower
	  shower_percent = sqrt((as_base[0] - as_c[0]) * (as_base[0] - as_c[0]) + (as_base[1] - as_c[1]) * (as_base[1] - as_c[1]) + (as_base[2] - as_c[2]) * (as_base[2] - as_c[2])) / s;
	  
	  cout << "Shower base truncated by horizon, shower tip truncated by upper FoV" << endl;
	}
      else //the base and tip are both truncated by the lower and upper fov's, respectively (the alt should just be lower vfov + upper vfov)
	{
	  Double_t *as = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovbelow, horizondist, true); //get the coordinates for the intersection point between the tau trajectory and the lower fov
	  
	  Double_t as_lower[3];
	  
	  as_lower[0] = as[0]; //load coordinates into stack memory
	  as_lower[1] = as[1];
	  as_lower[2] = as[2];
	  
	  as = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovabove, horizondist, false); //get the coordinates for the intersection point between the tau trajectory and the upper fov
	  
	  Double_t as_upper[3];
	  
	  as_upper[0] = as[0]; //load coordinates into stack memory
	  as_upper[1] = as[1];
	  as_upper[2] = as[2];
	  
	  //calculating the angular dimentions of the shower within the camera
	  ang_proj_alt = fabs(180. / pi * atan2(as_upper[2] - REarth, horizondist - as_upper[0]) - 180. / pi * atan2(as_lower[2] - REarth, horizondist - as_lower[0]));
	  ang_proj_hwidth = fabs(180. / pi * atan2(as_upper[1], horizondist - as_upper[0]) - 180. / pi * atan2(as_lower[1], horizondist - as_lower[0]));
	  
	  //setting the start and end coordinates of the shower within the camera; the middle of the camera is treated as the origin
	  image_start_x = 180. / pi * atan2(as_lower[1], horizondist - as_lower[0]) + hfov * 0.5;
	  image_start_y = 180. / pi * atan2(as_lower[2] - REarth, horizondist - as_lower[0]) + vfovbelow;
	  image_end_x = image_start_x + ang_proj_hwidth;
	  image_end_y = image_start_y + ang_proj_alt;
	  
	  //angular distance between the start and end of the shower in the camera
	  as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	  //percentage of the shower from the base to the topmost visible point of the shower
	  shower_percent = sqrt((as_base[0] - as_upper[0]) * (as_base[0] - as_upper[0]) + (as_base[1] - as_upper[1]) * (as_base[1] - as_upper[1]) + (as_base[2] - as_upper[2]) * (as_base[2] - as_upper[2])) / s;
	  
	  cout << "Shower base truncated by lower FoV, shower tip truncated by upper FoV" << endl;
	}
      
      ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth); //calculating the length of the shower within the camera
    }
  
  cout << "Projected angular altitude: " << ang_proj_alt << " deg" << endl;
  cout << "Projected angular width: " << ang_proj_hwidth << " deg" << endl;
  cout << "Projected angular length: " << ang_proj_length << " deg" << endl;
  
  
  printf("image start: (%f, %f) | image end: (%f, %f)\n", image_start_x, image_start_y, image_end_x, image_end_y);
  printf("h: %f | w: %f\n", image_end_y - image_start_y, image_end_x - image_start_x);
  
  Double_t showerslope = ang_proj_alt / ang_proj_hwidth; //calculating the slope of the air shower within the camera
  Double_t showerintercepty = image_start_y - showerslope * image_start_x; //calculating the y-intercept of the shower within the camera
  cout << "Shower formula: y = " << showerslope << "x + " << showerintercepty << endl;
  
  Double_t triggerthresh = 0.999999;
  Double_t oversample = 10; //amount of internal oversampling
  
  TString camtitle;
  camtitle.Form("Projected Air Shower Image In Camera (%.0f deg azimuth, %.0f deg elevation, %.0f km distance) with %.2f%% decay probability, %.1E GeV tau energy", azi, elv, l, p_decay * 100., E);
  
  TCanvas *pec = new TCanvas(camtitle,camtitle,1700,800);
  pec->cd(1);
  pec->SetRightMargin(0.12);
  
  TH2F *camera = new TH2F("camera", camtitle, nxpx, 0, nxpx, nypx, 0, nypx); //regular camera histogram
  camera->GetXaxis()->SetTitle("Pixel x");
  camera->GetYaxis()->SetTitle("Pixel y");
  camera->GetZaxis()->SetTitle("PE count");
  
  TH2F *cameraOS = new TH2F("camOS", "camOS", nxpx * oversample, 0, nxpx, nypx * oversample, 0, nypx); //oversampled camera histogram
  
  TRandom *r = new TRandom();
  
  Double_t xint_a, xint_b, x, y, totalcamPE = 0;
  Double_t xrange[2];
  Double_t int_range[2];
  
  //here we find what pixels the air shower intersects with in the camera
  for(int xbin = 1 ; xbin <= cameraOS->GetNbinsX(); xbin++) //loop through all the bins of the oversampled histogram
    {
      for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++) 
	{
	  xint_a = xint_b = -1; //the x-values for the intersection points that each pixel has with the air shower
	  
	  //check if the bottom of the pixel intersects with the air shower
	  y = (ybin - 1) * pxheight / oversample;
	  x = (y - showerintercepty) / showerslope;
	  if(x >= (xbin - 1) * pxwidth / oversample && x <= xbin * pxwidth / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
	    xint_a = x;
	  
	  //check if the top of the pixel intersects with the air shower
	  y = ybin * pxheight / oversample;
	  x = (y - showerintercepty) / showerslope;
	  if(x >= (xbin - 1) * pxwidth / oversample && x <= xbin * pxwidth / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
	    {
	      if(xint_a <= 0)
		xint_a = x;
	      else
		xint_b = x;
	    }
	  
	  //check if the left side of the pixel intersects with the air shower
	  x = (xbin - 1) * pxwidth / oversample;
	  y = showerslope * x + showerintercepty;
	  if(y >= (ybin - 1) * pxheight / oversample && y <= ybin * pxheight / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
	    {
	      if(xint_a <= 0)
		xint_a = x;
	      else
		xint_b = x;
	    }
	  
	  //check if the right of the pixel intersects with the air shower
	  x = xbin * pxwidth / oversample;
	  y = showerslope * x + showerintercepty;
	  if(y >= (ybin - 1) * pxheight / oversample && y <= ybin * pxheight / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
	    {
	      if(xint_a <= 0)
		xint_a = x;
	      else
		xint_b = x;
	    }
	  
	  //~ cout << xint_a << ", " << xint_b << " | " << xbin << ", " << ybin << endl;
	  
	  //checking the number of crossings that the air shower makes in a pixel (1 or 2)
	  if(xint_a > xint_b && xint_b < 0) //check to see if there is only a single crossing in the pixel
	    {
	      //~ cout << "Single crossing at x = " << xint_a << " for bin (" << xbin << ", " << ybin << ")" << endl;
	      if(abs(image_start_x - xint_a) < abs(image_end_x - xint_a)) //the single crossing starts from the base of the shower
		{
		  //the x range in camera coordinates for which the integration will be performed over
		  xrange[0] = image_start_x;
		  xrange[1] = xint_a;
		  //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
		  
		  //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
		  int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
				      ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
				      ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  
		  //integrate to get the average number of PEs and then distribute via the Poisson distribution
		  Double_t PEs = gh->Integral(int_range[0], int_range[1]);
		  Double_t poiPEs = r->PoissonD(PEs);
		  
		  //fill both the regular and oversampled histograms
		  cameraOS->SetBinContent(xbin, ybin, poiPEs);
		  camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), poiPEs);
		}
	      else //the single crossing ends at the tip of the shower
		{
		  //the x range in camera coordinates for which the integration will be performed over
		  xrange[0] = xint_a;
		  xrange[1] = image_end_x;
		  //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
		  
		  //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
		  int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
				      ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
				      ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  
		  //integrate to get the average number of PEs and then distribute via the Poisson distribution
		  Double_t PEs = gh->Integral(int_range[0], int_range[1]);
		  Double_t poiPEs = r->PoissonD(PEs);
		  
		  //fill both the regular and oversampled histograms
		  cameraOS->SetBinContent(xbin, ybin, poiPEs);
		  camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), poiPEs);
		}
	    }
	  else if(xint_a > xint_b) //now to consder the double crossings of the air shower through a pixel
	    {
	      //~ cout << "Double crossing at x1 = " << xint_b << " x2 = " << xint_a << " for bin (" << xbin << ", " << ybin << ")" << endl;
	      
	      //the x range in camera coordinates for which the integration will be performed over
	      xrange[0] = xint_b;
	      xrange[1] = xint_a;
	      //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
	      
	      //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
	      int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
				  ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
	      int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
				  ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
	      
	      //integrate to get the average number of PEs and then distribute via the Poisson distribution
	      Double_t PEs = gh->Integral(int_range[0], int_range[1]);
	      Double_t poiPEs = r->PoissonD(PEs);
	      
	      //fill both the regular and oversampled histograms	
	      cameraOS->SetBinContent(xbin, ybin, poiPEs);
	      camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), poiPEs);
	    }
	  else if(xint_a < xint_b)
	    {
	      //~ cout << "Double crossing at x1 = " << xint_a << " x2 = " << xint_b << " for bin (" << xbin << ", " << ybin << ")" << endl;
	      
	      //the x range in camera coordinates for which the integration will be performed over
	      xrange[0] = xint_a;
	      xrange[1] = xint_b;
	      //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
	      
	      //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
	      int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
				  ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
	      int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
				  ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
	      
	      //integrate to get the average number of PEs and then distribute via the Poisson distribution
	      Double_t PEs = gh->Integral(int_range[0], int_range[1]);
	      Double_t poiPEs = r->PoissonD(PEs);
	      
	      //fill both the regular and oversampled histograms	
	      cameraOS->SetBinContent(xbin, ybin, poiPEs);
	      camera->Fill(cameraOS->GetXaxis()->GetBinCenter(xbin), cameraOS->GetYaxis()->GetBinCenter(ybin), poiPEs);
	    }
	}
    }
  
  TH2F *cameraOSLat = new TH2F("camOSLat", "camOSLat", nxpx * oversample, 0, nxpx, nypx * oversample, 0, nypx); //creat a new oversampled histogram for the lateral distribution of the shower
  
  Double_t latsigma = 0; //the stdev for the lateral gaussian (to be properly implemented)
  
  for(int xbin = 1; xbin <= cameraOS->GetNbinsX(); xbin++) //loop through all the bins of the oversampled histogram
    {
      for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++)
	{
	  if(cameraOS->GetBinContent(xbin, ybin) > 0) //if the bin isn't empty
	    {
	      int numPEs = (int) cameraOS->GetBinContent(xbin, ybin); //store the number of PEs in that bin
	      Double_t dist;
	      Double_t latslope = -1. / showerslope; //lateral distribution is normal to the direction of the shower in the camera
	      Double_t angle = atan(latslope); //get the normal angle
	      
	      for(int i = 0; i < numPEs; i++) //loop through the number of PEs within that bin
		{
		  dist = r->Gaus(0, latsigma); //get a distance away from the pixel according to the gaussian distribution
		  
		  //decompose into x and y components
		  Double_t x = dist * cos(angle) + xbin;
		  Double_t y = dist * sin(angle) + ybin;
		  
		  if((x >= 1 && x <= cameraOS->GetNbinsX()) && (y >= 1 && y <= cameraOS->GetNbinsY())) //if the distibuted distance is within the camera
		    {
		      cameraOSLat->Fill(cameraOS->GetXaxis()->GetBinCenter(x), cameraOS->GetYaxis()->GetBinCenter(y), 1); //fill the bin with a single PE
		    }
		}
	    }
	}
    }
  
  //get the total PEs expected to arrive from the shower
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
  
  //if there is less than a single PE captured, terminate
  if(totalcamPE < 1)
    {
      cout << "Less than a single PE captured. No air shower can be displayed." << endl;
      exit(EXIT_SUCCESS);
    }
  
  cout.precision(17);
  cout << "Average total PEs produced: " << (int) totalPEs << endl;
  cout << "Total captured PEs: " << (int) totalcamPE << endl;
  
  camera->SetContour(99);
  camera->SetStats(0);
  camera->Draw("COLZ");
  
  free(as_base);
  free(as_tip);
  free(as_actual_ang);
  
  Double_t blursigma = 1 * oversample; //sigma for the point spread function of the camera
  
  TString camtitleb;
  camtitleb.Form("Projected Air Shower Image In Camera With %1.1f deg Gaussian Blur (%.0f deg azimuth, %.0f deg elevation, %.0f km distance) with %.2f%% decay probability, %.1E GeV tau energy", pxwidth * blursigma / oversample, azi, elv, l, p_decay * 100., E);
  
  TCanvas *pecBlur = new TCanvas(camtitleb,camtitleb,1200,1200);
  pecBlur->cd(1);
  pecBlur->SetRightMargin(0.12);
  
  TH2F *cameraBlur = new TH2F(camtitleb, camtitleb, nxpx, 0, nxpx, nypx, 0, nypx); //blured histogram of the shower
  cameraBlur->GetXaxis()->SetTitle("Pixel x");
  cameraBlur->GetYaxis()->SetTitle("Pixel y");
  cameraBlur->GetZaxis()->SetTitle("PE count");
  
  //~ Double_t kernrad = 6 * blursigma; //sigmas
  //~ Double_t GMatrix[(int)((1 + 2 * kernrad))][(int)((1 + 2 * kernrad))];
  
  //~ TF2 *G2D = new TF2("g2d", Gauss2D, -1. * kernrad, kernrad, -1. * kernrad, kernrad , 1);
  //~ G2D->SetParameter(0, blursigma);
  //~ Double_t gweight = 0;
  
  //~ for(int i = 0; i < (int)(1 + 2 * kernrad); i++)
  //~ {
  //~ for(int j = 0; j < (int)(1 + 2 * kernrad); j++)
  //~ {
  //~ GMatrix[i][j] = G2D->Eval(j - kernrad, kernrad - i);
  //~ gweight += G2D->Eval(j - kernrad, kernrad - i);
  //~ }
  //~ }
  
  //~ for(int i = 0; i < (int)(1 + 2 * kernrad); i++)
  //~ {
  //~ for(int j = 0; j < (int)(1 + 2 * kernrad); j++)
  //~ {
  //~ GMatrix[i][j] = GMatrix[i][j] / gweight;
  //~ }
  //~ }
  
  //~ for(int xbin = 1; xbin <= cameraOS->GetNbinsX(); xbin++)
  //~ {
  //~ for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++)
  //~ {
  //~ if(cameraOS->GetBinContent(xbin, ybin) > 0)
  //~ {
  //~ for(int x = xbin - kernrad; x <= xbin + kernrad; x++)
  //~ {
  //~ for(int y = ybin - kernrad; y <= ybin + kernrad; y++)
  //~ {
  //~ if((x >= 1 && x <= cameraOS->GetNbinsX()) && (y >= 1 && y <= cameraOS->GetNbinsY()))
  //~ {	
  //~ cameraBlur->Fill(cameraOS->GetXaxis()->GetBinCenter(x), cameraOS->GetYaxis()->GetBinCenter(y), 
  //~ GMatrix[(int)(y - ybin + kernrad)][(int)(x - xbin + kernrad)] * cameraOS->GetBinContent(xbin, ybin));
  //~ }
  //~ }
  //~ }
  //~ }
  //~ }
  //~ }
  
  for(int xbin = 1; xbin <= cameraOSLat->GetNbinsX(); xbin++) //loop through the bins of the oversampled and laterally distributed histogram
    {
      for(int ybin = 1; ybin <= cameraOSLat->GetNbinsY(); ybin++)
	{
	  if(cameraOSLat->GetBinContent(xbin, ybin) > 0) //if the bin is not empty
	    {
	      int numPEs = (int) cameraOSLat->GetBinContent(xbin, ybin); //get the number of PEs in that bin
	      Double_t x, y;
	      
	      for(int i = 0; i < numPEs; i++) //loop through the number of PEs
		{
		  x = r->Gaus(0, blursigma) + xbin; //generate a random x distance away from the bin according to a gaussian
		  y = r->Gaus(0, blursigma) + ybin; //generate a random y distance away from the bin according to a gaussian
		  
		  if((x >= 1 && x <= cameraOSLat->GetNbinsX()) && (y >= 1 && y <= cameraOSLat->GetNbinsY()))
		    {
		      cameraBlur->Fill(cameraOSLat->GetXaxis()->GetBinCenter(x), cameraOSLat->GetYaxis()->GetBinCenter(y), 1); //fill the bin with a single PE
		    }
		}
	    }
	}
    }
  
  totalcamPE = 0;
  
  //get the total of PEs after applying point spread function
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
  
  //formatting the aspect ratio of the histogram
  
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

void SimNAirShowers(Double_t azi, Double_t elv, Double_t l, Double_t E, int trials)
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
  int det_alt = 2;
  
  Double_t totalPEs = myPEfunction2(azi * TMath::DegToRad(), elv * TMath::DegToRad(), l) * eff_mirror_size * E;
  if(totalPEs == 0)
    {
      cout << "Azimuth angle too large." << endl;
      return;
    }
  
  TF1 *gh = new TF1("gh", GH, X_0, X_max, 4);
  gh->SetParameters(X_0, X_max, lmbda, N_0);
  Double_t norm = 1.0 / gh->Integral(X_0, X_max);
  N_0 = norm * totalPEs;
  
  Double_t hfov = 60; //deg
  Double_t vfov = 5; //deg
  Double_t pxwidth = 0.3; //deg
  Double_t pxheight = pxwidth; //deg
  Double_t vfovbelow = 3.0; //deg
  Double_t vfovabove = 2.0; //deg
  
  Double_t nxpx = ceil(hfov / pxwidth);
  Double_t nypx = ceil(vfov / pxheight);
  
  Double_t oversample = 10;
  
  Double_t horizondist = sqrt((DetectorAltitude[iConfig] + REarth) * (DetectorAltitude[iConfig] + REarth) - REarth * REarth); //km
  Double_t vfovbelowdist = DistanceFovBelow(vfovbelow, horizondist); //km
  Double_t DecayLength = E * c * DecayTime / Mtau; //km
  
  short tau_config = 0; //0 = beyond the horizon, 1 = in front of horizon but beyond detector fov, 2 = in front of horizon and in front of detector fov
  
  Double_t *tau_earth_coords = GetTauEarthCoords(azi, elv, l, horizondist);
  //~ printf("tau earth coords: %.3f, %.3f, %.3f, %1.0f | %.3f\n", horizondist - tau_earth_coords[0], tau_earth_coords[1], tau_earth_coords[2] - REarth, tau_earth_coords[3], DistanceThroughEarth(azi, elv, l));
  if(tau_earth_coords[3] < 0)
    {
      cout << "Tau trajectory does not intersect with Earth. No air shower will be created." << endl;
      return;
    }
  
  if(horizondist - tau_earth_coords[0] < vfovbelowdist && TauEmergeBelowVFoV(tau_earth_coords, vfovbelow, horizondist))
    tau_config = 2;
  else if(tau_earth_coords[0] >= 0)
    tau_config = 1;
  
  Double_t dprime, d, p_decay, p_decay_max = 0.9999999999999;
  
  if(tau_config == 0) //beyond the horizon
    {
      dprime = sqrt((horizondist - l - tau_earth_coords[0]) * (horizondist - l - tau_earth_coords[0]) + tau_earth_coords[1] * tau_earth_coords[1] + 
		    (REarth - tau_earth_coords[2]) * (REarth - tau_earth_coords[2]));
      
      //~ cout << "Tau emerges beyond the horizon." << endl;
      
      //~ cout << AirShowerAppearsInVFoV(tau_earth_coords, dprime, - 0.5 * dprime, azi, elv, l, vfovabove, vfovbelow, horizondist, s, tau_config) << endl; 
      //~ return;
    }
  else if (tau_config == 1) //in front of horizon and in view (d' = 0)
    {
      dprime = 0;
      //~ cout << "Tau emerges in front of horizon and within lower detector FoV." << endl;
    }
  else //in front of horizon and out of view
    {
      Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovbelow, horizondist, true);
      
      if(as_c[3] < 0)
	{
	  cout << "Tau trajectory intersects lower FoV beyond detector. No air shower will be visible." << endl;
	  return;
	}
      
      dprime = sqrt((as_c[0] - tau_earth_coords[0]) * (as_c[0] - tau_earth_coords[0]) + (as_c[1] - tau_earth_coords[1]) * (as_c[1] - tau_earth_coords[1]) + 
		    (as_c[2] - tau_earth_coords[2]) * (as_c[2] - tau_earth_coords[2]));
      
      //~ cout << "Tau emerges in front of horizon and out of lower detector FoV." << endl;
      
      //~ cout << AirShowerAppearsInVFoV(tau_earth_coords, dprime, 0.5 * dprime, azi, elv, l, vfovabove, vfovbelow, horizondist, s, tau_config) << endl; 
      //~ return;
    }
  
  bool is_visible, base_trunc, tip_trunc, fully_contained;
  Double_t D, totalCapturedPEs;
  int count, success = 0, full = 0;
  
  TH2F *cameraOS = new TH2F("camos", "camos", nxpx * oversample, 0, nxpx, nypx * oversample, 0, nypx);
  TH2F *cameraOSLat = new TH2F("camoslat", "camoslat", nxpx * oversample, 0, nxpx, nypx * oversample, 0, nypx);
  TH2F *cameraBlur = new TH2F("camblur", "camblur", nxpx, 0, nxpx, nypx, 0, nypx);
  
  TRandom *r = new TRandom();
  
  TString filetitle;
  filetitle.Form("%d_%.0f_%.0f_%.0f_%.0e_as.root", trials, azi, elv, l, E);
  TFile *f = new TFile(filetitle, "RECREATE");
  
  TTree *ASDataTree = new TTree("ASData", "Air Shower Trials");
  ASDataTree->Branch("trial_number", &count, "trials/I"); //
  ASDataTree->Branch("decay_prob", &p_decay, "p_decay/D"); //
  ASDataTree->Branch("decay_length", &D, "D/D"); //
  ASDataTree->Branch("is_visible", &is_visible, "is_visible/O"); //
  ASDataTree->Branch("fully_contained", &fully_contained, "fully_contained/O"); //
  ASDataTree->Branch("base_truncated", &base_trunc, "base_trunc/O"); //
  ASDataTree->Branch("tip_truncated", &tip_trunc, "tip_trunc/O"); //
  ASDataTree->Branch("total_captured_PEs", &totalCapturedPEs, "totalCapturedPEs/D"); //
  
  ASDataTree->Branch("camera_hist", "TH2F", &cameraBlur);
  
  for(count = 0; count < trials;)
    {
      count++;
      p_decay = rand() / double(RAND_MAX) * p_decay_max;
      //~ p_decay = 3.91597 / 100.; //DEBUGGING PROBABILITY
      D = -1. * log(1. - p_decay) * DecayLength;
      d = D - dprime;
      
      cameraOS->Reset();
      cameraOSLat->Reset();
      cameraBlur->Reset();
      
      if(!AirShowerAppearsInVFoV(tau_earth_coords, dprime, d, azi, elv, l, vfovabove, vfovbelow, horizondist, s, tau_config))
	{
	  is_visible = false;
	  fully_contained = false;
	  tip_trunc = true;
	  base_trunc = true;
	  totalCapturedPEs = 0;
	  
	  ASDataTree->Fill();
	  continue;
	}
      
      is_visible = true;
      fully_contained = false;
      success++;
      
      Double_t ang_proj_alt, ang_proj_hwidth, ang_proj_length;
      Double_t image_start_x, image_start_y, image_end_x, image_end_y;
      Double_t i, j, k;
      Double_t shower_percent;
      Double_t as_base[3];
      Double_t as_tip[3];
      Double_t as_actual_ang[3];
      
      i = cos(azi * pi / 180.) * cos(elv * pi / 180.); //getting the unit vector pointing in the direction of the air shower
      j = sin(azi * pi / 180.) * cos(elv * pi / 180.);
      k = sin(elv * pi / 180.);
      
      as_base[0] = tau_earth_coords[0] + i * D;
      as_base[1] = tau_earth_coords[1] + j * D;
      as_base[2] = tau_earth_coords[2] + k * D;
      
      as_tip[0] = tau_earth_coords[0] + i * (D + s);
      as_tip[1] = tau_earth_coords[1] + j * (D + s);
      as_tip[2] = tau_earth_coords[2] + k * (D + s);
      
      as_actual_ang[0] = 180. / pi * atan2(as_base[1], horizondist - as_base[0]) + hfov * 0.5; //the base x coordinate of the shower in the camera
      as_actual_ang[1] = 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]) + vfovbelow; //the base y coordinate of the shower in the camera
      
      bool *as_cutoff = ASCutoff(as_base, as_tip, vfovbelow, vfovabove, horizondist);
      base_trunc = as_cutoff[0];
      tip_trunc = as_cutoff[1];
      
      if(!as_cutoff[0] && !as_cutoff[1]) //both base and tip are in fov
	{	
	  fully_contained = true;
	  full++;
	  
	  ang_proj_alt = fabs(180. / pi * atan2(as_tip[2] - REarth, horizondist - as_tip[0]) - 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]));
	  ang_proj_hwidth = fabs(180. / pi * atan2(as_tip[1], horizondist - as_tip[0]) - 180. / pi * atan2(as_base[1], horizondist - as_base[0]));
	  ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
	  
	  image_start_x = 180. / pi * atan2(as_base[1], horizondist - as_base[0]) + hfov * 0.5;
	  image_start_y = 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]) + vfovbelow;
	  image_end_x = image_start_x + ang_proj_hwidth;
	  image_end_y = image_start_y + ang_proj_alt;
	  
	  as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	  shower_percent = 1;
	  
	  printf("Fully contained for trial %d\n", count + 1);
	}
      else if(as_cutoff[0] && !as_cutoff[1]) //only base is cut off
	{
	  if(as_base[0] < 0) //beyond horizon
	    {
	      ang_proj_alt = 180. / pi * atan2(as_tip[2] - REarth, horizondist - as_tip[0]);
	      ang_proj_hwidth = fabs(180. / pi * atan2(as_tip[1], horizondist - as_tip[0]));
	      
	      image_start_x = hfov * 0.5;
	      image_start_y = vfovbelow;
	      image_end_x = image_start_x + ang_proj_hwidth;
	      image_end_y = image_start_y + ang_proj_alt;
	      as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	      
	      shower_percent = 1;
	      
	      printf("Base truncated by horizon for trial %d\n", count + 1);
	    }
	  else //below lower vfov
	    {
	      Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovbelow, horizondist, true); //get the coordinates for the intersection point between the tau trajectory and the lower fov
	      
	      ang_proj_alt = fabs(180. / pi * atan2(as_tip[2] - REarth, horizondist - as_tip[0]) - 180. / pi * atan2(as_c[2] - REarth, horizondist - as_c[0]));
	      ang_proj_hwidth = fabs(180. / pi * atan2(as_tip[1], horizondist - as_tip[0]) - 180. / pi * atan2(as_c[1], horizondist - as_c[0]));
	      
	      image_start_x = 180. / pi * atan2(as_c[1], horizondist - as_c[0]) + hfov * 0.5;
	      image_start_y = 0;
	      image_end_x = image_start_x + ang_proj_hwidth;
	      image_end_y = image_start_y + ang_proj_alt;
	      as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	      
	      shower_percent = 1;
	      
	      printf("Based truncated by lower FoV for trial %d\n", count + 1);
	    }
	  ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
	}
      else if(!as_cutoff[0] && as_cutoff[1]) //only tip is cut off
	{
	  Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovabove, horizondist, false);
	  
	  ang_proj_alt = fabs(180. / pi * atan2(as_c[2] - REarth, horizondist - as_c[0]) - 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]));
	  ang_proj_hwidth = fabs(180. / pi * atan2(as_c[1], horizondist - as_c[0]) - 180. / pi * atan2(as_base[1], horizondist - as_base[0]));
	  ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
	  
	  image_start_x = 180. / pi * atan2(as_base[1], horizondist - as_base[0]) + hfov * 0.5;
	  image_start_y = 180. / pi * atan2(as_base[2] - REarth, horizondist - as_base[0]) + vfovbelow;
	  image_end_x = image_start_x + ang_proj_hwidth;
	  image_end_y = image_start_y + ang_proj_alt;
	  
	  as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	  shower_percent = sqrt((as_base[0] - as_c[0]) * (as_base[0] - as_c[0]) + (as_base[1] - as_c[1]) * (as_base[1] - as_c[1]) + (as_base[2] - as_c[2]) * (as_base[2] - as_c[2])) / s;
	  
	  printf("Tip truncated by upper FoV for trial %d\n", count + 1);
	}
      else //both base and tip are cut off
	{
	  if(as_base[0] < 0) //beyond horizon
	    {
	      Double_t *as_c = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovabove, horizondist, false);
	      ang_proj_alt = 180. / pi * atan2(as_c[2] - REarth, horizondist - as_c[0]);	
	      ang_proj_hwidth = fabs(180. / pi * atan2(as_c[1], horizondist - as_c[0]));
	      
	      image_start_x = hfov * 0.5;
	      image_start_y = vfovbelow;
	      image_end_x = image_start_x + ang_proj_hwidth;
	      image_end_y = image_start_y + ang_proj_alt;
	      
	      as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	      shower_percent = sqrt((as_base[0] - as_c[0]) * (as_base[0] - as_c[0]) + (as_base[1] - as_c[1]) * (as_base[1] - as_c[1]) + (as_base[2] - as_c[2]) * (as_base[2] - as_c[2])) / s;
	      
	      printf("Base truncated by horizon, tip truncated by upper FoV for trial %d\n", count + 1);
	    }
	  else //really close to the detector (the alt should just be lower vfov + upper vfov)
	    {
	      Double_t *as = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovbelow, horizondist, true);
	      
	      Double_t as_lower[3];
	      
	      as_lower[0] = as[0]; //load coordinates into stack memory
	      as_lower[1] = as[1];
	      as_lower[2] = as[2];
	      
	      as = TauTrajectoryFoV(tau_earth_coords, azi, elv, vfovabove, horizondist, false);
	      
	      Double_t as_upper[3];
	      
	      as_upper[0] = as[0]; //load coordinates into stack memory
	      as_upper[1] = as[1];
	      as_upper[2] = as[2];
	      
	      ang_proj_alt = fabs(180. / pi * atan2(as_upper[2] - REarth, horizondist - as_upper[0]) - 180. / pi * atan2(as_lower[2] - REarth, horizondist - as_lower[0]));
	      ang_proj_hwidth = fabs(180. / pi * atan2(as_upper[1], horizondist - as_upper[0]) - 180. / pi * atan2(as_lower[1], horizondist - as_lower[0]));
	      
	      image_start_x = 180. / pi * atan2(as_lower[1], horizondist - as_lower[0]) + hfov * 0.5;
	      image_start_y = 180. / pi * atan2(as_lower[2] - REarth, horizondist - as_lower[0]) + vfovbelow;
	      image_end_x = image_start_x + ang_proj_hwidth;
	      image_end_y = image_start_y + ang_proj_alt;
	      
	      as_actual_ang[2] = sqrt((image_end_x - as_actual_ang[0]) * (image_end_x - as_actual_ang[0]) + (image_end_y - as_actual_ang[1]) * (image_end_y - as_actual_ang[1]));
	      shower_percent = sqrt((as_base[0] - as_upper[0]) * (as_base[0] - as_upper[0]) + (as_base[1] - as_upper[1]) * (as_base[1] - as_upper[1]) + (as_base[2] - as_upper[2]) * (as_base[2] - as_upper[2])) / s;
	      
	      printf("Base truncated by lower FoV, tip truncated by upper FoV for trial %d\n", count + 1);
	    }
	  ang_proj_length = sqrt(ang_proj_alt * ang_proj_alt + ang_proj_hwidth * ang_proj_hwidth);
	}
      
      Double_t showerslope = ang_proj_alt / ang_proj_hwidth; //calculating the slope of the air shower within the camera
      Double_t showerintercepty = image_start_y - showerslope * image_start_x;
      
      Double_t triggerthresh = 0.999999;
      
      Double_t xint_a, xint_b, x, y, totalcamPE = 0;
      Double_t xrange[2];
      Double_t int_range[2];
      
      for(int xbin = 1 ; xbin <= cameraOS->GetNbinsX(); xbin++) //loop through all the bins of the oversampled histogram
	{
	  for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++) 
	    {
	      xint_a = xint_b = -1; //the x-values for the intersection points that each pixel has with the air shower
	      
	      //check if the bottom of the pixel intersects with the air shower
	      y = (ybin - 1) * pxheight / oversample;
	      x = (y - showerintercepty) / showerslope;
	      if(x >= (xbin - 1) * pxwidth / oversample && x <= xbin * pxwidth / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
		xint_a = x;
	      
	      //check if the top of the pixel intersects with the air shower
	      y = ybin * pxheight / oversample;
	      x = (y - showerintercepty) / showerslope;
	      if(x >= (xbin - 1) * pxwidth / oversample && x <= xbin * pxwidth / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
		{
		  if(xint_a <= 0)
		    xint_a = x;
		  else
		    xint_b = x;
		}
	      
	      //check if the left side of the pixel intersects with the air shower
	      x = (xbin - 1) * pxwidth / oversample;
	      y = showerslope * x + showerintercepty;
	      if(y >= (ybin - 1) * pxheight / oversample && y <= ybin * pxheight / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
		{
		  if(xint_a <= 0)
		    xint_a = x;
		  else
		    xint_b = x;
		}
	      
	      //check if the right of the pixel intersects with the air shower
	      x = xbin * pxwidth / oversample;
	      y = showerslope * x + showerintercepty;
	      if(y >= (ybin - 1) * pxheight / oversample && y <= ybin * pxheight / oversample && ((x >= image_start_x && x <= image_end_x) || (x <= image_start_x && x >= image_end_x)))
		{
		  if(xint_a <= 0)
		    xint_a = x;
		  else
		    xint_b = x;
		}
	      
	      //~ cout << xint_a << ", " << xint_b << " | " << xbin << ", " << ybin << endl;
	      
	      //checking the number of crossings that the air shower makes in a pixel (1 or 2)
	      if(xint_a > xint_b && xint_b < 0) //check to see if there is only a single crossing in the pixel
		{
		  //~ cout << "Single crossing at x = " << xint_a << " for bin (" << xbin << ", " << ybin << ")" << endl;
		  if(abs(image_start_x - xint_a) < abs(image_end_x - xint_a)) //the single crossing starts from the base of the shower
		    {
		      //the x range in camera coordinates for which the integration will be performed over
		      xrange[0] = image_start_x;
		      xrange[1] = xint_a;
		      //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
		      
		      //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
		      int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
					  ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		      int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
					  ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		      
		      //integrate to get the average number of PEs and then distribute via the Poisson distribution
		      Double_t PEs = gh->Integral(int_range[0], int_range[1]);
		      Double_t poiPEs = r->PoissonD(PEs);
		      
		      //fill both the regular and oversampled histograms
		      cameraOS->SetBinContent(xbin, ybin, poiPEs);
		    }
		  else //the single crossing ends at the tip of the shower
		    {
		      //the x range in camera coordinates for which the integration will be performed over
		      xrange[0] = xint_a;
		      xrange[1] = image_end_x;
		      //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
		      
		      //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
		      int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
					  ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		      int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
					  ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		      
		      //integrate to get the average number of PEs and then distribute via the Poisson distribution
		      Double_t PEs = gh->Integral(int_range[0], int_range[1]);
		      Double_t poiPEs = r->PoissonD(PEs);
		      
		      //fill both the regular and oversampled histograms
		      cameraOS->SetBinContent(xbin, ybin, poiPEs);
		    }
		}
	      else if(xint_a > xint_b) //now to consder the double crossings of the air shower through a pixel
		{
		  //~ cout << "Double crossing at x1 = " << xint_b << " x2 = " << xint_a << " for bin (" << xbin << ", " << ybin << ")" << endl;
		  
		  //the x range in camera coordinates for which the integration will be performed over
		  xrange[0] = xint_b;
		  xrange[1] = xint_a;
		  //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
		  
		  //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
		  int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
				      ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
				      ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  
		  //integrate to get the average number of PEs and then distribute via the Poisson distribution
		  Double_t PEs = gh->Integral(int_range[0], int_range[1]);
		  Double_t poiPEs = r->PoissonD(PEs);
		  
		  //fill both the regular and oversampled histograms	
		  cameraOS->SetBinContent(xbin, ybin, poiPEs);
		}
	      else if(xint_a < xint_b)
		{
		  //~ cout << "Double crossing at x1 = " << xint_a << " x2 = " << xint_b << " for bin (" << xbin << ", " << ybin << ")" << endl;
		  
		  //the x range in camera coordinates for which the integration will be performed over
		  xrange[0] = xint_a;
		  xrange[1] = xint_b;
		  //~ cout << "Integration range: [" << xrange[0] << ", " << xrange[1] << "]" << endl;
		  
		  //here we translate the camera x range into a domain of the gaisser hillas function to get the number of PEs in a pixel
		  int_range[0] = sqrt((xrange[0] - as_actual_ang[0]) * (xrange[0] - as_actual_ang[0]) + 
				      ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[0] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  int_range[1] = sqrt((xrange[1] - as_actual_ang[0]) * (xrange[1] - as_actual_ang[0]) + 
				      ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1]) * ((xrange[1] * showerslope + showerintercepty) - as_actual_ang[1])) / as_actual_ang[2] * shower_percent * (X_max - X_0) + X_0;
		  
		  //integrate to get the average number of PEs and then distribute via the Poisson distribution
		  Double_t PEs = gh->Integral(int_range[0], int_range[1]);
		  Double_t poiPEs = r->PoissonD(PEs);
		  
		  //fill both the regular and oversampled histograms	
		  cameraOS->SetBinContent(xbin, ybin, poiPEs);
		}
	    }
	}
      
      Double_t latsigma = 0; //the stdev for the lateral gaussian (to be properly implemented)
      
      for(int xbin = 1; xbin <= cameraOS->GetNbinsX(); xbin++) //loop through all the bins of the oversampled histogram
	{
	  for(int ybin = 1; ybin <= cameraOS->GetNbinsY(); ybin++)
	    {
	      if(cameraOS->GetBinContent(xbin, ybin) > 0) //if the bin isn't empty
		{
		  int numPEs = (int) cameraOS->GetBinContent(xbin, ybin); //store the number of PEs in that bin
		  Double_t dist;
		  Double_t latslope = -1. / showerslope; //lateral distribution is normal to the direction of the shower in the camera
		  Double_t angle = atan(latslope); //get the normal angle
		  
		  for(int i = 0; i < numPEs; i++) //loop through the number of PEs within that bin
		    {
		      dist = r->Gaus(0, latsigma); //get a distance away from the pixel according to the gaussian distribution
		      
		      //decompose into x and y components
		      Double_t x = dist * cos(angle) + xbin;
		      Double_t y = dist * sin(angle) + ybin;
		      
		      if((x >= 1 && x <= cameraOS->GetNbinsX()) && (y >= 1 && y <= cameraOS->GetNbinsY())) //if the distibuted distance is within the camera
			{
			  cameraOSLat->Fill(cameraOS->GetXaxis()->GetBinCenter(x), cameraOS->GetYaxis()->GetBinCenter(y), 1); //fill the bin with a single PE
			}
		    }
		}
	    }
	}
      
      Double_t blursigma = 1 * oversample;
      
      for(int xbin = 1; xbin <= cameraOSLat->GetNbinsX(); xbin++) //loop through the bins of the oversampled and laterally distributed histogram
	{
	  for(int ybin = 1; ybin <= cameraOSLat->GetNbinsY(); ybin++)
	    {
	      if(cameraOSLat->GetBinContent(xbin, ybin) > 0) //if the bin is not empty
		{
		  int numPEs = (int) cameraOSLat->GetBinContent(xbin, ybin); //get the number of PEs in that bin
		  Double_t x, y;
		  
		  for(int i = 0; i < numPEs; i++) //loop through the number of PEs
		    {
		      x = r->Gaus(0, blursigma) + xbin; //generate a random x distance away from the bin according to a gaussian
		      y = r->Gaus(0, blursigma) + ybin; //generate a random y distance away from the bin according to a gaussian
		      
		      if((x >= 1 && x <= cameraOSLat->GetNbinsX()) && (y >= 1 && y <= cameraOSLat->GetNbinsY()))
			{
			  cameraBlur->Fill(cameraOSLat->GetXaxis()->GetBinCenter(x), cameraOSLat->GetYaxis()->GetBinCenter(y), 1); //fill the bin with a single PE
			}
		    }
		}
	    }
	}
      
      totalCapturedPEs = 0;
      
      for(int i = 1; i <= cameraBlur->GetNbinsX(); i++)
	{
	  for(int j = 1; j <= cameraBlur->GetNbinsY(); j++)
	    {
	      if(cameraBlur->GetBinContent(i, j) > 0 && cameraBlur->GetBinContent(i, j) < triggerthresh)
		{
		  cameraBlur->SetBinContent(i, j, 0);
		}
	      totalCapturedPEs += cameraBlur->GetBinContent(i, j);
	    }
	}
      
      ASDataTree->Fill();
    }
  
  TTree *ASParamTree = new TTree("ASParams", "Air Shower Trial Params");
  ASParamTree->Branch("tau_azimuth", &azi, "azi/D");
  ASParamTree->Branch("tau_elevation", &elv, "elv/D");
  ASParamTree->Branch("tau_distance", &l, "l/D");
  ASParamTree->Branch("tau_energy", &E, "E/D");
  ASParamTree->Branch("shower_length", &s, "s/D");
  ASParamTree->Branch("detector_alt", &det_alt, "det_alt/I");
  ASParamTree->Branch("total_avg_generated_PEs", &totalPEs, "totalPEs/D");
  ASParamTree->Branch("n_trials", &trials, "trials/I");
  ASParamTree->Branch("n_success", &success, "success/I");
  ASParamTree->Branch("n_full", &full, "full/I");
  ASParamTree->Fill();
  
  printf("Number of visible trials: %d\n", success);
  printf("Number of fully contained trials: %d\n", full);
  printf("Total number trials: %d\n", trials);
  printf("Visible percentage: %.2f%%\n", (double) success / (double) trials * 100.);
  printf("Fully contained percentage: %.2f%%\n", (double) full / (double) trials * 100.);
  
  //~ ASDataTree->Print();
  f->Write();
  f->Close();
  delete(f);
}

void ReadASTrialSims(TString filetitle)
{
  TFile *f = new TFile(filetitle, "READ");
  
  if(f->IsZombie())
    {
      cout << "Error reading root file: " << filetitle << endl;
      f->Close();
      exit(-1);
    }
  
  TTree *ASDataTree = (TTree*) f->Get("ASData");
  TTree *ASParamTree = (TTree*) f->Get("ASParams");
  
  Double_t azi, elv, l, E, s, totalPEs;
  int det_alt, n_trials, n_success, n_full;
  
  ASParamTree->SetBranchAddress("tau_azimuth", &azi);
  ASParamTree->SetBranchAddress("tau_elevation", &elv);
  ASParamTree->SetBranchAddress("tau_distance", &l);
  ASParamTree->SetBranchAddress("tau_energy", &E);
  ASParamTree->SetBranchAddress("shower_length", &s);
  ASParamTree->SetBranchAddress("detector_alt", &det_alt);
  ASParamTree->SetBranchAddress("total_avg_generated_PEs", &totalPEs);
  ASParamTree->SetBranchAddress("n_trials", &n_trials);
  ASParamTree->SetBranchAddress("n_success", &n_success);
  ASParamTree->SetBranchAddress("n_full", &n_full);
  
  ASParamTree->GetEntry(0);
  
  Double_t decay_prob, decay_length, totalCapturedPEs;
  int trial_number;
  bool is_visible, fully_contained, tip_trunc, base_trunc;
  
  ASDataTree->SetBranchAddress("trial_number", &trial_number);
  ASDataTree->SetBranchAddress("decay_prob", &decay_prob);
  ASDataTree->SetBranchAddress("is_visible", &is_visible);
  ASDataTree->SetBranchAddress("fully_contained", &fully_contained);
  ASDataTree->SetBranchAddress("tip_truncated", &tip_trunc);
  ASDataTree->SetBranchAddress("base_truncated", &base_trunc);
  ASDataTree->SetBranchAddress("decay_length", &decay_length);
  ASDataTree->SetBranchAddress("total_captured_PEs", &totalCapturedPEs);
  
  string succ_trials = "Successful trials: ";
  
  for(int i = 0; i < n_trials; i++)
    {
      ASDataTree->GetEntry(i);
      
      printf("Statistics for trial %d:\n", trial_number);
      printf("Decay probability: %.3f%%\n", decay_prob * 100.);
      printf("Decay length: %.3f km\n", decay_length);
      
      if(is_visible)
	{
	  if(fully_contained)
	    {
	      printf("Air shower fully contained within detector FoV.\n");
	    }
	  else if(!base_trunc && tip_trunc)
	    {
	      printf("Air shower truncated by upper FoV.\n");
	    }
	  else if(base_trunc && !tip_trunc)
	    {
	      printf("Air shower truncated by lower FoV or horizon.\n");
	    }
	  else
	    {
	      printf("Air shower truncated by upper FoV and lower FoV or horizon.\n");
	    }
	  
	  succ_trials.append(to_string(trial_number) + ", ");
	}
      else
	{
	  printf("Air shower does not develop within detector FoV.\n");
	}
      
      printf("Total captured PEs: %d\n\n", (int) totalCapturedPEs);
    }
  
  printf("Global simulation parameters:\n");
  printf("Shower azimuth: %.1f deg\n", azi);
  printf("Shower elevation: %.1f deg\n", elv);
  printf("Shower distace: %.1f km\n", l);
  printf("Tau energy: %.0e GeV\n", E);
  printf("Shower length: %.3f km\n", s);
  printf("Detector altitude: %d km\n", det_alt);
  printf("Total average generated PEs: %d\n", (int) totalPEs);
  printf("Number of trials: %d\n", n_trials);
  printf("Number of visible trials: %d\n", n_success);
  printf("Number of fully contained trials: %d\n", n_full);
  printf("Visible percentage: %.2f%%\n", (double) n_success / (double) n_trials * 100.);
  printf("Fully contained percentage: %.2f%%\n\n", (double) n_full / (double) n_trials * 100.);
  
  succ_trials.pop_back();
  succ_trials.pop_back();
  cout << succ_trials << endl << endl;
  
  int user_event_number;
  TString input;
  
  TCanvas *c = new TCanvas("canvas", "canvas", 1200, 1200);
  c->cd(1);
  c->SetRightMargin(0.12);
  
  TH2F *cameraBlur = NULL;
  ASDataTree->SetBranchAddress("camera_hist", &cameraBlur);
  ASDataTree->GetEntry(0);
  
  TH2F *cameraBlurBigger = new TH2F("bigger", "bigger", cameraBlur->GetNbinsX(), 0, cameraBlur->GetNbinsX(), cameraBlur->GetNbinsX(), 0, cameraBlur->GetNbinsX());
  cameraBlurBigger->Draw("COLZ");
  
  while(1)
    {
      input = Getline("Enter the event number you would like to view, type 'q' to exit: ");
      
      if(input.EqualTo("q\n") || input.EqualTo(".q\n"))
	{
	  break;
	}
      
      user_event_number = input.Atoi();
      
      if(user_event_number < 1 || user_event_number > n_trials)
	{
	  printf("Invalid input.\n");
	  continue;
	}
      
      ASDataTree->GetEntry(user_event_number - 1);
      
      if(!is_visible)
	{
	  printf("Event number %d does not produce an observable shower.\n", user_event_number);
	  continue;
	}
      
      cameraBlurBigger->Reset();
      
      printf("Statistics for trial %d:\n", user_event_number);
      printf("Decay probability: %.3f%%\n", decay_prob * 100.);
      printf("Decay length: %.3f km\n", decay_length);
      
      if(is_visible)
	{
	  if(fully_contained)
	    {
	      printf("Air shower fully contained within detector FoV.\n");
	    }
	  else if(!base_trunc && tip_trunc)
	    {
	      printf("Air shower truncated by upper FoV.\n");
	    }
	  else if(base_trunc && !tip_trunc)
	    {
	      printf("Air shower truncated by lower FoV or horizon.\n");
	    }
	  else
	    {
	      printf("Air shower truncated by upper FoV and lower FoV or horizon.\n");
	    }
	}
      printf("Total captured PEs: %d\n\n", (int) totalCapturedPEs);
      
      TString camtitle;
      camtitle.Form("Projected air shower in camera for trial %d (%.0f deg azimuth, %.0f deg elevation, %.0f km distance) with %.2f%% decay probability, %.1E GeV tau energy", user_event_number, azi, elv, l, decay_prob * 100., E);
      
      cameraBlurBigger->GetXaxis()->SetTitle("Pixel x");
      cameraBlurBigger->GetYaxis()->SetTitle("Pixel y");
      cameraBlurBigger->GetZaxis()->SetTitle("PE count");
      cameraBlurBigger->SetContour(99);
      cameraBlurBigger->SetStats(0);
      cameraBlurBigger->SetTitle(camtitle);
      
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
	      cameraBlurBigger->SetBinContent(i, j, cameraBlur->GetBinContent(i, j));
	    }
	}
      
      int range = maxbinx - minbinx;
      
      cameraBlurBigger->GetXaxis()->SetRangeUser(minbinx, maxbinx);
      cameraBlurBigger->GetYaxis()->SetRangeUser((maxbiny + minbiny) / 2. - range / 2., (maxbiny + minbiny) / 2. + range / 2.);
      
      c->SetTitle(camtitle);
      c->Modified();
      c->Update();
      
      cout << succ_trials << endl << endl;
    }
  
  f->Close();
  delete(f);
  delete(c);
}

int main()
{
  return 0;
}
