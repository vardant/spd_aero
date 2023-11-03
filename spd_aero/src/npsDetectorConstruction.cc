#include "npsDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4MultiUnion.hh"
#include <iomanip>

#include "G4NistManager.hh"

#include <fstream>

using namespace std;

// Single module of the SPD's aerogel detector.

npsDetectorConstruction::npsDetectorConstruction() {

  //Read in reflector parameters.

  ifstream fin;
  fin.open("reflector.inp");

  istringstream iss;
  G4String line;
  
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> air_gap;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> gapRefInd;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> refFlag;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> refName;
  getline(fin,line); iss.clear(); iss.str(line);
  iss >> refNumData;
  refWL = new G4double[refNumData];

  //For the specular reflector, read real and imaginary parts of refractive
  //index. For the diffuse reflector, read reflectivity.
  if (refFlag==1) {
    refReIndex = new G4double[refNumData];
    refImIndex = new G4double[refNumData];
    for (G4int i=refNumData-1; i>-1; i--)
      fin >> refWL[i] >> refReIndex[i] >> refImIndex[i];
  }
  else {
    refRefl = new G4double[refNumData];
    for (G4int i=refNumData-1; i>-1; i--)
      fin >> refWL[i] >> refRefl[i];
  }
  
  getline(fin,line);
    
  //Read refractive index of substrate of reflector.
  getline(fin,line); iss.clear(); iss.str(line);
  //  cout << "line: : " <<  line << endl;
  //  getchar();
  iss >>  subRefrIndex;
  getline(fin,line);  iss.clear(); iss.str(line);
  //  cout << "line: " << line << endl;
  //  getchar();
  iss >> phDetFlag;
  fin.close();

  air_gap *= mm;

  for (G4int i=0; i<refNumData; i++) refWL[i] *= nanometer;

  //Print out parameters of reflector.
    
  G4cout << "npsDetectorConstruction::npsDetectorConstruction: input data:"
	  << G4endl;
  G4cout << "   Air gap = " << air_gap/mm << " mm" << G4endl;
  G4cout << "   Gap ref.index = " << gapRefInd << G4endl;
  G4cout << "   Reflector: " << refName << ", refFlag = " << refFlag << ", ";
  if (refFlag==0)
    G4cout << "diffuse reflector";
  else
    G4cout << "specular reflector";
  G4cout << "." << G4endl;

  G4cout << "   Reflector data:" << G4endl;
  for (G4int i=refNumData-1; i>-1; i--) {
    G4cout << "   " << refWL[i]/nanometer << " ";
    if (refFlag==1)
      G4cout  << refReIndex[i] << " " << refImIndex[i];
    else
      G4cout << refRefl[i];
    G4cout << " " << i << G4endl;
  };

  G4cout << "   Substrate refr. index = " << subRefrIndex;
  if (subRefrIndex == 0.)
    G4cout << ", no substrate between crystal and reflector.";
  else
    G4cout << ", substrate layer between crystal and reflector.";
  G4cout << G4endl;

  G4cout << "   phDetFlag = " << phDetFlag << ", ";
  if (phDetFlag == 0)
    G4cout << "PMT to detect photons" << G4endl;
  else
    G4cout << "SiPM to detect photons" << G4endl;

  //
    
  tedlar_thick = 0.040*mm;   //40um Tedlar
  mylar_thick = 0.025*mm;    // + 25um Mylar
  ////  air_gap = 0.035*mm;        //guess
  glue_thick = 0.035*mm;     //guess

  //Photon tracking test
  //  tedlar_thick = 1*mm;   //40um Tedlar
  //  mylar_thick = 1*mm;    // + 25um Mylar
  //  air_gap = 1*mm;        //guess
  //  glue_thick = 1*mm;     //test

  ///  PMT_diameter = 1.86*cm;
  PMT_diameter = 1.4*cm;
  PMTWin_thick = 1*mm;     //??

  ///  Cathode_diam = 1.5*cm;
  Cathode_diam = 1.2*cm;
  Cathode_X =  6*mm;      //2 6x6 mm^2 SiPM-s
  Cathode_Y = 12*mm;
  ///  Cathode_X = 20*mm;
  ///  Cathode_Y = 20*mm;
  Cathode_thick = 0.1*mm;

  wls_x = 3.*mm;
  wls_y = 14*mm;
  wls_z = 470*mm;

  wls_gap_size = 0.5*mm;
  n_wls_slabs = 5;
  
  block_x = 20.*cm;
  block_y =  8.*cm;
  block_z = 47.*cm;

  mylar_x = block_x + 2*air_gap + 2*mylar_thick;
  mylar_y = block_y + 2*air_gap + 2*mylar_thick;
  mylar_z = block_z + 2*air_gap + 2*mylar_thick;

  tedlar_x = mylar_x + 2*tedlar_thick;
  tedlar_y = mylar_y + 2*tedlar_thick;
  tedlar_z = mylar_z + 2*tedlar_thick;

  counter_x = tedlar_x;
  counter_y = tedlar_y;
  counter_z = tedlar_z + 2*glue_thick +  2*PMTWin_thick;

  expHall_x = counter_x * 1.5;
  expHall_y = counter_y * 1.5;
  expHall_z = counter_z * 1.5;
}

npsDetectorConstruction::~npsDetectorConstruction(){;}

//Construct the module.

G4VPhysicalVolume* npsDetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);

  //  G4Material* Vac    = man->FindOrBuildMaterial("G4_Galactic");

  G4double density;
  G4int nelements;
  G4int ncomponents;
  //  G4double fractionmass;

  G4Element* H  = man->FindOrBuildElement("H");
  G4Element* Si = man->FindOrBuildElement("Si");
  G4Element* C  = man->FindOrBuildElement("C");
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* K  = man->FindOrBuildElement("K");
  G4Element* N  = man->FindOrBuildElement("N");
  G4Element* Cs = man->FindOrBuildElement("Cs");
  G4Element* F  = man->FindOrBuildElement("F");

    //SiO2
  G4Material* SiO2 = new G4Material ("Si Oxyde", density = 2.634*g/cm3,
				     nelements = 2);
  SiO2->AddElement(Si, 1);
  SiO2->AddElement(O, 2);
  
  //H2O
  G4Material* H2O = new G4Material("Water", 1.000*g/cm3, 2);
  H2O->AddElement(H, 2);
  H2O->AddElement(O, 1);

  //Polymethylmethacrylate (PMMA), doped with BBQ. Neglect BBQ content for now.
  G4Material* WLS = new G4Material("WLS", 1.18*g/cm3, 3);
  WLS->AddElement(C, 5);
  WLS->AddElement(H, 8);
  WLS->AddElement(O, 2);

  //Refractive index for wl>400nm: https://refractiveindex.info/
  //?shelf=organic&book=poly(methyl_methacrylate)&page=Szczurowski

  const G4double hc = 1.239841857E-6*m*eV;   //(PDG)

  //---Aerogel optics

  /*
  const G4int nEntries = 20;
  G4double RIndex_Aerogel[nEntries];
  G4double Absorption[nEntries];
  G4double PhotonEnergy[nEntries];
  
  const G4int nEntries_scat = 25;
  G4double Scattering[nEntries_scat];
  G4double PhotonEnergy_scat[nEntries_scat];
  
  for(int i = 0; i < nEntries; i++){
  RIndex_Aerogel[i] = 1.02;   //SPD Tech. Design Report
  }
  */
    
  //Absorption length ---

  vector<G4double> Absorption;
  vector<G4double> PhotonEnergy;

  ifstream Indata_abs;
  Indata_abs.open("./aero_inputs/buzykaev/SAN96-3.7-absl.dat");
  if(!Indata_abs) { // file couldn't be opened
    G4cerr << "Error: Cannot read aerogel attenuation parameters" << G4endl;
    exit(1);
  }

  G4double wl, leng;

  while (Indata_abs >> wl >> leng) {
    wl *= nanometer;
    leng *= cm;
    Absorption.push_back(leng);
    PhotonEnergy.push_back(hc/wl);
  }

  reverse(PhotonEnergy.begin(), PhotonEnergy.end());
  reverse(Absorption.begin(), Absorption.end());

  Indata_abs.close();

  cout << "Aerogel absorption:" << endl;
  for(int i = 0; i < PhotonEnergy.size(); i++){
    G4cout << "wl = " << hc/PhotonEnergy[i]/nm << "nm ";
    G4cout << "PhEn = " << PhotonEnergy[i]/eV << "eV  ";
    G4cout << "Abs = " << Absorption[i]/cm << "cm " << endl;
  }
  //  getchar();
  
  //Scattering length ---

  vector<G4double> Scattering;
  vector<G4double> PhotonEnergy_scat;

  ifstream Indata_scat;
  Indata_scat.open("./aero_inputs/buzykaev/SAN96-scatl.dat");
  if(!Indata_scat) { // file couldn't be opened
    G4cerr << "Error: Cannot read aerogel scattering parameters" << G4endl;
    exit(1);
  }
  
  while (Indata_scat >> wl >> leng) {
    wl *= nanometer;
    leng *= cm;
    Scattering.push_back(leng);
    PhotonEnergy_scat.push_back(hc/wl);
  }

  reverse(PhotonEnergy_scat.begin(), PhotonEnergy_scat.end());
  reverse(Scattering.begin(), Scattering.end());

  Indata_scat.close();

  cout << "Aerogel scattering:" << endl;
  for(int i = 0; i < Scattering.size(); i++){
    G4cout << "PhEn = " << PhotonEnergy_scat[i]/eV << "eV  ";
    G4cout << "Scat = " << Scattering[i]/cm <<"cm  " << G4endl;
  }
  //  getchar();

  // Refractive index -----
  
  vector<G4double> RIndex_Aerogel;
  for (int i=0; i<Absorption.size(); i++)
    RIndex_Aerogel.push_back(1.02);

  const G4int nEntries = Absorption.size();
  const G4int nEntries_scat = Scattering.size();

  G4MaterialPropertiesTable* MPT_Aerogel = new G4MaterialPropertiesTable();
  MPT_Aerogel->AddProperty("RINDEX", PhotonEnergy, RIndex_Aerogel);
  MPT_Aerogel->AddProperty("ABSLENGTH", PhotonEnergy, Absorption);
  MPT_Aerogel->AddProperty("RAYLEIGH", PhotonEnergy_scat, Scattering);

  //Aerogel material
  
  G4Material* Aerog = new G4Material("Aerogel", 0.0922*g/cm3, 3);   //for n=1.02
  Aerog->AddMaterial(SiO2, 62.5*perCent);
  Aerog->AddMaterial(H2O , 37.4*perCent);
  Aerog->AddElement (C , 0.1*perCent);

  Aerog -> SetMaterialPropertiesTable(MPT_Aerogel);

  G4cout << "===== MPT_Aerogel: ============================" << G4endl;
  MPT_Aerogel->DumpTable();
  //  getchar();
  
  //Wavelength shifter, PMMA base material ---

  //Refractive index of PMMA for wl>400nm: https://refractiveindex.info/
  //?shelf=organic&book=poly(methyl_methacrylate)&page=Szczurowski
  //Use dispersive relation.
  /*
  const double wls_ri_wl_min = 0.200;
  const double wls_ri_wl_max = 0.600;
  const double wls_ri_step =   0.010;
  const int n_wls_ri = 41;
  G4double PhotonEnergy_WLS[n_wls_ri];
  G4double RIndex_WLS[n_wls_ri];
  for (G4int k = n_wls_ri - 1; k > -1; k-- ) {
    double wl = wls_ri_wl_max - k*wls_ri_step;
    double n2 = 1. + 0.99654*wl*wl/(wl*wl-0.00787) +
      0.18964*wl*wl/(wl*wl-0.02191) + 0.00411*wl*wl/(wl*wl-3.85727);
    RIndex_WLS[k] = sqrt(n2);
    PhotonEnergy_WLS[k] = hc/(wl*micrometer);
  }
  */

  //PMMA refractive index from Roychowdhury et al. at Brgham Young University

  vector<G4double> PhotonEnergy_WLS;
  vector<G4double> RIndex_WLS;

  ifstream Indata_ri;
  Indata_ri.open("./wls_inputs/buzykaev/ri_pmma.dat");
  if(!Indata_ri) { // file couldn't be opened
    G4cerr << "Error: Cannot read WLS refractive indexes" << G4endl;
    exit(1);
  }

  G4double ri;
  
  while (Indata_ri >> wl >> ri) {
    wl *= nanometer;
    RIndex_WLS.push_back(ri);
    PhotonEnergy_WLS.push_back(hc/wl);
  }

  reverse(PhotonEnergy_WLS.begin(), PhotonEnergy_WLS.end());
  reverse(RIndex_WLS.begin(), RIndex_WLS.end());

  Indata_ri.close();

  const G4int n_wls_ri = RIndex_WLS.size();

  cout << "WLS ref. index:" << endl;
  for(int i = 0; i < n_wls_ri; i++) {
    double wl = hc/PhotonEnergy_WLS[i];
    cout << PhotonEnergy_WLS[i]/eV << " ev  " << wl/nm << " nm  "
  	 << RIndex_WLS[i] << endl;
  }
  //  getchar();

  //WLS-BBQ absorption length.

  vector<double> wl_absl_bbq;
  vector<double> absl_bbq;
  vector<double> PhotonEnergy_absl_bbq;

  ifstream Indata_absl_bbq;
  Indata_absl_bbq.open("./wls_inputs/buzykaev/absl-bbq.dat");
  if(!Indata_absl_bbq) { // file couldn't be opened
    G4cerr << "Error: Cannot read WLS-BBQ absorption parameters" << G4endl;
    exit(1);
  }

  while (Indata_absl_bbq >> wl >> leng) {
    wl *= nanometer;
    leng *= cm;
    absl_bbq.push_back(leng);
    PhotonEnergy_absl_bbq.push_back(hc/wl);
  }

  reverse(PhotonEnergy_absl_bbq.begin(), PhotonEnergy_absl_bbq.end());
  reverse(absl_bbq.begin(), absl_bbq.end());

  Indata_absl_bbq.close();

  const int n_absl_bbq = absl_bbq.size();
  
  cout << "WLS-BBQ absorption:" << endl;
  for(int i = 0; i < n_absl_bbq; i++){
    G4cout << "PhEn = " << PhotonEnergy_absl_bbq[i]/eV << "eV  ";
    G4cout << "Absl = " << absl_bbq[i]/cm << "cm" << endl;
  }
  //  getchar();

  //PMMA absorption length.

  vector<double> wl_absl_pmma;
  vector<double> absl_pmma;
  vector<double> PhotonEnergy_absl_pmma;

  ifstream Indata_absl_pmma;
  Indata_absl_pmma.open("./wls_inputs/buzykaev/absl-pmma.dat");
  if(!Indata_absl_pmma) { // file couldn't be opened
    G4cerr << "Error: Cannot read PMMA absorption parameters" << G4endl;
    exit(1);
  }

  while (Indata_absl_pmma >> wl >> leng) {
    wl *= nanometer;
    leng *= cm;
    absl_pmma.push_back(leng);
    PhotonEnergy_absl_pmma.push_back(hc/wl);
  }

  reverse(PhotonEnergy_absl_pmma.begin(), PhotonEnergy_absl_pmma.end());
  reverse(absl_pmma.begin(), absl_pmma.end());

  Indata_absl_pmma.close();

  const int n_absl_pmma = absl_pmma.size();
  
  cout << "PMMA absorption:" << endl;
  for(int i = 0; i < n_absl_pmma; i++){
    G4cout << "PhEn = " << PhotonEnergy_absl_pmma[i]/eV << "eV  ";
    G4cout << "Absl = " << absl_pmma[i]/cm << "cm" << endl;
  }
  //  getchar();

  //WLS absorption length, estimated as absl_bbq*absl_pmma/(absl_bbq+absl_pmma).

  vector<double> wl_absl_wls;
  vector<double> absl_wls;
  vector<double> PhotonEnergy_absl_wls;

  ifstream Indata_absl_wls;
  Indata_absl_wls.open("./wls_inputs/buzykaev/absl-wls.dat");
  if(!Indata_absl_wls) { // file couldn't be opened
    G4cerr << "Error: Cannot read WLS absorption parameters" << G4endl;
    exit(1);
  }

  while (Indata_absl_wls >> wl >> leng) {
    wl *= nanometer;
    leng *= cm;
    absl_wls.push_back(leng);
    PhotonEnergy_absl_wls.push_back(hc/wl);
  }

  reverse(PhotonEnergy_absl_wls.begin(), PhotonEnergy_absl_wls.end());
  reverse(absl_wls.begin(), absl_wls.end());

  Indata_absl_wls.close();

  const int n_absl_wls = absl_wls.size();

  cout << "WLS absorption:" << endl;
  for(int i = 0; i < n_absl_wls; i++){
    G4cout << "PhEn = " << PhotonEnergy_absl_wls[i]/eV << "eV  ";
    G4cout << "Absl = " << absl_wls[i]/cm << "cm" << endl;
  }
  //  getchar();

  //WLS emission.

  vector<double> wl_emission;
  vector<double> emission;
  vector<double> PhotonEnergy_emission;

  ifstream Indata_emission;
  Indata_emission.open("./wls_inputs/buzykaev/emission-bbq.dat");
  if(!Indata_emission) { // file couldn't be opened
    G4cerr << "Error: Cannot read WLS emission parameters" << G4endl;
    exit(1);
  }

  while (Indata_emission >> wl >> leng) {
    wl *= nanometer;
    leng *= cm;
    emission.push_back(leng);
    PhotonEnergy_emission.push_back(hc/wl);
  }

  reverse(PhotonEnergy_emission.begin(), PhotonEnergy_emission.end());
  reverse(emission.begin(), emission.end());

  Indata_emission.close();

  const int n_emission = emission.size();

  cout << "WLS emission:" << endl;
  for(int i = 0; i < n_emission; i++){
    G4cout << "PhEn = " << PhotonEnergy_emission[i]/eV << "eV  ";
    G4cout << "Emission = " << emission[i] << endl;
  }
  //  getchar();

  G4MaterialPropertiesTable* MPT_WLS = new G4MaterialPropertiesTable();
  MPT_WLS->AddProperty("RINDEX", PhotonEnergy_WLS, RIndex_WLS);
  ///  MPT_WLS->AddProperty("ABSLENGTH", PhotonEnergy_absl_wls, absl_wls,
  ///		       n_absl_wls);
  MPT_WLS->AddProperty("ABSLENGTH", PhotonEnergy_absl_pmma, absl_pmma);
  MPT_WLS->AddProperty("WLSABSLENGTH", PhotonEnergy_absl_bbq, absl_bbq);
  ///MPT_WLS->AddProperty("WLSABSLENGTH",&PhotonEnergy_absl_wls[0],&absl_wls[0],
  ///		       n_absl_wls);
  MPT_WLS->AddProperty("WLSCOMPONENT", PhotonEnergy_emission, emission);
  MPT_WLS->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);

  WLS->SetMaterialPropertiesTable(MPT_WLS);

  G4cout << "===== MPT_WLS: ============================" << G4endl;
  MPT_WLS->DumpTable();
  //  getchar();
  
  // Air
  // 
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  vector<G4double> rindAir;
  for (G4int i=0; i<nEntries; i++) {
    rindAir.push_back(gapRefInd);   //Air @ STP normally
  };
  G4MaterialPropertiesTable *AirMPT = new G4MaterialPropertiesTable();
  AirMPT -> AddProperty("RINDEX", PhotonEnergy, rindAir);
  Air -> SetMaterialPropertiesTable(AirMPT);

  G4cout << "===== AirMPT: ============================" << G4endl;
  AirMPT->DumpTable();
  //  getchar();
  
  // Glass
  //

  density = 2.23*g/cm3;   //Borosilicate glass (wikipedia)
  G4Material* Glass = new G4Material("Glass", density, ncomponents=2);
  Glass->AddElement(Si, 1);
  Glass->AddElement(O,  2);

  vector<G4double> rindGlass;
  for (G4int i=0; i<nEntries; i++) {
    rindGlass.push_back(1.525);              //average of 1.51-1.54
  };

  G4MaterialPropertiesTable *GlassMPT = new G4MaterialPropertiesTable();
  GlassMPT -> AddProperty("RINDEX", PhotonEnergy ,rindGlass);
  Glass -> SetMaterialPropertiesTable(GlassMPT);

  G4cout << "===== GlassMPT: ============================" << G4endl;
  GlassMPT->DumpTable();
  //  getchar();
  
  //Optical glue, or grease ---------
  
  G4Material* OpticalGlue;
  vector<G4double> rindGlue;

  if (phDetFlag == 0) {

    // Optical grease BC630 from Bicron
    //
    density = 1.06*g/cm3;
    OpticalGlue = new G4Material("Silgard", density, ncomponents=1);
    OpticalGlue->AddElement(Si, 1); //not known

    for (G4int i=0; i<nEntries; i++) {
      rindGlue.push_back(1.465);
    };

  }
  else {
  
    // Dow Corning 3145 RTV-CLEAR MIL-A-44146 adhesive seelant
    // (polydimethylsiloxane)
    density = 0.965*g/cm3;
    OpticalGlue = new G4Material("Dow Corning 3145", density, ncomponents=4);
    OpticalGlue->AddElement(Si, 1);
    OpticalGlue->AddElement(O , 1);
    OpticalGlue->AddElement(C , 2);
    OpticalGlue->AddElement(H , 6);
  
    for (G4int i=0; i<nEntries; i++) {
      rindGlue.push_back(1.4);
    };
  }

  G4MaterialPropertiesTable *GlueMPT = new G4MaterialPropertiesTable();
  GlueMPT -> AddProperty("RINDEX", PhotonEnergy, rindGlue);
  OpticalGlue -> SetMaterialPropertiesTable(GlueMPT);

  G4cout << "===== GlueMPT: ============================" << G4endl;
  GlueMPT->DumpTable();
  //  getchar();
  
  // Optical insulation
  density = 1.38;   //from goodfellow
  G4Material* PVF = new G4Material("Polyvinyl fluoride",density,ncomponents=3);
  PVF->AddElement(C, 2);
  PVF->AddElement(H, 3);
  PVF->AddElement(F, 1);

  //Mylar, reflector substrate material.
  //

  G4Material* Mylar = man->FindOrBuildMaterial("G4_MYLAR");

  if (subRefrIndex != 0.) {
    
    //Mylar refractive index.
    vector<G4double> rindMylar;
    for (G4int i=0; i<nEntries; i++) {
      rindMylar.push_back(subRefrIndex);
    };

    G4MaterialPropertiesTable *MylarMPT = new G4MaterialPropertiesTable();
    MylarMPT -> AddProperty("RINDEX", PhotonEnergy, rindMylar);
    Mylar -> SetMaterialPropertiesTable(MylarMPT);

    G4cout << "===== MylarMPT: ============================" << G4endl;
    MylarMPT->DumpTable();
    //    getchar();
  }

  // Bialcali, the photochathode material
  //

  density = 1*g/cm3;   //Does not matter
  G4Material* Bialcali = new G4Material("Bialcali", density, ncomponents=2);
  Bialcali->AddElement(Cs, 1);
  Bialcali->AddElement(K,  1);

//
//	------------- Volumes --------------

  // Aerogel block
  //
  G4Box* block_box = new G4Box("Block_box",block_x/2,block_y/2,block_z/2);

  //Groove to put in WLS slab(s).
  ///  G4Box* groove_box = new G4Box("grove_box",wls_x/2,wls_y/2,wls_z/2);
  G4Box* groove_box = new G4Box("grove_box",wls_x/2,block_y/2,wls_z/2);

  G4RotationMatrix rot;

  G4SubtractionSolid* form_box = new G4SubtractionSolid("Aerogel", block_box,
  groove_box, G4Transform3D(rot, G4ThreeVector(0., 0., 0.)));
  ///groove_box, G4Transform3D(rot, G4ThreeVector(0., block_y/2-wls_y/2, 0.)));

  ///  block_log = new G4LogicalVolume(block_box,Aerog,"Block_log",0,0,0);
  block_log = new G4LogicalVolume(form_box,Aerog,"Block_log",0,0,0);

  // WLS slab
  //
  G4Box* wls_box = new G4Box("WLS_box",wls_x/2,wls_y/2,wls_z/2);

  wls_log = new G4LogicalVolume(wls_box,WLS,"WLS_log",0,0,0);

  // Optical insulation
  //
  G4Box* tedlar_outer =
    new G4Box("Tedlar_solid",tedlar_x/2,tedlar_y/2,tedlar_z/2);

  G4Box* tedlar_inner = new G4Box("Tedlar_cavity",
				 tedlar_x/2-tedlar_thick,
				 tedlar_y/2-tedlar_thick,
				 tedlar_z/2-tedlar_thick);

  G4SubtractionSolid* tedlar_box = new G4SubtractionSolid("Tedlar",
				  tedlar_outer, tedlar_inner);

  // Optical insulation with holes on right and left side.

  G4MultiUnion* tedlar_holes = new G4MultiUnion("tedlar_holes");
  for (int k=0; k<n_wls_slabs; k++) {
    double y = block_y-wls_gap_size-(k+0.5)*(wls_y+wls_gap_size) - block_y/2;
    G4Transform3D trans_tedlar_hole(rot, G4ThreeVector(0., y, 0.));
    G4Tubs*  tedlar_hole = new G4Tubs("tedlar_hole",
				      0., PMT_diameter/2, tedlar_z/2+1,
				      0.*deg, 360.*deg);
    tedlar_holes->AddNode(*tedlar_hole, trans_tedlar_hole);
  }
  tedlar_holes->Voxelize();
  
  G4SubtractionSolid* tedlar_holed = new G4SubtractionSolid("Tedlar_holed",
						    tedlar_box, tedlar_holes);
							    
  tedlar_log = new G4LogicalVolume(tedlar_holed,PVF,"Tedlar",0,0,0);
    
  //Mylar, reflector.
  
  G4Box* mylar_outer = new G4Box("Mylar_solid",mylar_x/2,mylar_y/2,mylar_z/2);

  G4Box* mylar_inner = new G4Box("Mylar_cavity",
				 mylar_x/2-mylar_thick,
				 mylar_y/2-mylar_thick,
				 mylar_z/2-mylar_thick);

  G4SubtractionSolid* mylar_box = new G4SubtractionSolid("Mylar",
				  mylar_outer, mylar_inner);

  G4MultiUnion* mylar_holes = new G4MultiUnion("mylar_holes");
  for (int k=0; k<n_wls_slabs; k++) {
    double y = block_y-wls_gap_size-(k+0.5)*(wls_y+wls_gap_size) - block_y/2;
    G4Transform3D trans_mylar_hole(rot, G4ThreeVector(0., y, 0.));
    G4Tubs*  mylar_hole = new G4Tubs("mylar_hole",
				      0., PMT_diameter/2, mylar_z/2+1,
				      0.*deg, 360.*deg);
    mylar_holes->AddNode(*mylar_hole, trans_mylar_hole);
  }
  mylar_holes->Voxelize();
  
  G4SubtractionSolid* mylar_holed = new G4SubtractionSolid("Mylar",
						   mylar_box, mylar_holes);

  ///  mylar_log=new G4LogicalVolume(mylar_frame,Mylar,"Mylar",0,0,0);
  mylar_log=new G4LogicalVolume(mylar_holed,Mylar,"Mylar",0,0,0);
    
  // PMT Window
  //
  G4Tubs*  PMTWin_tube =
  new G4Tubs("PMTWindow", 0., PMT_diameter/2, PMTWin_thick/2,0.*deg, 360.*deg);

  PMTWin_right_log = new G4LogicalVolume(PMTWin_tube,Glass, "PMTWindow");

  if (phDetFlag == 0) {
    G4Tubs*  Cathode_tube =
      new G4Tubs("Cathode",0., Cathode_diam/2, Cathode_thick/2,0.*deg,360.*deg);
    Cathode_log = new G4LogicalVolume(Cathode_tube, Bialcali, "Cathode");
  }
  else {
    G4Box*  Cathode_box =
      new G4Box("Cathode", Cathode_X/2., Cathode_Y/2, Cathode_thick/2.);
    Cathode_log = new G4LogicalVolume(Cathode_box, Bialcali, "Cathode");
  }

  // Optical glue
  //
  G4Tubs*  glue_tube =
    new G4Tubs("glue", 0., PMT_diameter/2, glue_thick/2, 0.*deg, 360.*deg);

  glue_log = new G4LogicalVolume(glue_tube,OpticalGlue, "Glue");

  // Counter
  //
  G4Box* counter_box = new G4Box("Counter",counter_x/2,counter_y/2,counter_z/2);

  counter_log = new G4LogicalVolume(counter_box,Air,"Counter",0,0,0);

  // The experimental Hall
  //
  G4Box* expHall_box = new G4Box("World",expHall_x/2,expHall_y/2,expHall_z/2);

  expHall_log = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);


  // Place constituents, construct physical volumes.
  //

  G4VPhysicalVolume* expHall_phys =
    new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

  //  G4VPhysicalVolume* counter_phys =
  new G4PVPlacement(0, //no rotation
		    G4ThreeVector(),
		    counter_log, //its logical volume
		    "Counter",   //its name
		    expHall_log,     //its mother  volume
		    false,         //no boolean operation
		    0);  //copy number

  // Insulation for counter
  new G4PVPlacement(0,  //no rotation
  		    G4ThreeVector(),
  		    tedlar_log,      //its logical volume
  		    "Tedlar",          //its name
  		    counter_log,      //its mother  volume
  		    false,              //no boolean operation
  		    0);                 //copy number
  
  //  G4VPhysicalVolume* mylar_phys =
  new G4PVPlacement(0,  //no rotation
  		    G4ThreeVector(),
		    mylar_log,        //its logical volume
		    "Mylar_phys",       //its name
		    counter_log,    //its mother  volume
		    false,            //no boolean operation
		    0);               //copy number
  
  // Glass block for counter

  new G4PVPlacement(0,  //no rotation
		    G4ThreeVector(),
		    block_log,        //its logical volume
		    "Block_phys",     //its name
		    counter_log,        //its mother  volume
		    false,            //no boolean operation
		    0);               //copy number

  for (int k=0; k<n_wls_slabs; k++) {
    
    double y = block_y-wls_gap_size-(k+0.5)*(wls_y+wls_gap_size) - block_y/2;

    new G4PVPlacement(0,  //no rotation
    		      G4ThreeVector(0., y, 0.),
		      wls_log,        //its logical volume
		      "WLS_phys",     //its name
		      counter_log,        //its mother  volume
		      false,            //no boolean operation
		      k+1);               //copy number

    // Glue, window and cathode for counter on both sides.

    G4double x = 0.;
    ///  G4double y = block_y/2-wls_y/2;
    G4double z = block_z/2 + glue_thick/2;
    new G4PVPlacement(0, //no rotation
		      G4ThreeVector(x,y,z),
		      glue_log,    //its logical volume
		      "Glue",      //its name
		      counter_log, //its mother  volume
		      false,         //no boolean operation
		      (k+1)*10);            //copy number
    new G4PVPlacement(0, //no rotation
		      G4ThreeVector(x,y,-z),
		      glue_log,    //its logical volume
		      "Glue",      //its name
		      counter_log, //its mother  volume
		      false,         //no boolean operation
		      (k+1)*10+1);            //copy number

    z = block_z/2 + glue_thick + PMTWin_thick/2;

    new G4PVPlacement(0,    //no rotation
		      G4ThreeVector(x,y,z),
		      PMTWin_right_log, //its logical volume
		      "PMTWindow",      //its name
		      counter_log,    //its mother  volume
		      false,            //no boolean oper.
		      (k+1)*10);               //copy number
    new G4PVPlacement(0,    //no rotation
		      G4ThreeVector(x,y,-z),
		      PMTWin_right_log, //its logical volume
		      "PMTWindow",      //its name
		      counter_log,    //its mother  volume
		      false,            //no boolean oper.
		      (k+1)*10+1);               //copy number

      z = block_z/2 + glue_thick + PMTWin_thick + Cathode_thick/2;
      new G4PVPlacement(0, //no rotation
			G4ThreeVector(x,y,z),
			Cathode_log,  //its logical volume
			"Cathode", //its name
			counter_log, //its mother  volume
			false,       //no boolean operation
			(k+1)*10);          //copy number
      new G4PVPlacement(0, //no rotation
			G4ThreeVector(x,y,-z),
			Cathode_log,  //its logical volume
			"Cathode", //its name
			counter_log, //its mother  volume
			false,       //no boolean operation
			(k+1)*10+1);          //copy number
  }   //k=1,n_wls_slabs

//	------------- Surfaces --------------
//

  G4MaterialPropertiesTable* ReflectorMPT = new G4MaterialPropertiesTable();
  G4OpticalSurface* Reflector = new G4OpticalSurface("Reflector");

  G4double* refKphot;   //Momenta of optical photons in eV units.
  refKphot = new G4double[refNumData];
  for (G4int i=0; i<refNumData; i++) refKphot[i] = hc/refWL[i];

  if (refFlag == 1) {

    ReflectorMPT->AddProperty("REALRINDEX",refKphot,refReIndex,refNumData);
    ReflectorMPT->AddProperty("IMAGINARYRINDEX",refKphot,refImIndex,
			      refNumData);

    Reflector -> SetType(dielectric_metal);
    Reflector -> SetFinish(polished);
    Reflector -> SetModel(glisur);
  }
  else {
    // Non metallic reflector, PTFE (Teflon), VM2000.

    ReflectorMPT -> AddProperty("REFLECTIVITY",refKphot,refRefl,refNumData);
    Reflector -> SetType(dielectric_dielectric);
    Reflector -> SetModel(unified);
    if (refFlag == 0)
      Reflector -> SetFinish(groundfrontpainted); //Purely Lambertian reflection
    else
      Reflector -> SetFinish(polishedfrontpainted); //Purely specular reflection
  }

  Reflector -> SetMaterialPropertiesTable(ReflectorMPT);

  G4cout << "===== ReflectorMPT: ============================" << G4endl;
  ReflectorMPT->DumpTable();
  Reflector->DumpInfo();
  //  getchar();
  
  if (subRefrIndex == 0.) {
    // Reflective front surface of Mylar.
    new G4LogicalSkinSurface("Reflector",mylar_log,Reflector);
  }
  else {
    // Reflective back surface of Mylar.
    // Tedlar borders Mylar from back. Making it reflective, makes effectively
    // Mylar back surface reflective.
    new G4LogicalSkinSurface("Reflector",tedlar_log,Reflector);
    G4cout << "   subRefrIndex = " << subRefrIndex
	   << ", substarate between crystal and reflector" << G4endl;
  }
  
  //Quantum efficiencies.
  //

  int nCat;
  if (phDetFlag == 0)
    nCat = 101;
  else
    nCat = 53;
  
  vector<G4double> wlCat(nCat);
  vector<G4double> kphotCat(nCat);   //Momenta of optical photons in eV units.
  vector<G4double> effCat(nCat);
  vector<G4double> reflCat(nCat);

  if (phDetFlag == 0) {
     wlCat = {675.,670.,665.,660.,655.,650.,645.,640.,635.,630.,
	      625.,620.,615.,610.,605.,600.,595.,590.,585.,580.,
	      575.,570.,565.,560.,555.,550.,545.,540.,535.,530.,
	      525.,520.,515.,510.,505.,500.,495.,490.,485.,480.,
	      475.,470.,465.,460.,455.,450.,445.,440.,435.,430.,
	      425.,420.,415.,410.,405.,400.,395.,390.,385.,380.,
	      375.,370.,365.,360.,355.,350.,345.,340.,335.,330.,
	      325.,320.,315.,310.,305.,300.,295.,290.,285.,280.,
	      275.,270.,265.,260.,255.,250.,245.,240.,235.,230.,
	      225.,220.,215.,210.,205.,200.,195.,190.,185.,180.,
	      175.};
     // Hamamatsu R4125 quantum efficiency (bialcali photocathode, borosilicate
     // window). Taken from the Hamamatsu booklet, p.65.
     effCat = {
       0.0030,0.0035,0.0040,0.0046,0.0052,0.0060,0.0068,0.0077,0.0087,0.0099,
       0.0112,0.0126,0.0141,0.0159,0.0177,0.0198,0.0221,0.0245,0.0272,0.0301,
       0.0332,0.0365,0.0401,0.0440,0.0481,0.0525,0.0572,0.0621,0.0673,0.0728,
       0.0785,0.0846,0.0908,0.0973,0.1041,0.1110,0.1181,0.1255,0.1329,0.1405,
       0.1482,0.1560,0.1638,0.1716,0.1793,0.1870,0.1946,0.2020,0.2092,0.2162,
       0.2229,0.2293,0.2354,0.2411,0.2463,0.2511,0.2554,0.2592,0.2625,0.2651,
       0.2673,0.2688,0.2697,0.2700,0.2688,0.2653,0.2595,0.2517,0.2419,0.2305,
       0.2177,0.2038,0.1891,0.1740,0.1586,0.1434,0.1285,0.1141,0.1004,0.0877,
       0.0758,0.0650,0.0553,0.0466,0.0389,0.0322,0.0264,0.0215,0.0173,0.0138,
       0.0110,0.0086,0.0067,0.0052,0.0040,0.0030,0.0023,0.0017,0.0012,0.0009,
       0.0007};
  }
  else {
     wlCat = {
       883.49,868.0,850.59,835.1,819.62,804.13,773.16,759.61,744.12,730.57,
       717.01,703.46,693.78,682.15,666.66,655.03,645.34,633.72,624.02,614.33,
       604.64,587.2,577.51,567.81,560.06,550.36,540.67,529.04,519.34,511.58,
       503.82,499.95,490.25,482.5,455.4,436.06,414.82,403.24,393.6,385.89,
       378.18,372.4,343.45,333.8,324.16,318.38,312.6,306.84,304.93,301.15,
       299.25,287.92,280.21};

     // Hamamatsu S14160-6050HS MPPC PDE (photon detection efficiency). Data
     // digitized from graph in booklet.
     effCat = {
       4.739,5.572,6.405,7.238,8.237,9.236,11.234,12.233,13.398,14.563,
       15.728,16.893,17.725,19.056,20.72,22.217,23.548,24.879,26.376,27.707,
       29.204,31.533,33.03,34.527,35.858,37.355,38.852,41.014,42.511,44.007,
       45.504,46.169,47.666,49.163,50.828,50.498,48.839,47.51,45.85,44.355,
       42.694,41.198,38.21,36.715,35.221,33.559,32.064,29.738,28.076,22.758,
       21.096,3.646,2.151};
       for (G4int i=0; i<nCat; i++)
	 effCat[i] /= 100.;
  }

  for (G4int i=0; i<101; i++) {
    wlCat[i] *= nanometer;
  };

  for (G4int i=0; i<nCat; i++) kphotCat[i] = hc/wlCat[i];
  
  for (G4int i = 0; i < nCat; i++) {
    reflCat[i] = 0.;
  }

  G4OpticalSurface* surfCat = new G4OpticalSurface("Cathode");

  surfCat -> SetType(dielectric_metal);
  surfCat -> SetFinish(polished);
  surfCat -> SetModel(glisur);

  G4MaterialPropertiesTable* surfCatMPT = new G4MaterialPropertiesTable();
  surfCatMPT -> AddProperty("REFLECTIVITY",kphotCat,reflCat,nCat);
  surfCatMPT -> AddProperty("EFFICIENCY",kphotCat,effCat,nCat);

  surfCat -> SetMaterialPropertiesTable(surfCatMPT);

  new G4LogicalSkinSurface("Cathode",Cathode_log,surfCat);

  G4cout << "===== surfCatMPT: ============================" << G4endl;
  surfCatMPT->DumpTable();
  //  getchar();

  expHall_log->SetVisAttributes (G4VisAttributes::GetInvisible());
  counter_log->SetVisAttributes (G4VisAttributes::GetInvisible());

  // Print the table of materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //  G4double maxStep = 1*meter, maxLength = 1*meter, maxTime = 5*ns,
  //    minEkin = 0*MeV;
  //  block_log->SetUserLimits(new G4UserLimits(maxStep, maxLength, maxTime,
  //					    minEkin));

//
//always return the physical World
//
  return expHall_phys;
}
