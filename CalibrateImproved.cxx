//To compile:
//Note: GRSISort looks for .cxx extensions when compiling (for example it looks in the myAnalysis directory)
//Alternatively you may use the following to compile:
//g++ myanalysis.cxx -o myanalysis -std=c++0x -I$GRSISYS/GRSISort/include/ `grsi-config --cflags --all-libs --root`

#include <fstream>    
#include <sstream> 
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <functional>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>
#include <math.h>       /* round, floor, ceil, trunc */
#include <time.h>
using namespace std;

#include "TF1.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TPPG.h"
#include "TScaler.h"
#include "TApplication.h"
#include "TLeaf.h"
#include "TStyle.h"
#include "TGainMatch.h"
#include "TSpectrum.h"
#include "TList.h"
#include "TPolyMarker.h"

#ifndef __CINT__ 
#include "TGriffin.h"
#include "TSceptar.h"
#endif

void
Calibrate(const char*& path, const char*& pref, const char*& kname, const char*& lname, const char*& calfile);
void
WriteCombo(const char*& path, const char*& pref, const char*& kname, const char*& lname, const char*& calfile, const char*& calfileb);

int main() {
	//CAL FILE
	const char* calfile = "./calImproved.cal";
	//DIAGNOSTIC TXT FILES
	const char* kname = "./fpeaks.txt";
	const char* lname = "./fpeaks_matched.txt";
	//DATA FILES
	const char* path = "./Data/";
	const char* pref = "fragment10591";

	Calibrate(path, pref, kname, lname, calfile);					//Finds peaks and gainmatches to known energies
	//WriteCombo(path, pref, kname, lname, calfile, calfileb);			//Writes several source calibrations to a single .cal file
	printf(DGREEN "Done." RESET_COLOR "\n");
  	return EXIT_SUCCESS;
}
void
Calibrate(const char*& path, const char*& pref, const char*& kname, const char*& lname, const char*& calfile){

	int nfile=1;
	int ebins=2e3, emax=2e3;
	double thresh=60;
	int stop=1;										//define a subset of entries to analyse (100==1%)
	FILE *ascii = fopen(kname,"w");
	FILE *asciib = fopen(lname,"w");

	//output file
	TFile *of = new TFile(Form("%sTSpectrum_%s.root",path,pref),"recreate");
	//histograms
	TH2F *ensum = new TH2F("grif_energy","GRIFFIN energy sum spectrum; Channel; Energy [keV]",90,0,90,ebins,0,emax);
	TH1F *energy[64];
	for(int i=0;i<64;i++){
		energy[i]= new TH1F(Form("grif_energy_chan%i",i),Form("GRIFFIN energy channel%i; Channel; Energy [keV]",i),ebins,0,emax);
	}
	//

	for(int k=0;k<nfile;k++){
	TFile *rf = new TFile(Form("%s%s_%03i.root",path,pref,k),"read");
	TTree *cypress=(TTree*)rf->Get("FragmentTree");
    	TFragment* frag = 0;
    	cypress->SetBranchAddress("TFragment", &frag);
   	TChannel::ReadCalFromTree(cypress);

	//nullify the current gains BEFORE measuring peak positions
	TChannel* channel = NULL;
	for(int i=0;i<64;i++){
   		channel = TChannel::GetChannelByNumber(i);
		channel->DestroyENGCal();
	}

	//defs
	int nentries = cypress->GetEntries(), grif;
	double frac, block = ((double)nentries/(double)1e5);
		
	for (int i=0;i<nentries/stop;i++) {
		cypress->GetEntry(i);
		grif = frag->GetChannelNumber();
			
			for(int j=0;j<64;j++){
				if(grif==j){
					ensum->Fill(j,frag->GetEnergy());
					if(frag->GetEnergy()>thresh){
						energy[j]->Fill(frag->GetEnergy());
					}
				}
			}

		frac = 100*((double)i/(double)nentries);
		if(i%int(block)==0 && i<nentries){
			printf(DRED "Reading entry %i / %i [%.1f percent done]" RESET_COLOR "\r",i,nentries,frac);	
		}
	}
	rf->Close();
	} //end file loop
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//Energy data
	//**MAIN REQUIREMENT FOR 'INTELLIGENT' CALIBRATION**
	//"Peak X" must appear in the same fp[i][j] array position in all channels being gainmatched! Look at fpeaks.txt to check this.

	//const int speaks = 5;	//Ba-133
	const int speaks = 8;	//Eu-152
	//Ba-133, uncomment as appropriate
	//double Sourcee[Ba133]={80.999,276.404,302.858,356.014,383.859};
	//double Sourcer[Ba133]={0.225,0.777,0.850,1.000,1.077};				//found peak ratios (ch) in 'matching' crystal, default==0
	//Eu-152
	double Sourcee[speaks]={121.782,244.692,344.275,778.903,964.131,1085.914,1112.116,1408.011};
	double Sourcer[speaks]={0.0865,0.1738,0.2445,0.5532,0.6847,0.7712,0.7898,1.0000};	//found peak ratios (ch) in 'matching' crystal, default==0
	double gain = 0.01;									//gain tolerance
	double threshold = 0.05;								//TSpectrum thresh
	int peakX = 10;
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//peak find using TSpectrum
	int ngam=10;			//Eu+X-rays
	int npeaks=10;			//<--should be some number greater than ngam
	int nchan=64;			//default==64
	double fp[15][64]={0,0};	//rows=ngam, cols=channel number
	double fpm[15][128]={0,0};	//rows=ngam, cols=channel number <--ORDERED, ENERGY-MATCHED data
	double temp, xp;	
	int nfound;
	double *xpeaks;			//ROOT6 seems to want this to be a double (ROOT5 demands float)
	TSpectrum* spec[64];
	double tcounts=100;

	for(int i=0;i<nchan;i++){
		//We have to limit the search to under 10 peaks to stop the Peak buffer warnings! 5% threshold seems to work for 152Eu
		spec[i] = new TSpectrum(npeaks);
		if((energy[i]->Integral(1,4000))<tcounts){
			printf(DYELLOW "Very few counts in channel %i, skipping ..." RESET_COLOR "\n",i);
			continue;
		}
		nfound = spec[i]->Search(energy[i],1.2,"",threshold);	//arg: spec,sigma,opt,threshold
	   	printf(DBLUE "Found %d candidate peaks to fit in channel %i" RESET_COLOR "\n",nfound,i);
		if(nfound>ngam){
			printf(DMAGENTA "Found more than %i peaks in channel %i" RESET_COLOR "\n",ngam,i);	
		} else if(nfound==0){
			printf(DRED "Warning: No peaks found in channel %i" RESET_COLOR "\n",i);	
		}		
		//Loop on all found peaks.
		xpeaks = spec[i]->GetPositionX();
		for(int p=0;p<nfound;p++){
	      		xp = xpeaks[p];
			fp[p][i] = xp;
			if(i==0){
				printf("xpeaks == %.1f\n",xp);
			}
	   	}
		//TSpectrum mixes the order of the peaks - correct here with code stolen from internet
		for(int p=0;p<npeaks;p++){
 			for(int q=0;q<(npeaks-p);q++){
 				if(fp[q][i]>fp[q+1][i]){
 					temp=fp[q][i];
 					fp[q][i]=fp[q+1][i];
 					fp[q+1][i]=temp;
 				}
 			}
 		}
	}
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//Here we assign the 'found' channel numbers to their corresponding energies - EXPERIMENTAL PART OF CODE!!
	for(int i=0;i<15;i++){
		for(int j=0;j<nchan;j++){

		bool match = false;
		double erat = fp[i][j]/fp[peakX][j];

			for(int k=0;k<speaks;k++){
				if(erat<=(Sourcer[k]+(Sourcer[k]*gain)) && erat>=(Sourcer[k]-(Sourcer[k]*gain))){
					fpm[i][(j*2)] = fp[i][j]; 
					fpm[i][(j*2)+1] = Sourcee[k];
					match=true;
					break;			
				}
			}
			if(match==false){
				fpm[i][(j*2)] = 0.;
				fpm[i][(j*2)+1] = 0.;
			}
		}
	}
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//Calculate slope/offset and write to .cal file
	Float_t slope[64]={0};
	Float_t offset[64]={0};
	int nm=0;
	double tmp1=0.,tmp2=0.,tmp3=0,tmp4=0.,tmp5=0.,tmp6=0.,tmp7=0.;
	double plimit = 1.E-05;
	TChannel* channel = NULL;

	for(int i=0;i<nchan;i++){

	//linear regression
	for(int j=0;j<15;j++){
		if(fpm[j][(i*2)]==0. || fpm[j][(i*2)+1]==0.) continue;
	
		tmp1+=((fpm[j][(i*2)])*(fpm[j][(i*2)+1]));	//sum(x.y)
		tmp2+=(fpm[j][(i*2)]);				//sum(x)
		tmp3+=(fpm[j][(i*2)+1]);			//sum(y)
		tmp4+=pow((fpm[j][(i*2)]),2);			//sum(x^2)

		nm+=1;		
		}
	double xave=tmp2/double(nm);								//x-mean
	double yave=tmp3/double(nm);								//y-mean
	double grad=((double(nm)*tmp1)-(tmp2*tmp3))/((double(nm)*tmp4)-pow(tmp2,2));		//gradient
	double icpt=(yave-(grad*xave));								//intercept

	//Pearson correlation coefficient (linearity of fit)
	for(int j=0;j<15;j++){
		if(fpm[j][(i*2)]==0. || fpm[j][(i*2)+1]==0.) continue;

		tmp5+=((fpm[j][(i*2)]-xave)*(fpm[j][(i*2)+1]-yave));	//sum(x-xave)(y-yave)
		tmp6+=pow((fpm[j][(i*2)]-xave),2.);			//sum(x-xave)^2
		tmp7+=pow((fpm[j][(i*2)+1]-yave),2.);			//sum(y-yave)^2

		}
	double pearson = tmp5/(sqrt(tmp6)*sqrt(tmp7));
	double residual = abs(1.-pearson);

	offset[i]=float(icpt);
	slope[i]=float(grad);

	//return to defaults if a problem occurs or a detector isn't present
	if(std::isnan(slope[i]) || std::isnan(offset[i])){
		offset[i]=0.;
		slope[i]=1.;
	}

	if(residual>plimit){
		printf(DRED "Writing gain and offset to channel %i [m=%.3f, c=%.3f, N=%i, RSQ=%.3f, Res=%.2e]" RESET_COLOR "\n",i,slope[i],offset[i],nm,pearson,residual);		
	} else {
		printf(DGREEN "Writing gain and offset to channel %i [m=%.3f, c=%.3f, N=%i, RSQ=%.3f, Res=%.2e]" RESET_COLOR "\n",i,slope[i],offset[i],nm,pearson,residual);
	}
   	channel = TChannel::GetChannelByNumber(i);
	channel->AddENGCoefficient((Float_t)offset[i]);
	channel->AddENGCoefficient((Float_t)slope[i]);
	TChannel::WriteCalFile(calfile);
	//initialise
	tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=tmp7=nm=0;	
	}
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//write out fp in nice format for Excel
	for(int i=0;i<nchan;i++){
		//if(i>15) continue;
		fprintf(ascii,"chan%i\t",i);
	}
	fprintf(ascii,"\n");
	for(int i=0;i<15;i++){
		for(int j=0;j<nchan;j++){
			fprintf(ascii,"%.1f\t",fp[i][j]);	
		}
		fprintf(ascii,"\n");
	}
	for(int i=0;i<15;i++){
		for(int j=0;j<int(nchan*2);j++){
			fprintf(asciib,"%.1f\t",fpm[i][j]);	
		}
		fprintf(asciib,"\n");
	}
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//write spectra
	of->cd();
	ensum->Write();
	for(int i=0;i<nchan;i++){
		energy[i]->Write();
	}
	ensum->Reset();
	for(int i=0;i<nchan;i++){
		energy[i]->Reset();
	}
	fclose(ascii);
	fclose(asciib);
	of->Close();
  return;
}
void WriteCombo(const char*& path, const char*& pref, const char*& kname, const char*& lname, const char*& calfile, const char*& calfileb){
	//user must provide data in "fpeaks_matched.txt" for each source together in a single file
	//TFile input
	TFile *rf = new TFile(Form("%s%s_%03i.root",path,pref,0),"read");
	TTree *cypress=(TTree*)rf->Get("FragmentTree");
    	TFragment* frag = 0;
    	cypress->SetBranchAddress("TFragment", &frag);
   	TChannel::ReadCalFromTree(cypress);
	TChannel* channel = NULL;

	ifstream ina;
	ina.open("/data1/griffin/S1607/fpeaks_matched_Ba133_Eu152_runs08593_4.txt");
	//ina.open("/data1/griffin/S1607/nothing");

	if (rf->IsZombie()) {
  		printf(DRED "No ROOT file provided. Check path." RESET_COLOR "\n");
  		return;
  	} else {
  		printf(DYELLOW "ROOT file %s opened." RESET_COLOR "\n",pref);
  	}
	if (not (ina.is_open())) {
		printf(DRED "No combined calibration data provided. Check path." RESET_COLOR "\n");
    		return;
	} else if(ina.is_open()){
		printf(DYELLOW "I have a combined calibration file I've been working on ..." RESET_COLOR "\n");
	}

	const int rows = 50;
	const int cols = 32;
	string line[rows][cols];

	size_t i=0;
	int nline=0;
	while(!ina.eof() && i < rows) {

		string tmp; 		
		getline(ina, tmp);
		stringstream ss(tmp); 	
                      			
		size_t j=0; 		
		while (ss) {
  		ss >> line[i][j];
  		++j;
		}	
  	++i;
	nline++;
	}
	//the above works ...
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
	//Calculate slope/offset and write to .cal file
	Float_t slope[64]={0};
	Float_t offset[64]={0};
	int nm=0;
	int nchan = cols/2;
	double tmp1=0.,tmp2=0.,tmp3=0,tmp4=0.,tmp5=0.,tmp6=0.,tmp7=0.;
	double plimit = 1.E-05;

	for(int k=0;k<nchan;k++){
	//linear regression
	for(int l=0;l<rows;l++){

    		std::stringstream ss,tt;	//problem here
    		ss<<line[l][(k*2)];
		tt<<line[l][(k*2)+1];
	
    		double x=0.0;
		double y=0.0;
    		ss>>x;
		tt>>y;

		if(x==0. || y==0.) continue;
	
		tmp1+=x*y;	
		tmp2+=x;	
		tmp3+=y;
		tmp4+=pow(x,2);

		nm+=1;
		ss.str("");
		tt.str("");	
		}

	double xave=tmp2/double(nm);								//x-mean
	double yave=tmp3/double(nm);								//y-mean
	double grad=((double(nm)*tmp1)-(tmp2*tmp3))/((double(nm)*tmp4)-pow(tmp2,2));		//gradient
	double icpt=(yave-(grad*xave));								//intercept
	offset[k]=float(icpt);
	slope[k]=float(grad);

	//Pearson correlation coefficient (linearity of fit)
	for(int l=0;l<rows;l++){
		
		std::stringstream ss,tt;
    		ss<<line[l][(k*2)];
		tt<<line[l][(k*2)+1];
	
    		double x=0.0;
		double y=0.0;
    		ss>>x;
		tt>>y;

		if(x==0. || y==0.) continue;

		tmp5+=((x-xave)*(y-yave));
		tmp6+=pow((x-xave),2.);
		tmp7+=pow((y-yave),2.);

		ss.str("");
		tt.str("");	
		}

	double pearson = tmp5/(sqrt(tmp6)*sqrt(tmp7));
	double residual = abs(1.-pearson);
	
	if(residual>plimit){
		printf(DRED "Writing gain and offset to channel %i [m=%.3f, c=%.3f, N=%i, RSQ=%.3f, Res=%.2e]" RESET_COLOR "\n",k,grad,icpt,nm,pearson,residual);		
	} else {
		printf(DGREEN "Writing gain and offset to channel %i [m=%.3f, c=%.3f, N=%i, RSQ=%.3f, Res=%.2e]" RESET_COLOR "\n",k,grad,icpt,nm,pearson,residual);
	}
   	channel = TChannel::GetChannelByNumber(k);
	channel->DestroyENGCal();
	channel->AddENGCoefficient((Float_t)offset[k]);	//something breaking here ... 
	channel->AddENGCoefficient((Float_t)slope[k]);
	TChannel::WriteCalFile(calfileb);
	//initialise
	tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=tmp7=nm=0;	
	}
		printf(DBLUE "Calibration written to %s" RESET_COLOR "\n",calfileb);
	//*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
   	ina.close();
 return;
}
