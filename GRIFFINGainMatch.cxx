///////////////////////////////////////////////////////////////
// Gain matching script across multiple energies
// Connor Natzke Summer 2018
// 
// Compile with 
// g++ myAnalysis.cxx -o myAnalysis -std=c++0x -I$GRSISYS/GRSISort/include/ `grsi-config --cflags --all-libs --root`
///////////////////////////////////////////////////////////////

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
#include "TSystem.h"
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

/////////////////////////////////////
// Definitions
/////////////////////////////////////
void Calibrate(const char*& dataFile, const char*& kName, const char*& lName, const char*& calFile);

/////////////////////////////////////
// Command line argument parsing (stolen from online, credit to iain)
/////////////////////////////////////
class InputParser{
    public:
        InputParser(int &argc, char** argv){
            for(int i=1;i<argc;i++)
                this->tokens.push_back(std::string(argv[i]));
        }

        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr = std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }

        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

/////////////////////////////////////
// Main call to parse command line arguments
/////////////////////////////////////
int main(int argc, char **argv){
    
    InputParser input(argc, argv);

    // Calibration file
    const char* calFile = input.getCmdOption("-c").c_str();
    // Data file
    const char* dataFile = input.getCmdOption("-f").c_str();

    if (input.cmdOptionExists("-h")){
        cout << "A program used to generate calibration files for the GRIFFIN spectrometer. Based on a script initially written by Mike Bowry. Adapted by Connor Natzke\n" << endl;
        cout << "-c <myCal.cal>\t\t Use calibration file myCal" << endl;
        cout << "-f <myFragment.root>\t Use sorted fragment myFragment" << endl;
        cout << "-h\t\t\t Displays this help message" << endl;
        cout << "\nExample usage:" << endl;
        cout << "\t ./GRIFFINGainMatch -c myCal.cal -f myFragment.root" << endl;
    }
    
    // ---------------Checking if requisite files have been passed------------------
    else if (calFile && !calFile[0]){
        printf(DRED "Calibration file and data file required \nUse -h for help" RESET_COLOR "\n");
    }
    else if (!calFile || calFile[0]){
        cout << DGREEN "Found cal file: " RESET_COLOR << calFile << endl;
        cout << DGREEN "Found data file: " RESET_COLOR << dataFile << endl;

        //Diagnostic text files
        const char* kName = "./fPeaks.txt";
        const char* lName = "./fPeaksMatched.txt";
        
        // Creates energy calibration file 
        Calibrate(dataFile, kName, lName, calFile);
        cout << DGREEN "Finished energy calibrations" RESET_COLOR << endl;
        
    }
        
    return EXIT_SUCCESS;
}

/////////////////////////////////////
// Calibrating energy
/////////////////////////////////////
void Calibrate(const char*& dataFile, const char*& kName, const char*& lName, const char*& calFile){

    int nFile = 1;
    int eBins =8e3, eMax = 8e3;
    double thresh = 60;
    int stop = 1;                       // Define subset of entries to analyze (100=1%)
    FILE *ascii = fopen(kName, "w");
    FILE *asciib = fopen(lName, "w");

    // Output file 
    const char* path = gSystem->pwd();
    TFile *of = new TFile(Form("%sTSpectrum.root", path),"recreate");

    // Histograms for output file
    TH2F *gE = new TH2F("gE", "GRIFFIN energy spectrum; Channel; Energy (keV)", 90, 0, 90, eBins, 0, eMax);
    TH1F *energy[64];
    // Loop through each channel
    for (int i=0;i<64;i++){
        energy[i] = new TH1F(Form("gE_%i", i),Form("GRIFFIN energy channel%i;Channel;Energy (keV)",i),eBins, 0, eMax);
    }
// NOTE: use loop when there are more than 1 data file. Also need to add prefix variable, see CalibrateImproved.cxx
//    for (int k=0;k<nfile;k++){
//        TFile *rf = new TFile(Form("%s%s_%03i.root",path, prefix, k),"read");
//        // etc
//    }
        
    TFile *rf = new TFile(dataFile);
    cout << DGREEN "Read in data file" RESET_COLOR << endl;
    TTree *cypress = (TTree*)rf->Get("FragmentTree");
        TFragment *frag = 0;
        cypress->SetBranchAddress("TFragment",&frag);
    TChannel::ReadCalFromTree(cypress);

    // Destroy current calibration before finding peaks
    TChannel *channelN = NULL;
    for (int i=0;i<64;i++){
        channelN = TChannel::GetChannelByNumber(i);
        channelN->DestroyENGCal();
    }

    // Definitions
    int nEntries = cypress->GetEntries(), grif;
    double frac, block = ((double)nEntries/(double)1e5);

    for (int i=0;i<nEntries/stop;i++){
        cypress->GetEntry(i);
        grif = frag->GetChannelNumber();
        
        // Filling histograms
        for (int j=0;j<64;j++){
            if (grif==j){
                gE->Fill(j,frag->GetEnergy());
                if (frag->GetEnergy()>thresh){
                    energy[j]->Fill(frag->GetEnergy());
                }
            }
        }

        frac = 100*((double)i/(double)nEntries);
        if (i%int(block)==0 && i< nEntries){
            printf(DRED "Reading entry %i / %i [%.1f percent done]" RESET_COLOR "\r",i,nEntries,frac);
        }
    }
    rf->Close();
    cout << DRED "\nRead in data, starting analysis" RESET_COLOR << endl;
// } End file loop here

    //------------------ Energy Data --------------------------------
    // 152Eu
    const int sPeaks = 8;
    // Energies
    vector<double> sourceE{121.782,244.692,344.275,778.903,964.131,1085.914,1112.116,1408.011};

    // Energy Peak ratios (ch) in 'matching' crystal, default =0
    vector<double> sourceR(sourceE.size());
    for (int i=0;i<sourceE.size();i++){
        sourceR[i] = sourceE[i]/sourceE.back();
    }
    
    // Gain tolerance
    double gain = 0.01;
    double threshold = 0.05;
    int peakX = 10;

    //------------------ Peak Search with TSpectrum --------------------------------

    int nGam = 10;      // Eu gamma+X-rays
    int nPeaks = 10;    // Should be >= nGam
    int nChan = 64;     // Default = 64
    double fp[15][64] = {0,0};      // rows = nGam, cols = channel
    double fpm[15][128] = {0,0};    // rows = nGam, cols = channel, ordered and calibrated
    double temp, xp;
    int nFound;
    double *xPeaks;
    TSpectrum *spec[64];
    double tCounts = 100;

	for(int i=0;i<nChan;i++){
		//We have to limit the search to under 10 peaks to stop the Peak buffer warnings! 5% threshold seems to work for 152Eu
		spec[i] = new TSpectrum(nPeaks);
		if((energy[i]->Integral(1,eMax))<tCounts){
			printf(DYELLOW "Very few counts in channel %i, skipping ..." RESET_COLOR "\n",i);
			continue;
		}
		nFound = spec[i]->Search(energy[i],1.2,"",threshold);	//arg: spec,sigma,opt,threshold
	   	printf(DBLUE "Found %d candidate peaks to fit in channel %i" RESET_COLOR "\n",nFound,i);
		if(nFound>nGam){
			printf(DMAGENTA "Found more than %i peaks in channel %i" RESET_COLOR "\n",nGam,i);	
		} else if(nFound==0){
			printf(DRED "Warning: No peaks found in channel %i" RESET_COLOR "\n",i);	
		}		
		//Loop on all found peaks.
		xPeaks = spec[i]->GetPositionX();
		for(int p=0;p<nFound;p++){
	      		xp = xPeaks[p];
			fp[p][i] = xp;
			if(i==0){
				printf("xpeaks == %.1f\n",xp);
			}
	   	}
		//TSpectrum mixes the order of the peaks - correct here with code stolen from internet
		for(int p=0;p<nPeaks;p++){
 			for(int q=0;q<(nPeaks-p);q++){
 				if(fp[q][i]>fp[q+1][i]){
 					temp=fp[q][i];
 					fp[q][i]=fp[q+1][i];
 					fp[q+1][i]=temp;
 				}
 			}
 		}
	}

    //------------------- Assigning found peaks corresponding energies -------------------
    
    for (int i=0;i<15;i++){
        for (int j=0;j<nChan;j++){
            bool match = false;
            double eRat = fp[i][j]/fp[peakX][j];

            for (int k=0;k<sPeaks;k++){
                if (eRat<=(sourceR[k]+(sourceR[k]*gain)) && eRat>=(sourceR[k]-(sourceR[k]*gain))){
                    fpm[i][(j*2)] = fp[i][j];
                    fpm[i][(j*2)+1] = sourceE[k];
                    match = true;
                    break;
                }
            }
            if (match==false){
                fpm[i][(j*2)] = 0;
                fpm[i][(j*2)+1] = 0;
            }
        }
    }

    //------------------- Calculating slope/offset to write to file ---------------------------

    Float_t slope[64] = {0};
    Float_t offset[64] = {0};
    int nm = 0;
    double tmp1 = 0., tmp2 = 0., tmp3 = 0., tmp4 = 0., tmp5 = 0., tmp6 = 0., tmp7 = 0.;
    double pLimit = 1.E-05;
    TChannel *channel = NULL;

    for (int i=0;i<nChan;i++){
        // Linear regression
        for (int j=0;j<15;j++){
            if (fpm[j][i*2]==0. || fpm[j][(i*2)+1] == 0.) continue;

            tmp1+=((fpm[j][(i*2)])*(fpm[j][(i*2)+1]));  //sum(x,y)
            tmp2+=(fpm[j][(i*2)]);                      //sum(x)
            tmp3+=(fpm[j][(i*2)+1]);                    //sum(y)
            tmp4+=pow((fpm[j][(i*2)]),2);               //sum(x^2)

            nm+=1;
        }
        double xave = tmp2/double(nm);
        double yave = tmp3/double(nm);
        double grad = ((double(nm)*tmp1)-(tmp2*tmp3))/((double(nm)*tmp4)-pow(tmp2,2));
        double icpt = (yave-(grad*xave));

        // Pearson correlation coeff (linearity of fit)
        for (int j=0;j<15;j++){
            if (fpm[j][i*2]==0. || fpm[j][(i*2)+1] == 0.) continue;

            tmp5+=((fpm[j][(i*2)]-xave)*(fpm[j][(i*2)+1]-yave)); //sum(x-xave,y-yave)
            tmp6+=pow((fpm[j][(i*2)]-xave),2.); //sum(x-xave)^2
            tmp7+=pow((fpm[j][(i*2)+1]-yave),2.); //sum(y-yave)^2
        }
        double pearson = tmp5/(sqrt(tmp6)*sqrt(tmp7));
        double residual = abs(1.-pearson);

        offset[i] = float(icpt);
        slope[i] = float(grad);

        // Defaults if detector is not present or a problem occurs
        if (std::isnan(slope[i]) || std::isnan(offset[i])){
            offset[i] = 0.;
            slope[i] = 1.;
        }

	    if(residual>pLimit){
		    printf(DRED "Writing gain and offset to channel %i [m=%.3f, c=%.3f, N=%i, RSQ=%.3f, Res=%.2e]" RESET_COLOR "\n",i,slope[i],offset[i],nm,pearson,residual);		
	    } else {
		    printf(DGREEN "Writing gain and offset to channel %i [m=%.3f, c=%.3f, N=%i, RSQ=%.3f, Res=%.2e]" RESET_COLOR "\n",i,slope[i],offset[i],nm,pearson,residual);
	    }      
   	    
        channel = TChannel::GetChannelByNumber(i);
	    channel->AddENGCoefficient((Float_t)offset[i]);
	    channel->AddENGCoefficient((Float_t)slope[i]);
	    TChannel::WriteCalFile(calFile);
	    
        //initialise
	    tmp1=tmp2=tmp3=tmp4=tmp5=tmp6=tmp7=nm=0;	
    }

    //------------------- Write out peak matrix in format for Excel ---------------------

    for (int i=0;i<nChan;i++){
        fprintf(ascii, "chan%i\t",i);
    }
    fprintf(ascii,"\n");
    for (int i=0;i<15;i++){
        for (int j=0;j<nChan;j++){
            fprintf(ascii, "%.1f\t",fp[i][j]);
        }
        fprintf(ascii,"\n");
    }
    for (int i=0;i<15;i++){
        for (int j=0;j<int(nChan*2);j++){
            fprintf(asciib,"%.1f\t",fpm[i][j]);
        }
        fprintf(asciib, "\n");
    }
    
    //------------------- Write spectra ---------------------

    of->cd();
    gE->Write();
    for (int i=0;i<nChan;i++){
        energy[i]->Write();
    }
    gE->Reset();
    for (int i=0;i<nChan;i++){
        energy[i]->Reset();
    }
    fclose(ascii);
    fclose(asciib);
    of->Close();
    
    return;  
}

