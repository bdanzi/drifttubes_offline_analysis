/*
 
 To run it, do:
 
 - Crate a file test.dat via the "Save" button in the browser connected
 to the wds server
 - start ROOT
 root [0] .L read_binary.C+
 root [1] decode("test.dat");
 
 */


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Getline.h"
#include <WaveData.h>
#include <RunHeader.h>
#include <boost/format.hpp>
#include <iostream>

#define NBOARDS 16

typedef struct {
   char tag[3];
   char version;
} FHEADER;

typedef struct {
   char time_header[4];
} THEADER;

typedef struct {
   char bn[2];
   unsigned short board_serial_number;
} BHEADER;

typedef struct {
   char event_header[4];
   unsigned int event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;
} EHEADER;

typedef struct {
   char tc[2];
   unsigned short trigger_cell;
} TCHEADER;

typedef struct {
   char c[1];
   char cn[3];
} CHEADER;

/*-----------------------------------------------------------------------------*/

void decode(const char *filename,int max_Ev,__uint32_t runnum) {


	FHEADER fh;
	THEADER th;
	BHEADER bh;
	EHEADER eh;
	TCHEADER tch;
	CHEADER ch;


	__uint32_t evnum;
	__uint64_t btime;
	__uint8_t ndevs;
	__uint16_t chnum;
	__uint8_t chpdev[2];
	__uint8_t devids[2];
	__uint8_t devtypes[2];


	unsigned int scaler;
	unsigned short voltage[1024];
	double voltage_2[1024];

	double waveform[NBOARDS][18][1024], time[NBOARDS][18][1024];
	float adc_waveform[NBOARDS][16][2048];
	unsigned char tdc_waveform[NBOARDS][16][512];
	unsigned long trg_data[NBOARDS][512];
	float bin_width[NBOARDS][18][1024];
	int i, j, b, chn, n, chn_index, n_boards;
	double t1, t2, dt;
	double amplitude[18];
	char rootfile[256];
	// open the binary waveform file
	FILE *f = fopen(Form("%s", filename), "r");
	if (f == NULL) {
		printf("Cannot find file \'%s\'\n", filename);
		return;
	}
	//open the root file
	strcpy(rootfile, Form("%s", filename));
	if (strchr(rootfile, '.'))
		*strchr(rootfile, '.') = 0;
	strcat(rootfile, ".root");
	TFile *outfile = new TFile(// rootfile
				   (boost::format("run_%1%.root") % runnum).str().c_str()
				   , "RECREATE");

	// define the rec tree



	// create canvas
	//TCanvas *c1 = new TCanvas();
	//c1->Divide(4,4);
	// create graph
	TGraph *g[16];
	TTree *rec;
	// read file header
	fread(&fh, sizeof(fh), 1, f);
	if (fh.tag[0] != 'D' || fh.tag[1] != 'R' || fh.tag[2] != 'S') {
		printf("Found invalid file header in file \'%s\', aborting.\n", filename);
		return;
	}

	if (fh.version != '8') {
		printf("Found invalid file version \'%c\' in file \'%s\', should be \'8\', aborting.\n", fh.version, filename);
		return;
	}

	// read time header
	fread(&th, sizeof(th), 1, f);
	if (memcmp(th.time_header, "TIME", 4) != 0) {
		printf("Invalid time header in file \'%s\', aborting.\n", filename);
		return;
	}

	for (b = 0;; b++) {
		// read board header
		fread(&bh, sizeof(bh), 1, f);
		if (memcmp(bh.bn, "B#", 2) != 0) {
			// probably event header found
			fseek(f, -4, SEEK_CUR);
			break;
		}

		printf("Found data for board #%d\n", bh.board_serial_number);

		// read time bin widths
		memset(bin_width[b], sizeof(bin_width[0]), 0);
		for (chn = 0; chn < 18; chn++) {
			fread(&ch, sizeof(ch), 1, f);
			if (ch.c[0] != 'C') {
				// event header found
				fseek(f, -4, SEEK_CUR);
				break;
			}
			i = (ch.cn[1] - '0') * 10 + ch.cn[2] - '0';
			printf("Found timing calibration for channel #%d\n", i);
			fread(&bin_width[b][i][0], sizeof(float), 1024, f);
		}
	}
	n_boards = b;
	RunHeader runHeader(
			n_boards, 18,/////
			runnum, eh.second);
	WaveData waveData(runHeader);
	__uint16_t nX742PointsPerCh =
	  runHeader.getNumberOfX742PointsPerChannel();
	__uint16_t* chdata0 =
	  new __uint16_t[nX742PointsPerCh];
	//	  __uint32_t* chdata1;
	rec = new TTree("data", "");
	rec->Branch("x", &waveData, 64000);
	waveData.setEventNumber(evnum);

	// loop over all events in data file
	for (n = 0; n<max_Ev ; n++) {

		// read event header
		i = fread(&eh, sizeof(eh), 1, f);
		if (i < 1)
			break;
		//if( n > 1000) break;



		printf("Found event #%d\r", eh.event_serial_number);
		fflush(stdout);
//		rec = new TTree(Form("Event%d", eh.event_serial_number), Form("Event%d", eh.event_serial_number));
     
		// loop over all boards in data file
		for (b = 0; b < n_boards; b++) {

			// read board header
			fread(&bh, sizeof(bh), 1, f);
			if (memcmp(bh.bn, "B#", 2) != 0) {
				printf("Invalid board header in file \'%s\', aborting.\n", filename);
				return;
			}

			if (n_boards > 1)
				printf("Found data for board #%d\n", bh.board_serial_number);

			// reach channel data
			for (chn=0 ; chn<18 ; chn++) {

				// read channel header
				fread(&ch, sizeof(ch), 1, f);
				if (ch.c[0] == 'E') {
					// event header found
					fseek(f, -4, SEEK_CUR);
					break;
				}
				chn_index = (ch.cn[1] - '0')*10 + ch.cn[2] - '0';

				if (ch.c[0] == 'C') {

					// Read DRS data
					fread(&scaler, sizeof(int), 1, f);

					// read trigger cell
					fread(&tch, sizeof(tch), 1, f);
					if (memcmp(tch.tc, "T#", 2) != 0) {
						printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
						return;
					}

					fread(voltage, sizeof(unsigned short), 1024, f);
					for (i = 0; i < 1024; i++) {
						// convert data to volts
						waveform[b][chn_index][i] = (voltage[i] / 65536. + eh.range / 1000.0 - 0.5);
//						voltage_2[i]=(voltage[i]/ 65536. + eh.range / 1000.0 - 0.5);
//						std::cout<<"wf value "<<waveform[b][chn_index][i]<<std::endl;
//						std::cout <<eh.range<<std::endl;
						// calculate time for this cell
						for (j = 0, time[b][chn_index][i] = 0; j < i; j++)
							time[b][chn_index][i] += bin_width[b][chn_index][(j + tch.trigger_cell) % 1024];
					}
				} else if (ch.c[0] == 'A') {

					// Read ADC data
					short adc_voltage[2048];
					fread(adc_voltage, sizeof(short), 2048, f);
					for (i = 0; i < 2048; i++) {
						// convert data to volts
						adc_waveform[b][chn_index][i] = (adc_voltage[i] / 65536. + eh.range / 1000.0 - 0.5);
					}

				} else if (ch.c[0] == 'T') {

					if (ch.cn[0] == '0') {
						// Read TDC
						fread(tdc_waveform[b][chn_index], sizeof(unsigned char), 512, f);
					} else {
						// Read TRG
						fread(trg_data[b], sizeof(unsigned long), 512, f);
					}
				}
				waveData.setX742ChannelData(
						chn, voltage);
			}

			// align cell #0 of all channels
			t1 = time[b][0][(1024 - tch.trigger_cell) % 1024];
			for (chn = 1; chn < 18; chn++) {
				t2 = time[b][chn][(1024 - tch.trigger_cell) % 1024];
				dt = t1 - t2;
				for (i = 0; i < 1024; i++){
					time[b][chn][i] += dt;
					voltage_2[i]=(waveform[b][chn][i]);
//					voltage_2[i]-=0.45;
//					voltage_2[i]*=-1;

//					std::cout<<voltage_2[i]<<" "<<voltage_2[i]-0.45<<std::endl;
//					std::cout<<"chn "<<chn<<" wf value "<<waveform[0][chn][i]<<std::endl;
				}
//				waveData.setX742ChannelData(
//						chn, voltage_2);
//				waveData.setX742ChannelData(
//						chn, voltage);

			}


			/*
         You need to insert here the splitting for the binary file. 
         In this way the C streamer can navigate through the file and be able to find the right event number.
			 */

			// draw graph and wait for user click

			//c1->Update();
			//gPad->WaitPrimitive();

		}

		//if(eh.event_serial_number  < (10000 * (n_batch-1) ) ) continue;
		//else if (eh.event_serial_number >= (10000 * n_batch)) break;

//		for (chn=0; chn<16; chn++)
//		{
////			//c1->cd(chn+1);
////			//g[chn]=new TGraph(1024, (double *) time[0][chn], (double *) waveform[0][chn]);
////			//
////			//for (i = 0; i < 1024; i++)
////			//{
////			//   g[chn]->SetPoint(i, time[0][chn][i], waveform[0][chn][i]);
////			//}
////			//
////			//g[chn]->SetTitle(Form("Wafeform %d",chn));
////			//g[chn]->GetXaxis()->SetTitle("Time");
////			//g[chn]->GetYaxis()->SetTitle("Voltage(V)");
////			//g[chn]->Draw("ACP");
////			//gianl
//////			rec->Branch(Form("Wafeform%d",chn),&waveform[0][chn], Form("Wafeform [1024]%d/D",chn));
//////			rec->Branch(Form("Time%d",chn), &time[0][chn], Form("Time[1024]%d/D",chn));
////			//
////			//rec->Branch(Form("Graph%d",chn),"TGraph",g[chn],32000,0);
////			//rec->Draw(Form("Graph%d",chn));
//	        waveData.setX742ChannelData(
//	        		chn, voltage);
////
//		}
		rec->Fill();
	}
	rec->SetDirectory(0);

	printf("\n");
	rec->Write();
	outfile->Close();
	fclose(f);
	return;
	// project tree
	//rec->Draw("amp0");
}

int main(int argc, char* argv[]) {
  int runnum = atoi(argv[1]);
  int Nev=atoi(argv[2]);
  std::cout<<runnum<<" "<<Nev<<std::endl;
  decode(argv[3], Nev, runnum);
  return 0;
}
