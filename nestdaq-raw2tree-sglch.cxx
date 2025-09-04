#include <iostream>
#include <inttypes.h>
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>

#include "FileSinkHeader.h"
#include "FileSinkTrailer.h"
#include "AmQStrTdcData.h"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "FilterHeader.h"
#include "TFile.h"
#include "TTree.h"

int hbf_sorting (std::ifstream &ifs, uint64_t max_num_read_hbf,
		 std::map<uint32_t, uint64_t>& hbf_sorted_address, int64_t &read_pos){
  if(read_pos == -1){ // read_pos == -1: the whole of the file was fully read.
    return 0;
  }
  ifs.seekg(read_pos, std::ios_base::beg);
  uint32_t currTimeFrameId = 0;
  uint32_t file_read_flag = 1;
  while(file_read_flag == 1){
    uint64_t magic;
    ifs.read((char*)&magic, sizeof(magic));
    if (ifs.eof()) {
      read_pos = -1;
      ifs.clear();
      break;
    }
    ifs.seekg(-sizeof(magic), std::ios_base::cur);
    switch (magic) {
    case TimeFrame::MAGIC: {
      TimeFrame::Header tfbHeader;
      ifs.read((char*)&tfbHeader, sizeof(tfbHeader));
      if (tfbHeader.timeFrameId != currTimeFrameId) {
	if (hbf_sorted_address.size() == max_num_read_hbf) {
	  read_pos = (uint64_t) ifs.tellg() - sizeof(TimeFrame::Header);
	  file_read_flag = 0;
	  break;
	}
	currTimeFrameId = tfbHeader.timeFrameId;
	hbf_sorted_address[currTimeFrameId] = (uint64_t) ifs.tellg() - sizeof(TimeFrame::Header);
	// std::cout << "TimeFrameId: " << std::dec << tfbHeader.timeFrameId << std::hex << " 0x" << tfbHeader.timeFrameId << std::endl;
      }
      break;}
    case SubTimeFrame::MAGIC: {
      SubTimeFrame::Header stfHeader;
      ifs.read((char*)&stfHeader, sizeof(stfHeader));
      unsigned int nword = (stfHeader.length - sizeof(stfHeader)) / 8;
      ifs.seekg(nword * sizeof(AmQStrTdc::Data::Bits), std::ios_base::cur);
      break;}
    case FileSinkTrailer::MAGIC: {
      ifs.seekg(sizeof(FileSinkTrailer::Trailer), std::ios_base::cur);
      break;}
    case Filter::MAGIC: {
      ifs.seekg(sizeof(Filter::Header), std::ios_base::cur);
      break;}
    default: {
      ifs.read((char*)&magic, sizeof(magic));
      uint32_t length;
      ifs.read((char*)&length, sizeof(length));
      ifs.seekg(length - sizeof(length) - sizeof(magic), std::ios_base::cur);
      break;}
    }
  }
  //for (auto ite = hbf_sorted_address.begin(); ite != hbf_sorted_address.end(); ite++) {
  //  std::cout << "Key = " << ite->first << ", Value = " << ite->second << std::endl;
  //}
  return 0;
}

int main(int argc, char* argv[]){
  if (argc <= 2) {
    std::cout << "Usage: ./nestdaq-raw2tree-sglch input-filename output-filename" << std::endl;
    return 1;
  }
  std::string   filename = argv[1];
  std::ifstream ifs(filename.c_str(), std::ios::binary);
  ifs.seekg(0, std::ios::end);
  std::ifstream::pos_type fsize = ifs.tellg();
  ifs.seekg(0, std::ios_base::beg);
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Data file: "<< filename
	    << ", size: " << fsize / 1024 / 1024
	    << " MB" << std::endl;
  
  std::string   rootfile = argv[2];
  TFile * of = TFile::Open(rootfile.c_str(),"RECREATE");
  std::cout << "ROOT file: "<< rootfile << std::endl;
  TTree *tr = new TTree("tr","tr");
  Long64_t rawtdc, rawtot, rawhbfn, hbfn;
  Double_t tdc;
  tr->Branch("rawtdc",  &rawtdc,  "rawtdc/L");
  tr->Branch("rawtot",  &rawtot,  "rawtot/L");
  tr->Branch("rawhbfn", &rawhbfn, "rawhbfn/L");
  tr->Branch("hbfn",    &hbfn,    "hbfn/L");
  tr->Branch("tdc",     &tdc,     "tdc/D");
  
  FileSinkHeader::Header fileHeader;
  ifs.read((char*)&fileHeader,sizeof(fileHeader));
  std::map<uint32_t, uint64_t> hbf_sorted_address;
  int64_t read_pos = ifs.tellg();
  uint64_t max_num_read_hbf = 20000;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Number of HBF for initial sorting: " << max_num_read_hbf << std::endl;
  std::cout << "Now initial HBF sorting is starting..." << std::endl;
  hbf_sorting(ifs, max_num_read_hbf, hbf_sorted_address, read_pos);
  std::cout << "Initial sorting was finished!" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Then, starting analysis..." << std::endl;
  
  uint64_t selectedFemId = 0xc0a802a9; // IP 192.168.2.169 was selected
  uint64_t selectedCh    = 33;         // IP 192.168.2.169, ch33 was selected.
  uint64_t currentFemId  = 0xffffffff;
  uint64_t currentCh     = 0xffff;
  
  uint32_t timeFrameChangedFlag  = 1;
  uint32_t currTimeFrameId       = 0;
  int64_t  savedCurrPos          = 0;
  int64_t  hbfn0                 = -1;
  int64_t  hbfnPrev              = -1;
  int64_t  rawtdcPrev            = -1;
  uint64_t hbfnCarryFlag         = 0;
  uint64_t hit_counter           = 0;
  uint64_t hit_counter_no_double = 0;
  
  uint64_t i = 0;
  //while( true && i < 1000000){
  while( true ){
    i++;
    savedCurrPos = ifs.tellg();
    if (timeFrameChangedFlag == 1){
      if (!hbf_sorted_address.empty()) {
	savedCurrPos = hbf_sorted_address.begin()->second;
	currTimeFrameId = hbf_sorted_address.begin()->first;
	hbf_sorted_address.erase(hbf_sorted_address.begin());
	hbf_sorting(ifs, max_num_read_hbf, hbf_sorted_address, read_pos);
      }else{
	break;
      }
      timeFrameChangedFlag = 0;
    }
    ifs.seekg(savedCurrPos, std::ios_base::beg);
    
    uint64_t magic;
    ifs.read((char*)&magic, sizeof(magic));
    if (ifs.eof()) {break;}
    ifs.seekg(-sizeof(magic), std::ios_base::cur);
    switch (magic) {
    case TimeFrame::MAGIC: {
      TimeFrame::Header tfbHeader;
      ifs.read((char*)&tfbHeader, sizeof(tfbHeader));
      std::ifstream::pos_type pos = ifs.tellg();
      static int prev_read_ratio = 0;
      int curr_read_ratio = (double) pos / fsize * 10;
      if (  curr_read_ratio > prev_read_ratio ) {
	std::cout << "Read: " <<  curr_read_ratio * 10 << " % ("<< pos/1024/1024 << " MB / " << fsize/1024/1024 << " MB)" << std::endl;
	prev_read_ratio = curr_read_ratio;
      }
      if (hbfn0 == -1) {
	hbfn0 = tfbHeader.timeFrameId;
      }
      if (currTimeFrameId != tfbHeader.timeFrameId) {
	timeFrameChangedFlag = 1;
	break;
      }
      break;}
    case SubTimeFrame::MAGIC: {
      SubTimeFrame::Header stfHeader;
      ifs.read((char*)&stfHeader, sizeof(stfHeader));
      unsigned int nword = (stfHeader.length - sizeof(stfHeader)) / 8;
      for(unsigned int i=0; i< nword; i++){
	AmQStrTdc::Data::Bits idata;
	ifs.read((char*)&idata, sizeof(idata));
	if (idata.head == AmQStrTdc::Data::Heartbeat) {
	  if ((stfHeader.femId == selectedFemId)
	      && (currentFemId == selectedFemId)
	      && (currentCh == selectedCh)){
	    if ((hbfnPrev - idata.hbframe) > 0x10000) {
	      hbfnCarryFlag++;
	    }
	    if ( (rawtdc != rawtdcPrev)
		 || (idata.hbframe != hbfnPrev) ) {
	      rawhbfn = idata.hbframe;
	      hbfn = idata.hbframe + hbfnCarryFlag * 0xffffff - hbfn0; // hbfn0: first heart beat frame number 
	      tdc = hbfn * 524288.0 + rawtdc / 1024.;                  // 1 hbf = 0.524288 msec, unit of the parameter "tdc" is nsec
	      tr->Fill();
	      hit_counter_no_double++;
	      //std::cout << "FemId: 0x" << std::hex << std::setw(8) << std::setfill('0') << stfHeader.femId << std::setfill(' ') << std::dec;
	      //std::cout << ", hbfn: " << hbfn
	      //	      << ", rawhbfn: " << rawhbfn
	      //	      << ", rawtdc:  " << rawtdc
	      //	      << ", rawtot:  " << rawtot
	      //	      << ", hbfnCarryFlag:  " << hbfnCarryFlag
	      //	      << ", tdc: " << tdc << std::endl;
	    }
	    rawtdcPrev   = rawtdc;
	    hbfnPrev     = idata.hbframe;
	    currentFemId = 0xffffffff;
	    currentCh    = 0xffff;
	  }
	}else if (idata.head == AmQStrTdc::Data::Data){
	  if ( stfHeader.femType == 2 || stfHeader.femType == 5 ) { // HRTDC
	    if ((stfHeader.femId == selectedFemId)
		&& (idata.hrch == selectedCh)){
	      rawtdc = idata.hrtdc;
	      rawtot = idata.hrtot;
	      hit_counter++;
	    }
	    currentFemId = stfHeader.femId;
	    currentCh    = idata.hrch;
	  }else if ( stfHeader.femType == 3 || stfHeader.femType == 6 ) { // LRTDC
	    if ((stfHeader.femId == selectedFemId)
		&& (idata.ch == selectedCh)){
	      rawtdc = idata.tdc;
	      rawtot = idata.tot;
	      hit_counter++;
	    }
	    currentFemId = stfHeader.femId;
	    currentCh = idata.ch;
	  }
	}
      }
      break;}
    case FileSinkTrailer::MAGIC: {
      ifs.seekg(sizeof(FileSinkTrailer::Trailer), std::ios_base::cur);
      break;}
    case Filter::MAGIC: {
      ifs.seekg(sizeof(Filter::Header), std::ios_base::cur);
      break;}
    default: {
      ifs.read((char*)&magic, sizeof(magic));
      uint32_t length;
      ifs.read((char*)&length, sizeof(length));
      ifs.seekg(length - sizeof(length) - sizeof(magic), std::ios_base::cur);
      break;}
    }
  }
  tr->Write();
  of->Close();
  std::cout << "For IP address: 0x" << std::hex << std::setw(8) << std::setfill('0') << selectedFemId << std::setfill(' ') << std::dec
	    << ", ch: " << selectedCh << std::endl;
  std::cout << "Hit count: " << hit_counter << std::endl;
  std::cout << "Hit count (no double count): " << hit_counter_no_double << std::endl;
  return 0;
}
