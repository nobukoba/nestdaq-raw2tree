#include <iostream>
#include <inttypes.h>
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <map>
#include <chrono>

#include "FileSinkHeader.h"
#include "FileSinkTrailer.h"
#include "AmQStrTdcData.h"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"
#include "FilterHeader.h"
#include "TFile.h"
#include "TTree.h"

int read_tf (std::ifstream &ifs, uint64_t max_num_read_tf,
	     std::map<uint32_t, std::vector<char> >& sorted_time_frame_data){
  uint32_t currTimeFrameId = 0x1000000;
  uint32_t file_read_flag = 1;
  while(file_read_flag == 1){
    TimeFrame::Header tfbHeader;
    ifs.read((char*)&tfbHeader, sizeof(tfbHeader));
    if (ifs.eof()) { break; }
    switch (tfbHeader.magic) {
    case TimeFrame::MAGIC: {
      if (tfbHeader.timeFrameId != currTimeFrameId) {
	if (sorted_time_frame_data.size() == max_num_read_tf) {
	  ifs.seekg(-sizeof(TimeFrame::Header), std::ios_base::cur);
	  file_read_flag = 0;
	  break;
	}
	currTimeFrameId = tfbHeader.timeFrameId;
	sorted_time_frame_data[currTimeFrameId] = std::vector<char>(tfbHeader.length);
	std::memcpy(sorted_time_frame_data[currTimeFrameId].data(), &tfbHeader, sizeof(TimeFrame::Header));
	ifs.read(sorted_time_frame_data[currTimeFrameId].data() + sizeof(TimeFrame::Header),
		 tfbHeader.length - sizeof(TimeFrame::Header));
	//std::cout << "TimeFrameId: " << std::dec << tfbHeader.timeFrameId << std::hex << " 0x" << tfbHeader.timeFrameId << std::endl;
      }else{
	std::ifstream::pos_type offset = sorted_time_frame_data[currTimeFrameId].size();
	sorted_time_frame_data[currTimeFrameId].resize(offset + (std::ifstream::pos_type) tfbHeader.length);
	std::memcpy(sorted_time_frame_data[currTimeFrameId].data() + offset, &tfbHeader, sizeof(TimeFrame::Header));
	ifs.read(sorted_time_frame_data[currTimeFrameId].data() + offset + sizeof(TimeFrame::Header),
		 tfbHeader.length - sizeof(TimeFrame::Header));
      }
      break;}
    default: {
      ifs.seekg( -sizeof(TimeFrame::Header) + sizeof(FileSinkTrailer::Trailer), std::ios_base::cur);
      break;}
    }
  }
  return 0;
}

int fill_tree( std::map<int, Long64_t>& map_rawtdc,
	       std::map<int, Long64_t>& map_rawtot,
	       uint64_t* selectedChs,
	       int64_t& hbfn0, Long64_t& rawhbfn, Long64_t& hbfn,
	       Long64_t* rawtdc, Long64_t* rawtot,
	       Double_t* tdc, Double_t*tot,
	       TTree* tr){
  static Long64_t prev_rawhbfn = 0;
  static Long64_t prev_rawtdc  = 0;
  if(map_rawtdc.count(selectedChs[0])==0){ return 0; }
  if ((map_rawtdc[selectedChs[0]] == prev_rawtdc) &&
      (rawhbfn == prev_rawhbfn)) {
    return 0;
  }
  prev_rawtdc = map_rawtdc[selectedChs[0]];
  prev_rawhbfn = rawhbfn;
  
  for (int ich = 0; ich < 9; ich++) {
    uint64_t ch = selectedChs[ich];
    if (map_rawtdc.count(ch) != 0){
      rawtdc[ich] = map_rawtdc[ch];
      rawtot[ich] = map_rawtot[ch];
      tdc[ich] = hbfn * 524288.0 + rawtdc[ich] / 1024.;       // 1 hbf = 0.524288 msec, unit of the parameter "tdc" is nsec
      tot[ich] = rawtot[ich] / 1024.;                         // 1 LSB = 1 ns / 1024 ~ 0.9766 ps
    }else{
      rawtdc[ich] = -1000;
      tdc[ich]    = -1000;
      rawtot[ich] = -1000;
      tot[ich]    = -1000;
    }
  }
  tr->Fill();
  map_rawtdc.clear();
  map_rawtot.clear();
  hbfn = -1;
  rawhbfn = -1;
  return 0;
}

int print_read_ratio(std::ifstream &ifs){
  static bool first_call_flag = true;
  static std::ifstream::pos_type fsize = 0;
  static std::chrono::time_point first_call_time = std::chrono::system_clock::now();
  std::ifstream::pos_type pos = ifs.tellg();
  if (first_call_flag == true) {
    ifs.seekg(0, std::ios::end); 
    fsize = ifs.tellg();
    ifs.seekg(pos, std::ios_base::beg);
    first_call_flag = false;
  }
  static int prev_read_ratio = 0;
  int curr_read_ratio = (double) pos / fsize * 10;
  if (  curr_read_ratio > prev_read_ratio ) {
    std::chrono::time_point curr_time = std::chrono::system_clock::now();
    std::chrono::seconds elapsed_time = std::chrono::duration_cast<std::chrono::seconds> (curr_time - first_call_time);
    std::cout << std::dec << "Read: " <<  curr_read_ratio * 10 << " % ("<< pos/1024/1024 << " MB / " << fsize/1024/1024 << " MB), "
	      << std::dec << "elapsed time: " <<  elapsed_time.count() << " sec"
	      << std::endl;
    prev_read_ratio = curr_read_ratio;
  }
  return 0;
}

int main(int argc, char* argv[]){
  if (argc <= 2) {
    std::cout << "Usage: ./nestdaq-raw2tree-sglch input-filename output-filename" << std::endl;
    return 1;
  }
  std::string   filename = argv[1];
  std::ifstream ifs(filename.c_str(), std::ios::binary | std::ios::ate); // std::ios::ate is "at end"?, equivalent to ifs.seekg(0, std::ios::end) 
  std::ifstream::pos_type fsize = ifs.tellg();
  ifs.seekg(0, std::ios_base::beg);
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Data file: "<< filename
	    << ", size: " << fsize / 1024 / 1024
	    << " MB" << std::endl;
  
  std::string rootfile = argv[2];
  TFile * of = TFile::Open(rootfile.c_str(),"RECREATE");
  std::cout << "ROOT file: "<< rootfile << std::endl;
  TTree *tr = new TTree("tr","tr");
  Long64_t rawhbfn, hbfn;
  Long64_t rawtdc[9], rawtot[9];
  Double_t tdc[9], tot[9];
  std::map<int, Long64_t> map_rawtdc;
  std::map<int, Long64_t> map_rawtot;
  tr->Branch("rawhbfn", &rawhbfn, "rawhbfn/L");
  tr->Branch("hbfn",    &hbfn,    "hbfn/L");
  tr->Branch("rawtdc",  rawtdc,   "rawtdc[9]/L");
  tr->Branch("rawtot",  rawtot,   "rawtot[9]/L");
  tr->Branch("tdc",     tdc,      "tdc[9]/D");
  tr->Branch("tot",     tot,      "tot[9]/D");
  
  FileSinkHeader::Header fileHeader;
  ifs.read((char*)&fileHeader,sizeof(fileHeader));
  std::map<uint32_t, std::vector<char> > sorted_time_frame_data;
  uint64_t max_num_read_tf = 2000;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Number of time frames for initial sorting: " << max_num_read_tf << std::endl;
  std::cout << "Now initial time frame sorting is starting..." << std::endl;
  read_tf(ifs, max_num_read_tf, sorted_time_frame_data);
  std::cout << "Initial sorting was finished! Num. of sorted tf: " << sorted_time_frame_data.size() << std::endl;
  std::cout << "sizeof(sorted_time_frame_data): " << sizeof(sorted_time_frame_data) << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  std::cout << "Then, starting analysis..." << std::endl;
  
  uint64_t selectedChs[9] = {33,10,11,12,13,42,43,44,45}; // IP 192.168.2.169, ch33 was selected.
  
  int64_t  hbfn0                  = -1;
  uint64_t hbfnCarryFlag          = 0;
  uint32_t currTimeFrameId        = 0x1000000;
  
  while( !sorted_time_frame_data.empty() ){
    print_read_ratio(ifs);
    sorted_time_frame_data.erase(sorted_time_frame_data.begin());
    if ( sorted_time_frame_data.empty() ) { break; }
    read_tf(ifs, max_num_read_tf, sorted_time_frame_data);
    char* ptr = sorted_time_frame_data.begin()->second.data();
    char* end = sorted_time_frame_data.begin()->second.data() + sorted_time_frame_data.begin()->second.size();
    while (ptr < end) {
      uint64_t magic = *reinterpret_cast<const uint64_t*>(ptr);
      switch (magic) {
      case TimeFrame::MAGIC: {
	TimeFrame::Header tfbHeader = *reinterpret_cast<const TimeFrame::Header*>(ptr);
	ptr += sizeof(TimeFrame::Header);
	if (hbfn0 == -1) {
	  hbfn0 = tfbHeader.timeFrameId;
	}
	if (currTimeFrameId == 0x1000000) {
	  currTimeFrameId = tfbHeader.timeFrameId;
	}else if (currTimeFrameId != tfbHeader.timeFrameId) {
	  if (tfbHeader.timeFrameId < currTimeFrameId) {
	    hbfnCarryFlag++;
	  }
	  currTimeFrameId = tfbHeader.timeFrameId;
	}
	//std::cout << "TimeFrameId: " << std::dec << tfbHeader.timeFrameId << std::hex << " 0x" << tfbHeader.timeFrameId << std::endl;
	break;}
      case SubTimeFrame::MAGIC: {
	rawhbfn = -1;
	hbfn = -1;
	map_rawtdc.clear();
	map_rawtot.clear();
	SubTimeFrame::Header stfHeader = *reinterpret_cast<const SubTimeFrame::Header*>(ptr);
	ptr += sizeof(SubTimeFrame::Header);
	unsigned int nword = (stfHeader.length - sizeof(stfHeader)) / 8;
	for(unsigned int i=0; i< nword; i++){
	  AmQStrTdc::Data::Bits idata = *reinterpret_cast<const AmQStrTdc::Data::Bits*>(ptr);
	  ptr += sizeof(AmQStrTdc::Data::Bits);
	  if (idata.head == AmQStrTdc::Data::Heartbeat) {
	    rawhbfn = idata.hbframe;
	    hbfn = idata.hbframe + hbfnCarryFlag * 0x1000000 - hbfn0; // hbfn0: first heart beat frame number 
	    //std::cout << "FemId: 0x" << std::hex << std::setw(8) << std::setfill('0') << stfHeader.femId << std::setfill(' ') << std::dec;
	    //std::cout << ", hbfn: " << hbfn
	    //	      << ", rawhbfn: " << rawhbfn
	    //	      << ", rawtdc:  " << rawtdc
	    //	      << ", rawtot:  " << rawtot
	    //	      << ", hbfnCarryFlag:  " << hbfnCarryFlag
	    //	      << ", tdc: " << tdc << std::endl;
	  }else if (idata.head == AmQStrTdc::Data::Data){
	    if ( stfHeader.femType == 2 || stfHeader.femType == 5 ) { // HRTDC
	      if (map_rawtdc.count(idata.hrch) == 0){
		map_rawtdc[idata.hrch] = idata.hrtdc;
	        map_rawtot[idata.hrch] = idata.hrtot;
	      }
	    }else if ( stfHeader.femType == 3 || stfHeader.femType == 6 ) { // LRTDC
	    }
	  }
	}
	fill_tree(map_rawtdc, map_rawtot, selectedChs,
		  hbfn0, rawhbfn, hbfn, rawtdc, rawtot, tdc, tot, tr);
	break;}
      case FileSinkTrailer::MAGIC: {
	ptr += sizeof(FileSinkTrailer::Trailer);
	break;}
      case Filter::MAGIC: {
	ptr += sizeof(Filter::Header);
	break;}
      default: {
	magic = *reinterpret_cast<const uint64_t*>(ptr);
	ptr += sizeof(uint64_t);
	uint32_t length = *reinterpret_cast<const uint32_t*>(ptr);
	ptr += length - sizeof(magic);
	break;}
      }
    }
  }
  tr->Write();
  of->Close();
  return 0;
}
