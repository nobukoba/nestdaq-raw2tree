#include <iostream>
#include <inttypes.h>
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
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
  static uint32_t currTimeFrameId = 0;
  uint32_t file_read_flag = 1;
  while(file_read_flag == 1){
    TimeFrame::Header tfbHeader;
    ifs.read((char*)&tfbHeader, sizeof(tfbHeader));
    if (ifs.eof() || (ifs.tellg() == -1)) { break; }
    switch (tfbHeader.magic) {
    case TimeFrame::MAGIC: {
      if (tfbHeader.timeFrameId != currTimeFrameId) {
	if (sorted_time_frame_data.size() == max_num_read_tf) {
	  ifs.seekg(-sizeof(TimeFrame::Header), std::ios_base::cur);
	  file_read_flag = 0;
	  break;
	}
	if ((tfbHeader.timeFrameId + 0x800000) < currTimeFrameId) {
	  uint32_t max_out_counter = 0;
	  while ((tfbHeader.timeFrameId + 0x1000000 * max_out_counter + 0x800000) < currTimeFrameId ) {
	    max_out_counter++;
	  }
	  currTimeFrameId = tfbHeader.timeFrameId + 0x1000000 * max_out_counter;
	}else{
	  currTimeFrameId = tfbHeader.timeFrameId;
	}
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
  ULong64_t rawhbfn, rawtdc, rawtot,  hbfn;
  Double_t tdc, tot;
  std::map<uint64_t, std::unordered_multimap<uint64_t, ULong64_t> > map_rawtdc;
  std::map<uint64_t, std::unordered_multimap<uint64_t, ULong64_t> > map_rawtot;
  tr->Branch("rawhbfn", &rawhbfn, "rawhbfn/l");
  tr->Branch("rawtdc",  &rawtdc,  "rawtdc/l");
  tr->Branch("rawtot",  &rawtot,  "rawtot/l");
  tr->Branch("hbfn",    &hbfn,    "hbfn/l");
  tr->Branch("tdc",     &tdc,     "tdc/D");
  tr->Branch("tot",     &tot,     "tot/D");
  
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
  
  uint64_t selectedFemId = 0xc0a802a9; // IP 192.168.2.169 was selected
  uint64_t selectedCh    = 33;         // IP 192.168.2.169, ch33 was selected.
  
  uint64_t hbfn0                  = 0;
  uint64_t rawhbfnPrev            = 0;
  uint64_t rawtdcPrev             = 0;
  //uint64_t hbfnCarryFlag          = 0;
  uint64_t hit_counter_all        = 0;
  uint64_t hit_counter            = 0;
  uint64_t hit_counter_no_double  = 0;
  //uint32_t currTimeFrameId = 0x1000000;
  
  while( !sorted_time_frame_data.empty() ){
    print_read_ratio(ifs);
    sorted_time_frame_data.erase(sorted_time_frame_data.begin());
    if ( sorted_time_frame_data.empty() ) { break; }
    read_tf(ifs, max_num_read_tf, sorted_time_frame_data);
    char* ptr = sorted_time_frame_data.begin()->second.data();
    char* end = sorted_time_frame_data.begin()->second.data() + sorted_time_frame_data.begin()->second.size();
    uint64_t max_out_counter_of_TimeFrameId = 0;
    while (ptr < end) {
      uint64_t magic = *reinterpret_cast<const uint64_t*>(ptr);
      switch (magic) {
      case TimeFrame::MAGIC: {
	TimeFrame::Header tfbHeader = *reinterpret_cast<const TimeFrame::Header*>(ptr);
	ptr += sizeof(TimeFrame::Header);
	if (hbfn0 == 0) {
	  hbfn0 = tfbHeader.timeFrameId;
	}
	while ( sorted_time_frame_data.begin()->first !=
		(tfbHeader.timeFrameId + max_out_counter_of_TimeFrameId * 0x1000000)) {
	  max_out_counter_of_TimeFrameId++;
	}
	//std::cout << "TimeFrameId: " << std::dec << tfbHeader.timeFrameId << std::hex << " 0x" << tfbHeader.timeFrameId << std::endl;
	break;}
      case SubTimeFrame::MAGIC: {
	rawhbfn = 0;
	hbfn = 0;
	for (auto it = map_rawtdc.begin(); it != map_rawtdc.end(); ++it){
	  it->second.clear();
	}
	for (auto it = map_rawtot.begin(); it != map_rawtot.end(); ++it){
	  it->second.clear();
	}
	SubTimeFrame::Header stfHeader = *reinterpret_cast<const SubTimeFrame::Header*>(ptr);
	ptr += sizeof(SubTimeFrame::Header);
	unsigned int nword = (stfHeader.length - sizeof(stfHeader)) / 8;
	for(unsigned int i=0; i< nword; i++){
	  AmQStrTdc::Data::Bits idata = *reinterpret_cast<const AmQStrTdc::Data::Bits*>(ptr);
	  ptr += sizeof(AmQStrTdc::Data::Bits);
	  if (idata.head == AmQStrTdc::Data::Heartbeat) {
	    rawhbfn = idata.hbframe;
	    hbfn = rawhbfn + max_out_counter_of_TimeFrameId * 0x1000000 - hbfn0; // hbfn	    //std::cout << "FemId: 0x" << std::hex << std::setw(8) << std::setfill('0') << stfHeader.femId << std::setfill(' ') << std::dec;
	    //std::cout << ", hbfn: " << hbfn
	    //	      << ", rawhbfn: " << rawhbfn
	    //	      << ", rawtdc:  " << rawtdc
	    //	      << ", rawtot:  " << rawtot
	    //	      << ", hbfnCarryFlag:  " << hbfnCarryFlag
	    //	      << ", tdc: " << tdc << std::endl;
	  }else if (idata.head == AmQStrTdc::Data::Data){
	    if (map_rawtdc.count(stfHeader.femId) == 0){
	      map_rawtdc[stfHeader.femId] = std::unordered_multimap<uint64_t, ULong64_t>{};
	      map_rawtot[stfHeader.femId] = std::unordered_multimap<uint64_t, ULong64_t>{};
	    }
	    if ( stfHeader.femType == 2 || stfHeader.femType == 5 ) { // HRTDC
	      if (map_rawtdc[stfHeader.femId].count(idata.hrch) == 0) {
		map_rawtdc[stfHeader.femId].insert({(uint64_t)idata.hrch, (ULong64_t)idata.hrtdc});
		map_rawtot[stfHeader.femId].insert({(uint64_t)idata.hrch, (ULong64_t)idata.hrtot});
	      }
	    }else if ( stfHeader.femType == 3 || stfHeader.femType == 6 ) { // LRTDC
	      if (map_rawtdc[stfHeader.femId].count(idata.ch) == 0) {
		map_rawtdc[stfHeader.femId].insert({(uint64_t)idata.ch, (ULong64_t)idata.tdc});
		map_rawtot[stfHeader.femId].insert({(uint64_t)idata.ch, (ULong64_t)idata.tot});
	      }
	    }
	  }
	}

//	for (auto it = map_rawtdc[selectedFemId].begin(); it != map_rawtdc[selectedFemId].end(); ++it){
//	  std::cout << "it->first: "<< it->first << std::endl;
//	  std::cout << "it->second: "<< it->second << std::endl;
//	}

	if ( (map_rawtdc.count(selectedFemId) > 0) &&
	     (map_rawtdc[selectedFemId].count(selectedCh) > 0) ){
	  rawtdc = map_rawtdc[selectedFemId].find(selectedCh)->second;
	  rawtot = map_rawtot[selectedFemId].find(selectedCh)->second;
	  hit_counter_all += map_rawtdc[selectedFemId].count(selectedCh);
	  hit_counter++;
	  if ((rawhbfn != rawhbfnPrev) ||
	      (rawtdc != rawtdcPrev)) {
	    tdc = hbfn * 524288.0 + rawtdc / 1024.;       // 1 hbf = 0.524288 msec, unit of the parameter "tdc" is nsec
	    tot = rawtot / 1024.;                         // 1 LSB = 1 ns / 1024 ~ 0.9766 ps
	    tr->Fill();
	    hit_counter_no_double++;
	    rawhbfnPrev = rawhbfn;
	    rawtdcPrev  = rawtdc;
	  }
	}
	
	break;}
      case FileSinkTrailer::MAGIC: {
	ptr += sizeof(FileSinkTrailer::Trailer);
	break;}
      case Filter::MAGIC: {
	ptr += sizeof(Filter::Header);
	break;}
      case 0x00454d4954475254: { /* TRGTIME */
	ptr += sizeof(magic);
	uint32_t length, hlength;
	ptr += (sizeof(length) + sizeof(hlength));
	magic = *reinterpret_cast<const uint64_t*>(ptr);
	while (true) { 
	  if (magic == SubTimeFrame::MAGIC) {
	    break;
	  }
	  ptr += sizeof(magic);
	  magic = *reinterpret_cast<const uint64_t*>(ptr);
	}
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
  std::cout << "For IP address: 0x" << std::hex << std::setw(8) << std::setfill('0') << selectedFemId << std::setfill(' ') << std::dec
	    << ", ch: " << selectedCh << std::endl;
  std::cout << "Hit count_all: " << hit_counter_all << std::endl;
  std::cout << "Hit count: " << hit_counter << std::endl;
  std::cout << "Hit count (no double count): " << hit_counter_no_double << std::endl;
  return 0;
}
