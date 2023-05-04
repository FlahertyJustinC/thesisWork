#include <boost/multiprecision/cpp_dec_float.hpp>

#include <boost/pool/pool_alloc.hpp>
#include <boost/pool/poolfwd.hpp>

using Mult_t = boost::multiprecision::cpp_dec_float_100;

void writeTimes(const char* in_path){

	ifstream stream(Form("%s.list", in_path));

	vector<string> run_list;
	string line;

	while(getline(stream, line))
		run_list.emplace_back(line);

	ROOT::RDataFrame df("eventTree", run_list);

	unordered_map<int, unordered_map<int, int>> pps_offset;

	df.Foreach([&](ULong64_t unixTime, UInt_t ppsNumber, Int_t run){
		++pps_offset[run][unixTime - ppsNumber];
	},{"unixTime", "ppsNumber", "run"});

	df.Define("eventTime", [&](UInt_t timeStamp, UInt_t ppsNumber, int run){
		ppsNumber += max_element(begin(pps_offset[run]), end(pps_offset[run]))->first;
		return Mult_t(ppsNumber + Mult_t(timeStamp)/1e+8);
	},{"timeStamp", "ppsNumber", "run"}).Snapshot("T", Form("%s.root", in_path), {"eventTime","run","eventNumber"});

}
