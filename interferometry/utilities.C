#include "./IceRayTracing/namespace/wROOT/IceRayTracing.cc"

// #include "AraEventCalibrator.h"
#include "/users/PAS0654/jflaherty13/source/AraRoot/AraRoot_build/include/AraEventCalibrator.h"
#include <glob.h>

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/pool/poolfwd.hpp>

#include "TTree.h"

#include <ROOT/RDataFrame.hxx>  //Adding to fix compile error - JCF 9/22/2022

using Mult_t = boost::multiprecision::cpp_dec_float_100;

double prop_time_s(double tx[3], double rx[3]){

        double x1 = sqrt(pow(tx[0] - rx[0], 2) + pow(tx[1] - rx[1], 2));
        return IceRayTracing::IceRayTracing(0, tx[2], x1, rx[2])[6];

}

std::vector<Mult_t>::iterator search(std::vector<Mult_t> &vec, double value){

        auto it = lower_bound(vec.begin(), vec.end(), value);

        if(it != vec.begin())
                if(abs(value - *(it - 1)) <  abs(value - *it))
                        --it;

        return it;

}


// ARA Stuff

// Given L1 directory, gives list of all run files. Unsafe.
vector<std::string> getRunTable(const char *input_path, const char* type = ""){

        glob_t glob_result;

        char *wild_card;
        asprintf(&wild_card, "%s/*/*/*/event0*.root", input_path);

        glob(wild_card, GLOB_TILDE, NULL, &glob_result);

        vector<std::string> table;

        for(int i = 0; i < glob_result.gl_pathc; ++i)
                table.push_back(glob_result.gl_pathv[i]);

        return table;

}

// Given station, run, returns path to run file. Unsafe.
std::string getRunFile(const char* L1_path, int station, int run){

        auto run_table = getRunTable(L1_path);

        for(auto& run_file: run_table)
                if(run_file.find(Form("ARA0%d", station)) != std::string::npos &&
                   run_file.find(Form("%d", run)) != std::string::npos)
                        return run_file;

}

//Adding a definition of rms function, as it seems to be undeclared.  Need to hear back fron Keith on the proper function definition. - JCF 9/20/2022
double rms(double* x, int startingIndex, int endingIndex)
{
	double sum = 0;
    int n = 0;

	for (int i = startingIndex; i < endingIndex; i++)
		sum += pow(x[i], 2);
        n++;

	return sqrt(sum / n);
}
//rms()

// Simple SNR calculation. Quarter waveform, define snr as Vpeak/Lowest Quarter RMS. Unsafe. Should check if no channels are empty.

vector<double> snr(vector<TGraph*>& waveform){

	vector<double> channels;

	for(auto& channel: waveform){

		// auto volt_array = channel->GetY();  //Commenting out for debugging - JCF 9/20/2022
        auto array = channel->GetY();  //Changing "volt_array" to "array" to try fixing a compile error. - JCF 9/20/2022
		int segment = channel->GetN()/4;

		double RMS = DBL_MAX;

		for(int s = 0; s < 4; ++s)
			RMS = std::min(RMS, rms(array, segment * s, segment * (s+1)));

		// chennels.push_back(*std::max_element(array, array + (segment * 4)) / RMS);  //Pretty sure this should be "channels" and not "chennels". - JCF 9/20/2022
		channels.push_back(*std::max_element(array, array + (segment * 4)) / RMS);        

	}

	return channels;

}
