#include "utilities.C"

void time_cut(const char* data_path, double tx[3], double rxi[3], double rxj[3], double rxk[3], double dt = 0){

	ROOT::RDataFrame depth_data("N","./SPICE_depth.root");
	auto pulser_depth = depth_data.Graph("unix","depth");

	ROOT::RDataFrame coincident_events("D", data_path);

	vector<vector<int>> ctr(160, vector<int>(3));

	double depth;
	int depth_bin;

	coincident_events.Foreach([&](double ai, double aj, double ak){

		depth = pulser_depth->Eval(ai);
		depth_bin = floor(abs(depth)/10);
		
		// double ai_pred = prop_time_s({tx[0], tx[1], depth}, rxi),
		//        aj_pred = prop_time_s({tx[0], tx[1], depth}, rxj),
		//        ak_pred = prop_time_s({tx[0], tx[1], depth}, rxk);
        
        //Reorganizing this for debugging - JCF 9/23/20222
        double txTemp[3] = {tx[0], tx[1], depth};
        
        double ai_pred = prop_time_s(txTemp, rxi),
		       aj_pred = prop_time_s(txTemp, rxj),
		       ak_pred = prop_time_s(txTemp, rxk);
        //debugging

		if(abs(ai-aj) < abs(ai_pred - aj_pred) + dt) 
			++ctr[depth_bin][0];

		if(abs(ai-ak) < abs(ai_pred - ak_pred) + dt) 
			++ctr[depth_bin][1];

		if(abs(aj-ak) < abs(aj_pred - ak_pred) + dt)
			++ctr[depth_bin][2];

		if(abs(ai-aj) < abs(ai_pred - aj_pred) + dt &&
		   abs(ai-ak) < abs(ai_pred - ak_pred) + dt &&
		   abs(aj-ak) < abs(aj_pred - ak_pred) + dt)
			++ctr[depth_bin][3];

	},{"ai","aj","ak"});

	// Do plotting

}
