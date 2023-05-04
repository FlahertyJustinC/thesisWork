#include <glob.h>
#include <iostream>
#include <fstream>

void recursive(const char* statn, const char* year, const char* out_path){

  glob_t glob_result;
 
  glob(Form("/data/exp/ARA/%s/unblinded/L1/ARA%s/*/run*/event0*.root", year,statn),GLOB_TILDE,NULL,&glob_result);

  ofstream outfile(out_path);

  for(int i = 0; i<glob_result.gl_pathc; ++i){
    outfile << glob_result.gl_pathv[i] << "\n";
    cout <<  glob_result.gl_pathv[i] << "\n";		
  }

}
