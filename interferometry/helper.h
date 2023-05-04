std::vector<double> get_detector_cog(Detector *detectorPtr){
    std::vector<double> average_position{0., 0., 0.};
    int n_ant = 0;
    for(int s=0; s<4; s++){
        cout << s << endl;
        for(int a=0; a<4; a++){
            cout << a << endl;
            n_ant++;
            average_position[0]+=detectorPtr->stations[0].strings[s].antennas[a].GetX();
            average_position[1]+=detectorPtr->stations[0].strings[s].antennas[a].GetY();
            average_position[2]+=detectorPtr->stations[0].strings[s].antennas[a].GetZ();
        }
    }
    for(int i=0; i<n_ant; i++){
        cout << i << endl;
        average_position[i]/=double(n_ant);
    }
    return average_position;
}


int guess_triggering_solution(Event *eventPtr, Report *reportPtr){

    // guess the triggering solution
    // start by getting all of the view angles
    double changle_deg = TMath::RadToDeg()*eventPtr->Nu_Interaction[0].changle;
    std::vector< std::vector<double> > view_angles;
    for(int s=0; s<4; s++){
        for(int a=0; a<4; a++){
            std::vector<double> these_view_angs;
            for(int v=0; v<reportPtr->stations[0].strings[s].antennas[a].view_ang.size(); v++){
                double the_view_angle = reportPtr->stations[0].strings[s].antennas[a].view_ang[v];
                these_view_angs.push_back(the_view_angle);
            }
            view_angles.push_back(these_view_angs);
        }
    }
    int count_dir = 0;
    int count_ref = 0;
    for(int ant=0; ant < view_angles.size(); ant++){
        int num_sols = view_angles[ant].size();
        if(num_sols!=0 && num_sols!=2){
            continue;
        }
        if(num_sols==2){
            double dir_angle = changle_deg - TMath::RadToDeg()*view_angles[ant][0];
            double ref_angle = changle_deg - TMath::RadToDeg()*view_angles[ant][1];
            // printf("  Ant %d, Dir Angle %.2f, Ref Angle %.2f\n",ant, dir_angle, ref_angle);
            if(abs(dir_angle) < abs(ref_angle)){
                count_dir++;
            }
            else if(abs(ref_angle) < abs(dir_angle)){
                count_ref++;
            }
        }
    }
    int likely_sol = 0;
    if(count_ref > count_dir){
        likely_sol = 1;
    }
    return likely_sol;
}

double get_value(string which, Report *reportPtr, int station, int string, int antenna, int ray){
    if(which =="theta"){
        return reportPtr->stations[station].strings[string].antennas[antenna].theta_rec[ray];
    }
    else if(which =="phi"){
        return reportPtr->stations[station].strings[string].antennas[antenna].phi_rec[ray];
    }
    else{
        throw std::invalid_argument("An unsupported argument for 'which' was supplied.");
    }
}

std::map<int, double> get_value_from_mc_truth(string which, int solution, Report *reportPtr){
    // now, we need to yank out the true arrival angles for the "likely" solution
    std::map<int, double> retvals;
    for(int s=0; s<4; s++){
        for(int a=0; a<4; a++){
            int idx = (s*4) + a;
            int num_sol = reportPtr->stations[0].strings[s].antennas[a].arrival_time.size();
            if(num_sol!=0 && num_sol!=2){
                std::cerr<<"Something is up with MC truth for index "<<idx<<". Skipping."<<std::endl;
                continue;
            }
            else if(num_sol==2){
                double retval = get_value(which, reportPtr, 0, s, a, solution);
                retvals[idx] = retval;
                // printf("String %d, Antenna %d, idx %d, retval %f \n", s, a, idx, retval);
            }
        }
    }
    return retvals;
}

void rerange_theta_phi(double &theta, double &phi){

        // fix the theta range
        if(theta < 0){
            theta = abs(theta) + 90;
        }
        else if(theta>0){
            theta = 90 - theta;
        }
        else{
            theta+=90;
        }

        if(phi < 0) phi =360.0 - abs(phi);
}