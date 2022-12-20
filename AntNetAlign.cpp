/***************************************************************************
                           AntNetAlign.cpp  -  description
                             -------------------
    begin                : Wed Oct 21 2020
    copyright            : (C) 2020 by Guillem Rodriguez
    email                : guillem.rodriguez.corominas@upc.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <sys/time.h>
#include <sys/resource.h>
#include <chrono>
#include <stack>
#include <queue>
using namespace std;

random_device rd;
std::mt19937 generator(rd());
std::uniform_real_distribution<double> dis01(0.0, 1.0);

//Parameters (and default values)
int n_of_ants = 10; // 1 <= n_ants
double tau_max = 0.999; // 0 < tau_min <= tau_max < 1
double tau_min = 0.001;
double positive_learning_rate = 0.3; // 0 < learning_rate <= 0.5
double d_rate_align = 0.9, d_rate_select = 0.8; //0 <= d_rate <= 1
int max_constructions = 1000;
double max_running_time = 3600; //Maximum running time in seconds (Default: 1h)
string score = "S3";
bool tuning = false;

//Fixed parameters (and default values)
double cf = 0.0;
int n_topo_iter = 4;
double epsilon = 1e-3;

//Instance data
int V1, V2;
int E1, E2;
vector< set<int> > G1, G2;
vector<vector<bool>> AdjMat_G1, AdjMat_G2; 

//Similarity metrics
vector<vector<double>> similarityMatrix;
vector<double> centralityG1;

//Mapping strings to integers
map<string,int> string2vertexG1;
vector<string> vertex2stringG1;
map<string,int> string2vertexG2;
vector<string> vertex2stringG2;

//Input files
bool save = false;
string g1_file, g2_file, similarity_file = "", output_file = "output", output_dir = "."; 



// ----- UTILITY ----- //

struct Option {
	int vertex;
	double value;
};

struct Solution {
	map<string,double> scores;
	double time;
	vector<int> alignment;
};

inline int stoi(string &s) {
	return atoi(s.c_str());
}

inline double stof(string &s) {
	return atof(s.c_str());
}

double getCurrentTime(){
	struct rusage res;
	getrusage(RUSAGE_SELF, &res);
  	return	(double) res.ru_utime.tv_sec + (double) res.ru_stime.tv_sec + 
			(double) res.ru_utime.tv_usec * 1.0E-6 + (double) res.ru_stime.tv_usec * 1.0E-6;
}

double magnitude(const vector<double> &v){
	double m = 0;
	for(int i = 0; i < v.size(); ++i) m += v[i]*v[i];
	return sqrt(m);
}

double dot(const vector<double> &v1, const vector<double> &v2){
	double d = 0;
	for(int i = 0; i < v1.size(); ++i) d += v1[i]*v2[i];
	return d;
}

int calculate_induced_G2(const vector<int>& alignment){
	int induced_g2 = 0;
	for(int i = 0; i < V1; ++i){
		int j = alignment[i];
		for(int i2 = i+1; i2 < V1; ++i2){
			int j2 = alignment[i2];
			if(G2[j].find(j2) != G2[j].end()) ++induced_g2;
		}
	}
	return induced_g2;
}

// ----- INPUT -----//

void read_parameters(int argc, char **argv) {
	int iarg=1;
	while (iarg < argc) {
		if (strcmp(argv[iarg],"-g1")==0) g1_file = (argv[++iarg]);
		else if (strcmp(argv[iarg],"-g2")==0) g2_file = (argv[++iarg]);
		else if (strcmp(argv[iarg],"-similarity")==0) similarity_file = (argv[++iarg]);
		else if (strcmp(argv[iarg],"-max_constructions")==0) max_constructions = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-n_ants")==0) n_of_ants = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pos_l_rate")==0) positive_learning_rate = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-d_rate_select")==0) d_rate_select = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-d_rate_align")==0) d_rate_align = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-score")==0) score = (argv[++iarg]);
		else if (strcmp(argv[iarg],"-max_t")==0) max_running_time = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-tuning")==0) tuning = true;
		else if (strcmp(argv[iarg],"-save")==0) save = true;
		else if (strcmp(argv[iarg],"-output_file")==0) output_file = (argv[++iarg]);
		else if (strcmp(argv[iarg],"-output_dir")==0) output_dir = (argv[++iarg]);
		iarg++;
	}
}

//Read input data
void read_input(){
	
	//Read G1
	ifstream indata;
	indata.open(g1_file);
	if(!indata) { //File couldn't be opened
		cout << "Error: G1 file could not be opened" << endl;
	}

	string2vertexG1 = map<string,int>();
	vertex2stringG1 = vector<string>();
	G1 = vector<set<int>>();

	string p1, p2;

	int current_vertex = 0;
	E1 = 0;
	while(indata >> p1 >> p2){
		if(string2vertexG1.find(p1) == string2vertexG1.end()) {
			string2vertexG1.insert({p1,current_vertex});
			vertex2stringG1.push_back(p1);
			G1.push_back(set<int>());
			current_vertex += 1;
		}
		if(string2vertexG1.find(p2) == string2vertexG1.end()) {
			string2vertexG1.insert({p2,current_vertex});
			vertex2stringG1.push_back(p2);
			G1.push_back(set<int>());
			current_vertex += 1;
		}
		G1[string2vertexG1[p1]].insert(string2vertexG1[p2]);
		G1[string2vertexG1[p2]].insert(string2vertexG1[p1]);
		E1++;
	}
	
	V1 = current_vertex;
	AdjMat_G1 = vector<vector<bool>> (V1,vector<bool> (V1,false));
	for(int i = 0; i < V1; ++i){
		for(set<int>::iterator it = G1[i].begin(); it != G1[i].end(); ++it){
			AdjMat_G1[i][*it] = true;
		}
	}
	
	if(!tuning) cout << "[G1] Nodes: " << V1 << " Edges: " << E1 << endl;
	
	indata.close();
	
	

	//Read G2
	indata.open(g2_file);
	if(!indata) { //File couldn't be opened
		cout << "Error: G2 file could not be opened" << endl;
	}

	string2vertexG2 = map<string,int>();
	vertex2stringG2 = vector<string>();
	G2 = vector<set<int>>();
	
	current_vertex = 0;
	E2 = 0;
	while (indata >> p1 >> p2){
		if(string2vertexG2.find(p1) == string2vertexG2.end()) {
			string2vertexG2.insert({p1,current_vertex});
			vertex2stringG2.push_back(p1);
			G2.push_back(set<int>());
			current_vertex += 1;
		}
		if(string2vertexG2.find(p2) == string2vertexG2.end()) {
			string2vertexG2.insert({p2,current_vertex});
			vertex2stringG2.push_back(p2);
			G2.push_back(set<int>());
			current_vertex += 1;
		}
		G2[string2vertexG2[p1]].insert(string2vertexG2[p2]);
		G2[string2vertexG2[p2]].insert(string2vertexG2[p1]);
		E2++;
	}
	
	V2 = current_vertex;
	AdjMat_G2 = vector<vector<bool>> (V2,vector<bool> (V2,false));
	for(int i = 0; i < V2; ++i){
		for(set<int>::iterator it = G2[i].begin(); it != G2[i].end(); ++it){
			AdjMat_G2[i][*it] = true;
		}
	}
	
	if(!tuning) cout << "[G2] Nodes: " << V2 << " Edges: " << E2 << endl << endl;
	

	if(V1 > V2) { //File couldn't be opened
		cout << "Error: first network must be smaller or equal than the second" << endl;
	}

	indata.close();
}

void read_similarity(){
	ifstream indata;
	indata.open(similarity_file);
	if(!indata) { //File couldn't be opened
		cout << "Error: Similarity file could not be opened" << endl;
	}
	for(int i = 0; i < V1; ++i){
		for(int j = 0; j < V2; ++j){
			indata >> similarityMatrix[i][j];
		}
	}
	indata.close();
}

// ----- OUTPUT ----- //

//Displays the alignment
void print_solution(Solution& sol){
	cout << "EC\t" << sol.scores["EC"];
	cout << "ICS\t" << sol.scores["ICS"];
	cout << "S3\t" << sol.scores["S3"];
 	for(int i = 0; i < V1; ++i){
		cout << vertex2stringG1[i] << "\t" << vertex2stringG2[sol.alignment[i]] << endl;
	}
}

//Write the results in the output files
void save_results(Solution& sol){
	ofstream outdata;
	outdata.open(output_dir+"/"+output_file+".info");
	outdata << "S3\t" << sol.scores["S3"] << endl;
	outdata << "EC\t" << sol.scores["EC"] << endl;
	outdata << "ICS\t" << sol.scores["ICS"] << endl;
	outdata << "Time\t" << sol.time << endl;
	outdata.close();
	outdata.open(output_dir+"/"+output_file+".alignment");
	for(int i = 0; i < V1; ++i){
		outdata << vertex2stringG1[i] << "\t" << vertex2stringG2[sol.alignment[i]] << endl;
	}
	outdata.close();
}

// ----- ACO ----- //

//Generates a solution given a positive pheromone matrix
void generate_aco_solution(Solution& sol, vector<vector<double>>& positive_pheromones) {

	sol.scores["EC"] = 0;
	sol.scores["ICS"] = 0;
	sol.scores["S3"] = 0;
	sol.alignment = vector<int>(V1);

	int num_unaligned_vertices = V1;

	//Set of non-aligned vertices from G1
	set<int> unaligned_g1;
	for (int i = 0; i < V1; ++i) unaligned_g1.insert(i);

	//Boolean vector telling if a vertex in G1 is aligned
	vector<bool> is_aligned_g1(V1,false);
    
	//Set of non-aligned vertices from G2
	set<int> unaligned_g2;
	for (int i = 0; i < V2; ++i) unaligned_g2.insert(i);
	
	//Boolean vector telling if a vertex in G2 is aligned
	vector<bool> is_aligned_g2(V2,false);

	//Information about induced and conserved edges
	vector<int> num_induced_G1(V1,0);
	vector<int> num_induced_G2(V2,0);
	vector<vector<int>> num_conserved(V1,vector<int>(V2,0));
		
	int total_conserved_edges = 0;
	int total_induced_edges = 0;

	double total_centrality_sum = 0;
	for(int i = 0; i < V1; ++i) total_centrality_sum += centralityG1[i];

	double total_induced_sum = total_centrality_sum;

	while (num_unaligned_vertices > 0) {

		//Find next vertex to align
		int current_vertex = 0;
		//Selection by roullette wheel 
		if (total_induced_sum == total_centrality_sum or dis01(generator) > d_rate_select) {
			set<int>::iterator it = unaligned_g1.begin();
			double wheel = 0.0;
			double rand = dis01(generator);
			while ((wheel < rand) and (it != unaligned_g1.end())) {
				wheel += ((double)num_induced_G1[*it]+centralityG1[*it])/total_induced_sum;
				++it;
			}
			it--;
			current_vertex = *it;
		}

		//Deterministic selection
		else {
			int max_val = -1.0;
			int n_max_can = 0;
			for(set<int>::iterator it = unaligned_g1.begin(); it != unaligned_g1.end(); ++it){
				int current_value = num_induced_G1[*it];
				if(current_value >= max_val){
					if(current_value > max_val){
						max_val = current_value;
						n_max_can = 0;
						current_vertex = *it;
					}
					n_max_can++;
					if(dis01(generator) <= 1.0/(double)n_max_can) current_vertex = *it;
				}	
			}
		}

		//Calcuating the candidates
		vector<Option> choices;
		double option_sum = 0.0;
		double max_value = -1.0;
		Option max_candidate;
		int number_max_candidates = 0;
		Option chosen_option;
		int chosen_vertex;
		
		for (set<int>::iterator cit = unaligned_g2.begin(); cit != unaligned_g2.end(); ++cit) {
			
			double improvement = ((double)num_conserved[current_vertex][*cit]+epsilon)/(double)(num_induced_G1[current_vertex]+num_induced_G2[*cit]-num_conserved[current_vertex][*cit]+epsilon);
			double h_func = improvement*similarityMatrix[current_vertex][*cit]; //Greedy function

			Option opt;
			opt.vertex = *cit;
			opt.value = h_func*positive_pheromones[current_vertex][*cit];
			option_sum += opt.value;
			if (opt.value >= max_value) {
				if (opt.value > max_value) {
					max_value = opt.value;
					number_max_candidates = 0;
					max_candidate = opt;
				}
				number_max_candidates++;
				if(dis01(generator) <= 1.0/(double)number_max_candidates) max_candidate = opt;
			}
			choices.push_back(opt);
		}

		//SELECT NODE TO ALIGN
		//Selection by roullette wheel 
		if (dis01(generator) > d_rate_align) {
			vector<Option>::iterator oit = choices.begin();
			double wheel = 0.0;
			double rand = dis01(generator);
			while ((wheel < rand) and (oit != choices.end())) {
				wheel += (*oit).value/option_sum;
				++oit;
			}
			oit--;
			chosen_option = *oit;
		}

		//Deterministic selection
		else chosen_option = max_candidate;

		chosen_vertex = chosen_option.vertex;

		//Update data
		num_unaligned_vertices -= 1;
		sol.alignment[current_vertex] = chosen_vertex;

		//Remove aligned nodes
		unaligned_g1.erase(current_vertex);
		is_aligned_g1[current_vertex] = true;
		total_induced_sum-=((double)num_induced_G1[current_vertex]+centralityG1[current_vertex]);
		total_centrality_sum-= centralityG1[current_vertex];

		unaligned_g2.erase(chosen_vertex);
		is_aligned_g2[chosen_vertex] = true;

		//Update alignment values (conserved and induced edges)
		total_conserved_edges += num_conserved[current_vertex][chosen_vertex];
		total_induced_edges += num_induced_G1[current_vertex]+num_induced_G2[chosen_vertex];

		if(num_unaligned_vertices == 0) break;

		for(set<int>::iterator nit1 = G1[current_vertex].begin(); nit1 != G1[current_vertex].end(); ++nit1) {
			num_induced_G1[*nit1]++;
			if(!is_aligned_g1[*nit1]) total_induced_sum++;
		}
		for(set<int>::iterator nit2 = G2[chosen_vertex].begin(); nit2 != G2[chosen_vertex].end(); ++nit2) num_induced_G2[*nit2]++;
		for(set<int>::iterator nit1 = G1[current_vertex].begin(); nit1 != G1[current_vertex].end(); ++nit1)
			for(set<int>::iterator nit2 = G2[chosen_vertex].begin(); nit2 != G2[chosen_vertex].end(); ++nit2)
				num_conserved[*nit1][*nit2]++;
	}
	
	//Calculate scores
	sol.scores["EC"] = (double)total_conserved_edges/(double)E1;
	sol.scores["ICS"] = (double)total_conserved_edges/(double)(total_induced_edges-E1);
	sol.scores["S3"] = (double)total_conserved_edges/(double)(total_induced_edges-total_conserved_edges);
}	

double compute_convergence_factor(vector<vector<double>>& pheromone) {
	double ret_val = 0.0;
	for (int i = 0; i < V1; ++i) {
		for(int j = 0; j < V2; ++j) {
        	if ((tau_max - pheromone[i][j]) > (pheromone[i][j] - tau_min)) ret_val += tau_max - pheromone[i][j];
        	else ret_val += pheromone[i][j] - tau_min;
		}
	}
	int ph_cnt = V1*V2;
	ret_val = ret_val / (double(ph_cnt) * (tau_max - tau_min));
	ret_val = (ret_val - 0.5) * 2.0;
	return ret_val;
}


// -----TOPOLOGICAL SIMILARITY-----//

//Compute similarity between v1 in G1 and v2 in G2 given S 
double vertex_similarity(int v1, int v2, const vector<vector<double>> & S){
	int d1 = G1[v1].size(); //Degree v1
	int d2 = G2[v2].size(); //Degree v2

	
	if(d1 == 0 and d2 == 0) return 1;
	if(d1 == 0 or d2 == 0) return 0;

	double sum = 0;
	if (d1 <= d2){
		set<int> candidates = G2[v2];

		for(set<int>::iterator it1 = G1[v1].begin(); it1 != G1[v1].end(); ++it1){
			set<int>::iterator chosen = candidates.begin();
			double best_similarity = S[*it1][*candidates.begin()];
			for(set<int>::iterator it2 = candidates.begin(); it2 != candidates.end(); ++it2){
				if(S[*it1][*it2] > best_similarity){
					best_similarity = S[*it1][*it2];
					chosen = it2;
				}
			}
			candidates.erase(chosen);
			sum += best_similarity;
		}
	}
	else {
		set<int> candidates = G1[v1];

		for(set<int>::iterator it2 = G2[v2].begin(); it2 != G2[v2].end(); ++it2){
			set<int>::iterator chosen = candidates.begin();
			double best_similarity = S[*candidates.begin()][*it2];
			for(set<int>::iterator it1 = candidates.begin(); it1 != candidates.end(); ++it1){
				if(S[*it1][*it2] > best_similarity){
					best_similarity = S[*it1][*it2];
					chosen = it1;
				}
			}
			candidates.erase(chosen);
			sum += best_similarity;
		}
	}
	
	if(d1 >= d2) return sum /= (double)d1;
	else return sum /= (double)d2;
}




// -----CENTRALITY-----//
void compute_centrality(){

	if(!tuning) cout << "-> Computing centrality..." << endl;
	
	//Brandes Algorithm
	centralityG1 = vector<double>(V1,0);
	
	for (int i = 0; i < V1; ++i){
		stack<int> S;
		vector<list<int>> P(V1, list<int>());
		vector<int> sigma(V1,0); sigma[i] = 1;
		vector<int> d(V1,-1); d[i] = 0;
		queue<int> Q; Q.push(i);
		while (!Q.empty()) {
			int j = Q.front();
			Q.pop();
			S.push(j);
			for(set<int>::iterator it = G1[j].begin(); it != G1[j].end(); ++it){
				int k = *it;
				if (d[k] < 0){
					Q.push(k);
					d[k] = d[j] + 1;
				}
				if(d[k] == d[j] + 1){
					sigma[k] += sigma[j];
					P[k].push_back(j);
				}
			}
		}
		
		vector<int> delta(V1,0);
		while(!S.empty()){
			int k = S.top();
			S.pop();
			for (list<int>::iterator it = P[k].begin(); it != P[k].end(); ++it){
				int j = *it;
				delta[j] += (sigma[j]/sigma[k])*(1+delta[k]);
			}
			if(k != i) centralityG1[k] += delta[k];
		}
	}

	//Normalize centralities
	double max = 0.0;
	for(int i = 0; i < V1; ++i)
		if(centralityG1[i] > max) max = centralityG1[i];
	for(int i = 0; i < V1; ++i)
		centralityG1[i] /= max;
}

// ----- SIMILARITY-----//

//Compute similarity matrix between G1 and G2
void compute_similarity(){
	
	if(!tuning) cout << "-> Computing topological similarity..." << endl;
	
	
	similarityMatrix = vector<vector<double>>(V1,vector<double>(V2));
	
	if(similarity_file != "") read_similarity();
	else{
		//Construct Simiarity Matrix
		vector<vector<double>> temporalMatrix(V1, vector<double>(V2));

		//Initialize with degree difference
		for(int i = 0; i < V1; ++i)
			for(int j = 0; j < V2; ++j)
				similarityMatrix[i][j] = (min((double)G1[i].size(),(double)G2[j].size())+1)/(max((double)G1[i].size(),(double)G2[j].size())+1);

		for(int it = 0; it < n_topo_iter; it++){
			for(int i = 0; i < V1; i++) {
				for(int j = 0; j < V2; j++) {
					if(it%2 == 0) temporalMatrix[i][j] = vertex_similarity(i,j,similarityMatrix);
					else similarityMatrix[i][j] = vertex_similarity(i,j,temporalMatrix);
				}
			}
		}
		if(n_topo_iter%2 != 0){
			for(int i = 0; i < V1; i++) {
				for(int j = 0; j < V2; j++) {
					similarityMatrix[i][j] = temporalMatrix[i][j];
				}
			}
		}

		//Normalize
		for(int i = 0; i < V1; i++){
			double max = 0.0;
			for(int j = 0; j < V2; j++)
				if(similarityMatrix[i][j] > max) max = similarityMatrix[i][j];
			if(max > 0) 
				for(int j = 0; j < V2; j++)
					similarityMatrix[i][j] /= max;
		}
	}
}


// ----- MAIN ----- //

int main( int argc, char **argv ) {

	if(argc < 3) {
		cout << "Use: na_aco_induced -input <input_file>  ..." << endl;
		exit(1);
	}
	else read_parameters(argc,argv);

	//Set decimal precision to 5
	std::cout << std::setprecision(5) << std::fixed; 
	
	if(!tuning) cout << "--- INITIALIZATION ---" << endl << endl;

	//Read data from the input file
	read_input();

	//Compute similarity matrix
	compute_similarity();
	
	//Compute node betweenness centrality on the source network
	compute_centrality();

	//Time starts running!
	double startingTime = getCurrentTime();

	//Initialization of the pheromone values. Note that, in the is version of ACO, all pheromone values tau move in 0 < tau_min <= tau <= tau_max < 1
	vector<vector<double>> positive_pheromones(V1,vector<double>(V2,0.5));

	//Initalization Global-Best Solution
	Solution bestSol;
	bestSol.scores = map<string,double>();
	bestSol.scores["EC"] = std::numeric_limits<double>::lowest();
	bestSol.scores["ICS"] = std::numeric_limits<double>::lowest();
	bestSol.scores["S3"] = std::numeric_limits<double>::lowest();

	//Initialization Restart-Best Solution
	Solution restartSol;
	restartSol.scores = map<string,double>();
	restartSol.scores["EC"] = std::numeric_limits<double>::lowest();
	restartSol.scores["ICS"] = std::numeric_limits<double>::lowest();
	restartSol.scores["S3"] = std::numeric_limits<double>::lowest();

	bool global_convergence = false;

	double result, time;
	int n_constructions = 0; //Number of constructed solutions
	int current_it = 0;
	
	//Main loop
	while (n_constructions < max_constructions && (getCurrentTime() - startingTime) < max_running_time) {
	
		if(!tuning) cout << endl << "--- ITERATION " << current_it << " ---" << endl << endl;
	
		double iteration_average = 0.0;
		//Initialization Iteration-Best Solution
		Solution iBestSol;
		iBestSol.scores = map<string,double>();
		iBestSol.scores["EC"] = std::numeric_limits<double>::lowest();
		iBestSol.scores["ICS"] = std::numeric_limits<double>::lowest();
		iBestSol.scores["S3"] = std::numeric_limits<double>::lowest();
		
		//Generate n_of_ants solutions
		n_of_ants = min(n_of_ants,max_constructions-n_constructions);
		for (int na = 0; na < n_of_ants; ++na) {
		
			//Initialization Iteration Solution
			Solution iSol;
			iSol.scores = map<string,double>();
			iSol.scores["EC"] = std::numeric_limits<double>::lowest();
			iSol.scores["ICS"] = std::numeric_limits<double>::lowest();
			iSol.scores["S3"] = std::numeric_limits<double>::lowest();
	
			//Generate ACO solution
			generate_aco_solution(iSol, positive_pheromones);
			if(!tuning) cout << "[" << na+1 << "]\tS3: " << iSol.scores["S3"] << "\tEC: " << iSol.scores["EC"] << "\tICS: " << iSol.scores["ICS"] << endl;
			iteration_average += iSol.scores[score];
			if (iSol.scores[score] > iBestSol.scores[score]) iBestSol = iSol;
			if (iSol.scores[score] > restartSol.scores[score]) restartSol = iSol;
			if (iSol.scores[score] > bestSol.scores[score]) {
				bestSol = iSol;
				result = bestSol.scores[score];
				time = getCurrentTime() - startingTime;
				bestSol.time = time;
			}
		}
		n_constructions += n_of_ants;
		iteration_average /= double(n_of_ants);

		//The following three weights determine the influence of iBestSol, restartSol and bestSol on the pheromone update
		double i_weight, r_weight, g_weight;
		if (global_convergence) {
			i_weight = 0.0;
			r_weight = 0.0;
			g_weight = 1.0;
		}
		else {
			if (cf < 0.2) {
				i_weight = 1.0;
				r_weight = 0.0;
				g_weight = 0.0;
			}
			else if (cf < 0.5) {
				i_weight = 2.0/3.0;
				r_weight = 1.0/3.0;
				g_weight = 0.0;
			}
			else if (cf < 0.8) {
				i_weight = 1.0/3.0;
				r_weight = 2.0/3.0;
				g_weight = 0.0;
			}
			else {
				i_weight = 0.0;
				r_weight = 1.0;
				g_weight = 0.0;
			}
		}

		//Calculate the contribution for the positive pheromone update 
		vector<vector<double>> contribution(V1,vector<double>(V2,0.0));
		if (i_weight > 0.0) {
			for (int i = 0; i < V1; ++i) contribution[i][iBestSol.alignment[i]] += i_weight;
		}
		if (r_weight > 0.0) {
			for (int i = 0; i < V1; ++i) contribution[i][restartSol.alignment[i]] += r_weight;
		}
		if (g_weight > 0.0) {
			for (int i = 0; i < V1; ++i) contribution[i][bestSol.alignment[i]] += g_weight;
		}

		//Positive pheromone update using the contributions and taking into account the positive learning rate
		for (int i = 0; i < V1; ++i) {
			for(int  j = 0; j < V2; ++j){
				positive_pheromones[i][j] += (positive_learning_rate * (contribution[i][j] - positive_pheromones[i][j]));
				if (positive_pheromones[i][j] > tau_max) positive_pheromones[i][j] = tau_max;
				if (positive_pheromones[i][j] < tau_min) positive_pheromones[i][j] = tau_min;
			}
		}

		//Compute convergence factor
		cf = compute_convergence_factor(positive_pheromones);
		if(!tuning) cout << endl << "[STATS]\tBest: " << result << "\tTime: " << time << "\tCF: " << cf << "\t" << "It.Average: " << iteration_average << endl;

		//Algorithm has converged (cf > 0.9999)
		if (cf > 0.9999) {
			cf = 0;
			if (global_convergence) {
				global_convergence = false;
				positive_pheromones = vector<vector<double>>(V1,vector<double>(V2,0.5));
				restartSol.scores["EC"] = std::numeric_limits<double>::lowest();
				restartSol.scores["ICS"] = std::numeric_limits<double>::lowest();
				restartSol.scores["S3"] = std::numeric_limits<double>::lowest();
			}
			else global_convergence = true;
		}
		
		current_it++;
	}
	if(tuning) cout << -result << endl;
	else {
		cout << endl << "--- BEST SOLUTION ---" << endl << endl;
		cout << "S3\t" << bestSol.scores["S3"] << endl;
		cout << "EC\t" << bestSol.scores["EC"] << endl;
		cout << "ICS\t" << bestSol.scores["ICS"] << endl;
		cout << "Time\t" << time << endl;
	}
	if(save) save_results(bestSol);
}

