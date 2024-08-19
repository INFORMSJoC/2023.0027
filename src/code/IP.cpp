#include <vector>
using namespace std;
#include "gurobi_c++.h"
#include <fstream>
#include <sstream>
#include <time.h>

#include <algorithm>
using namespace std;

class lazy_status
{
public:
	bool TL;
	double solution;
	double bound;
	double total_time_limit;
	clock_t start;
	double number_of_lazy;
	bool take_start;

	lazy_status(bool tl, double tt, clock_t st)
	{
		TL = tl;
		solution = 0;
		bound = 0;
		total_time_limit = tt;
		start = st;
		number_of_lazy = 0;
		take_start = true;
	}
	lazy_status()
	{
		TL = false;
		solution = 0;
		bound = 0;
		total_time_limit = 0;
		start = clock();
		number_of_lazy = 0;
		take_start = true;
	}
};

class vertex
{
public:
	double blocking_cost;
	double vertex_weight;

	vertex()
	{
		blocking_cost = 0;
		vertex_weight = 0;
	}
};

class edge
{
public:
	int head;
	int tail;
	int ID;
	double blocking_costs;
	edge(int hd, int tl, int id, double bc)
	{
		head = hd;
		tail = tl;
		ID = id;
		blocking_costs = bc;
	}

	edge()
	{
		head = 0;
		tail = 0;
		ID = 0;
		blocking_costs = 0;
	}
};

class WGCI_instance
{
public:
	vector <vertex> vertices;
	vector <vector<int>> two_d;
	vector <edge> one_d;
	double desired_clique_weight;
	double gamma;

	WGCI_instance()
	{
		vertices = vector<vertex>();
		two_d = vector <vector<int>>();
		one_d = vector<edge>();
		desired_clique_weight = -1;
		gamma = 0;
	}
};

class gamma_clique
{
public:
	vector <int> vertices;
	vector<int>edge;

	gamma_clique(vector<int>& v, vector<int>& e)
	{
		vertices = v;
		edge = e;
	}

	gamma_clique()
	{
		vertices = vector<int>();
		edge = vector<int>();
	}
};

class node_edge
{
public:
	vector<int> edge_list;
	int last_vertex;
	bool visited;
	node_edge(vector<int>& bs, int lv, bool v)
	{
		edge_list = bs;
		last_vertex = lv;
		visited = v;
	}
	node_edge()
	{
		last_vertex = -1;
		visited = false;
	}
};

bool get_gamma_clique(vector<int>& remain_vertices, vector<bool>& remain_edges, vector<int>& gamma_clique, WGCI_instance& instance, lazy_status& ls);

vector<int> get_A(WGCI_instance& instance, vector<int> vertices, vector<int> candidate);

int get_BIA(WGCI_instance& instance, vector<int> vertices, vector<bool>& remain_edges, vector<int> a);


gamma_clique get_EGC(WGCI_instance& instance, vector<int> vertices, vector<bool>& remain_edges);

class mycallback : public GRBCallback
{
public:
	int numx;
	int numy;
	GRBVar* varsx;
	GRBVar* varsy;
	WGCI_instance instance;
	lazy_status* lazy_s;
	mycallback(int xnumvars, int ny, GRBVar* xvarsx, GRBVar* xvarsy, WGCI_instance& xinstance, lazy_status* ls)
	{
		numx = xnumvars;
		numy = ny;
		varsx = xvarsx;
		varsy = xvarsy;
		instance = xinstance;
		lazy_s = ls;
	}
protected:
	void callback() {
		try {
			if (where == GRB_CB_MIPSOL) {
				// MIP solution callback
				double objbnd = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
				(*lazy_s).bound = objbnd;
				double* x = getSolution(varsx, numx);
				vector<int> remain_vertices;
				for (int i = 0; i < numx; i++)
				{
					if (x[i] < 0.5)
					{//not block
						remain_vertices.push_back(i);
					}
				}

				double* y = getSolution(varsy, numy);
				vector<bool> remain_edges;
				for (int i = 0; i < numy; i++)
				{
					if (y[i] < 0.5)
					{//not block
						remain_edges.push_back(true);
					}
					else
					{
						remain_edges.push_back(false);
					}
				}
				clock_t now = clock();

				vector<int> gamma_c;
				bool reach_time_limit = get_gamma_clique(remain_vertices, remain_edges, gamma_c, instance, *lazy_s);
				if (gamma_c.size() > 0)
				{
					gamma_clique EGC = get_EGC(instance, gamma_c, remain_edges);
					GRBLinExpr lhsx = NULL;
					for (int i = 0; i < EGC.vertices.size(); i++)
					{
						lhsx += varsx[EGC.vertices[i]];
					}
					for (int j = 0; j < EGC.edge.size(); j++)
					{
						lhsx += varsy[EGC.edge[j]];
					}
					addLazy(lhsx >= 1);
					((*lazy_s).number_of_lazy)++;
				}
				else if (reach_time_limit)
				{
					double objbnd = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
					(*lazy_s).bound = objbnd;
					(*lazy_s).TL = true;
					abort();
				}
				else
				{
					double objbst = getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
					(*lazy_s).solution = objbst;
				}
				delete[] x;
				delete[] y;
			}
			else if (where == GRB_CB_MESSAGE) {
				// Message callback
				string msg = getStringInfo(GRB_CB_MSG_STRING);
				if (msg == "User MIP start did not produce a new incumbent solution\n")
				{
					(*lazy_s).take_start = false;
				}
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};

int readdata(WGCI_instance& instance, string filename);

void runmodel(string);

vector<vector<int>> get_eta_sets(WGCI_instance&);

vector<int> get_edges(WGCI_instance& instance, vector<int>& eta_set, vector<bool>& remain_edges);

vector<int> get_candidate_set(WGCI_instance& edge, vector<int>& candidate_set, int target);

int main()
{
	ofstream fout;
	fout.open("WGCI_lazy_UB_E.csv");
	fout << "File Name,|V|,density,Eta,Gamma,Time Limit (s),number_of_lazy,Run Time (s),Gap,Best Sol, take_start, Best Bound, Nodes #,Sol's,removed vertices,removed edgs" << endl;
	fout.close();

	ifstream fin;
	fin.open("filename.txt");
	string str3;
	int number_of_files;
	fin >> number_of_files;
	//finish reading the current line
	getline(fin, str3);
	for (int i = 0; i < number_of_files; i++)
	{
		getline(fin, str3);
		runmodel(str3);
	}
	return 0;
}

void runmodel(string filename)
{
	double total_time_limit = 3600.00;
	WGCI_instance instance;
	int number_of_edges = readdata(instance, filename);

	clock_t pre_process_start = clock();
	//Defining Gurobi environment
	GRBEnv* env = 0;
	env = new GRBEnv();

	// Model
	GRBModel model = GRBModel(*env);
	model.set(GRB_StringAttr_ModelName, filename);

	// Must set LazyConstraints parameter when using lazy constraints
	model.set(GRB_IntParam_LazyConstraints, 1);
	// Decision variables
	GRBVar* x = 0;//edges
	x = model.addVars(instance.vertices.size(), GRB_BINARY);
	GRBVar* y = 0;//eta cliques
	y = model.addVars(instance.one_d.size(), GRB_BINARY);
	for (int i = 0; i < instance.vertices.size(); i++)
	{
		x[i].set(GRB_DoubleAttr_Obj, instance.vertices[i].blocking_cost);
	}
	for (int i = 0; i < instance.one_d.size(); i++)
	{
		y[i].set(GRB_DoubleAttr_Obj, instance.one_d[i].blocking_costs);
	}
	// The objective is minimization
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
	model.update();

	model.set(GRB_IntParam_Threads, 1);
	model.update();

	for (int i = 0; i < instance.one_d.size(); i++)
	{//3b
		model.addConstr(y[i] + x[instance.one_d[i].head], GRB_LESS_EQUAL, 1);
		model.addConstr(y[i] + x[instance.one_d[i].tail], GRB_LESS_EQUAL, 1);
	}

	for (int i = 0; i < instance.vertices.size(); i++)
	{
		x[i].set(GRB_DoubleAttr_Start, 0);
	}
	for (int i = 0; i < instance.one_d.size(); i++)
	{
		y[i].set(GRB_DoubleAttr_Start, 1);
	}

	//set mipfocus=1
	model.set(GRB_IntParam_MIPFocus, 1);

	// Create a callback object and associate it with the model
	lazy_status ls = lazy_status(false, total_time_limit, pre_process_start);
	mycallback cb = mycallback(instance.vertices.size(), instance.one_d.size(), x, y, instance, &ls);

	model.setCallback(&cb);
	clock_t now = clock();
	double limit = total_time_limit - static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC;
	model.set(GRB_DoubleParam_TimeLimit, limit);
	model.optimize();
	clock_t end_time = clock();

	ofstream fout;
	fout.open("WGCI_lazy_UB_E.csv", std::ios_base::app);
	fout << filename << ',' << instance.vertices.size();
	double density = number_of_edges / (static_cast<double>(instance.vertices.size() * (instance.vertices.size() - 1)) / 2);
	fout << ',' << density << ',' << instance.desired_clique_weight << ',' << instance.gamma;
	fout << ',' << total_time_limit << ",";
	fout << ls.number_of_lazy << ",";

	if (ls.TL)
	{//time limit reached in call back
		fout << "LTL" << ',';
		fout << (ls.solution - ls.bound) / ls.solution << ',';
		fout << ls.solution << ',';
		fout << ls.take_start << ',';
		fout << ls.bound << ',';
		fout << model.get(GRB_DoubleAttr_NodeCount) << ',';
		cout << model.get(GRB_DoubleAttr_ObjVal) << ',';
		cout << model.get(GRB_DoubleAttr_ObjBoundC) << ',';
		cout << model.get(GRB_DoubleAttr_NodeCount) << ',';
		fout << endl;
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		fout << static_cast<double>(end_time - pre_process_start) / CLOCKS_PER_SEC << ",";
		fout << "0,";
		fout << model.get(GRB_DoubleAttr_ObjVal) << ',';
		fout << ls.take_start << ',';
		fout << model.get(GRB_DoubleAttr_ObjBoundC) << ',';
		fout << model.get(GRB_DoubleAttr_NodeCount) << ',';
		vector<bool> solutionx, solutiony;
		solutionx.reserve(instance.vertices.size());
		solutiony.reserve(instance.one_d.size());
		int countv = 0;
		int counte = 0;
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			solutionx.push_back(x[i].get(GRB_DoubleAttr_X));
		}
		for (int i = 0; i < instance.vertices.size(); i++)
		{
			if (solutionx[i])
			{
				fout << i + 1 << ';';
				countv++;
			}
		}
		for (int i = 0; i < instance.one_d.size(); i++)
		{
			solutiony.push_back(y[i].get(GRB_DoubleAttr_X));
		}

		for (int i = 0; i < instance.one_d.size(); i++)
		{
			if (solutiony[i])
			{
				fout << '(' << instance.one_d[i].head + 1 << ';' << instance.one_d[i].tail + 1 << ')' << ' ';
				counte++;
			}
		}
		fout << ',' << countv << ',' << counte;
		fout << endl;
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
	{
		fout << "TL" << ',';
		fout << static_cast<double> (model.get(GRB_DoubleAttr_ObjVal) - model.get(GRB_DoubleAttr_ObjBoundC)) / model.get(GRB_DoubleAttr_ObjVal) << ',';
		fout << model.get(GRB_DoubleAttr_ObjVal) << ',';
		fout << ls.take_start << ',';
		fout << model.get(GRB_DoubleAttr_ObjBoundC) << ',';
		fout << model.get(GRB_DoubleAttr_NodeCount) << ',';
		fout << endl;
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
	{
		fout << "IN" << ',';
		fout << '-,';
		fout << '-,';
		fout << '-,';
		fout << model.get(GRB_DoubleAttr_NodeCount) << ',';
		fout << endl;
	}
	fout.close();
}

bool get_gamma_clique(vector<int>& remain_vertices, vector<bool>& remain_edges, vector<int>& gamma_clique, WGCI_instance& instance, lazy_status& ls)
{//return true if reach time limit, else return false
	vector<double> weight;
	for (int i = 0; i < remain_vertices.size(); i++)
	{
		weight.push_back(instance.vertices[remain_vertices[i]].vertex_weight);
	}
	sort(weight.begin(), weight.end());
	double temp_sum_weight = 0;
	int K = 1;
	for (int i = weight.size() - 1; i >= 0; i--)
	{
		temp_sum_weight += weight[i];
		if (temp_sum_weight > instance.desired_clique_weight)
		{
			K = remain_vertices.size() - i;
			break;
		}
	}
	//Defining Gurobi environment
	GRBEnv* sub_env = NULL;
	sub_env = new GRBEnv();

	// Model
	GRBModel sub_model = GRBModel(*sub_env);

	// Decision variables
	GRBVar* y = NULL;
	y = sub_model.addVars(remain_vertices.size(), GRB_BINARY);
	GRBVar* u = NULL;
	u = sub_model.addVars(remain_vertices.size());

	GRBVar* z = NULL;
	z = sub_model.addVars(remain_vertices.size());

	GRBLinExpr lh4b = 0, lh4c = 0, lh4f = 0, lh4g = 0;

	for (int i = 0; i < remain_vertices.size(); i++)
	{
		z[i].set(GRB_DoubleAttr_LB, 0);
		lh4b += instance.vertices[remain_vertices[i]].vertex_weight * y[i];
		lh4c += z[i];
		sub_model.addConstr(z[i] - instance.vertices.size() * y[i], GRB_LESS_EQUAL, 0);
		GRBLinExpr lh4e = z[i];
		for (int j = 0; j < remain_vertices.size(); j++)
		{
			int edge_id = instance.two_d[remain_vertices[i]][remain_vertices[j]];
			if (edge_id != -1 && remain_edges[edge_id])
			{
				lh4e -= y[j];
			}
		}
		sub_model.addConstr(lh4e, GRB_LESS_EQUAL, 0);
		lh4f += y[i];
	}
	for (int k = K; k < remain_vertices.size(); k++)
	{
		u[k].set(GRB_DoubleAttr_LB, 0);
		lh4c -= instance.gamma * k * (k - 1) * u[k];
		lh4f -= k * u[k];
		lh4g += u[k];
	}
	sub_model.addConstr(lh4b, GRB_GREATER_EQUAL, instance.desired_clique_weight + 0.001);
	sub_model.addConstr(lh4c, GRB_GREATER_EQUAL, 0);
	sub_model.addConstr(lh4f, GRB_EQUAL, 0);
	sub_model.addConstr(lh4g, GRB_EQUAL, 1);
	sub_model.set(GRB_IntParam_Threads, 1);
	clock_t build_end = clock();

	double passing_time = static_cast<double>(build_end - ls.start) / CLOCKS_PER_SEC;
	if (passing_time >= ls.total_time_limit)
	{
		return true;
	}
	double optimization_time_limit = ls.total_time_limit - passing_time;
	sub_model.set(GRB_DoubleParam_TimeLimit, optimization_time_limit);
	sub_model.update();

	sub_model.optimize();

	if (sub_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
	{
		for (int i = 0; i < remain_vertices.size(); i++)
		{
			double temp_solution_i = y[i].get(GRB_DoubleAttr_X);
			if (temp_solution_i > 0.5)
			{
				gamma_clique.push_back(remain_vertices[i]);
			}
		}
		delete[] y;
		delete[] u;
		delete sub_env;
		return false;
	}
	else if (sub_model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
	{
		delete[] y;
		delete[] u;
		delete sub_env;
		return true;
	}
	else
	{
		delete[] y;
		delete[] u;
		delete sub_env;
		return false;
	}
}

vector<int> get_A(WGCI_instance& instance, vector<int> vertices, vector<int> candidate)
{
	double weight = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		weight += instance.vertices[vertices[i]].vertex_weight;
	}
	vector<int> temp_a;
	for (int i = 0; i < candidate.size(); i++)
	{
		if (weight - instance.vertices[candidate[i]].vertex_weight >= instance.desired_clique_weight)
		{
			temp_a.push_back(candidate[i]);
		}
	}
	return temp_a;
}

int get_BIA(WGCI_instance& instance, vector<int> vertices, vector<bool>& remain_edges, vector<int> a)
{//if find a vertex in B intersect with a, then return position of the vertex in vertices, else return -1
	int potential_edge = 0;
	int Nedge = 0;
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			int Eid = instance.two_d[vertices[i]][vertices[j]];
			if (Eid != -1 && remain_edges[Eid])
			{
				Nedge++;
			}
		}
	}
	double Desired_edge = (vertices.size() - 1) * (vertices.size() - 2) / 2 * instance.gamma;
	for (int i = 0; i < vertices.size(); i++)
	{
		int degree = 0;
		for (int j = 0; j < vertices.size(); j++)
		{
			if (j == i)
			{
				continue;
			}
			int Eid = instance.two_d[vertices[i]][vertices[j]];
			if (Eid != -1 && remain_edges[Eid])
			{
				degree++;
			}
		}
		if ((Nedge - degree) > Desired_edge)
		{
			for (int j = 0; j < a.size(); j++)
			{
				if (a[j] == vertices[i])
				{
					return i;//position of the interesect vertex in vertices
				}
			}

		}
	}
	return -1;
}

gamma_clique get_EGC(WGCI_instance& instance, vector<int> vertices, vector<bool>& remain_edges)
{
	vector<int> A = get_A(instance, vertices, vertices);
	if (A.size() != 0)
	{
		int inter_pos = get_BIA(instance, vertices, remain_edges, A);
		while (inter_pos != -1)
		{
			vertices[inter_pos] = vertices.back();
			vertices.pop_back();
			A = get_A(instance, vertices, A);
			inter_pos = get_BIA(instance, vertices, remain_edges, A);
		}
	}
	vector<int> edges = get_edges(instance, vertices, remain_edges);
	return gamma_clique(vertices, edges);
}

vector<int> get_edges(WGCI_instance& instance, vector<int>& vertices, vector<bool>& remain_edges)
{//find edges in the eta set
	double desired_edge = vertices.size() * (vertices.size() - 1) / 2 * instance.gamma;
	int edge_number = ceil(desired_edge);
	int count_e = 0;
	vector<int> edges;
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			int temp_id = instance.two_d[vertices[i]][vertices[j]];
			if (temp_id != -1 && remain_edges[temp_id])
			{//exist edge
				edges.push_back(instance.one_d[temp_id].ID);
				count_e++;
				if (count_e == edge_number)
				{
					return edges;
				}
			}
		}
	}
}

vector<int> get_candidate_set(WGCI_instance& instance, vector<int>& candidate_set, int target)
{
	vector<int> temp_candidate;
	for (int i = target + 1; i < candidate_set.size(); i++)
	{
		int id = instance.two_d[candidate_set[target]][candidate_set[i]];
		if (id > -0.5)
		{
			temp_candidate.push_back(candidate_set[i]);
		}
	}
	return temp_candidate;
}

int readdata(WGCI_instance& instance, string filename)
{
	ifstream fin;
	fin.open(filename);
	string str3;
	int number_of_points, number_of_edges;
	fin >> number_of_points;
	fin >> number_of_edges;
	fin >> instance.desired_clique_weight;
	fin >> instance.gamma;
	instance.two_d.resize(number_of_points);
	instance.vertices.resize(number_of_points);
	for (int i = 0; i < number_of_points; i++)
	{
		instance.two_d[i].resize(number_of_points, -1);
	}
	getline(fin, str3);//finish current line
	int one_code = 0;//code ID of edges
	for (int i = 0; i < number_of_points; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		ss0 >> instance.vertices[i].vertex_weight;
		ss0 >> instance.vertices[i].blocking_cost;
	}
	for (int i = 0; i < number_of_edges; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		int h, t;
		ss0 >> h;
		ss0 >> t;
		h--;
		t--;
		if (instance.two_d[h][t] == -1)
		{
			double temp_cost;
			ss0 >> temp_cost;
			instance.one_d.push_back(edge(h, t, one_code, temp_cost));
			instance.two_d[h][t] = one_code++;
			instance.two_d[t][h] = instance.two_d[h][t];
		}
	}
	return number_of_edges;
}