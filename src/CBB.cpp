#include <iostream>

#include <vector>
using namespace std;
#include <list>
#include <utility> 
#include <fstream>
#include <sstream>
#include <time.h>

#include <algorithm>
class element
{//an element of a vertex or an edge
public:
	int e;
	double blocking_cost;
	double vertex_weight;

	element(int id, double bc)
	{
		blocking_cost = bc;
		vertex_weight = 0;
		e = id;
	}
	element(int id, double bc, double vw)
	{
		blocking_cost = bc;
		vertex_weight = vw;
		e = -1;
	}
};

class edge
{
public:
	int head;
	int tail;
	int IDm;//index in mix of instance
	edge(int hd, int tl, int id)
	{
		head = hd;
		tail = tl;
		IDm = id;
	}

	edge()
	{
		head = 0;
		tail = 0;
		IDm = 0;
	}
};

class WKCI_instance
{
public:
	int vertices_size;//idv is the same as idm
	vector <vector<int>> two_d;//store edge id correspond to one_d
	vector<vector<int>> neighbor;//store the neighbor of each vertex
	//==-1 means no edge
	vector <edge> one_d;
	vector <element> mix;//store all vertices and edges
	double desired_clique_weight;
	double gamma;
	double total_edge_costs;

	WKCI_instance()
	{
		vertices_size = 0;
		two_d = vector <vector<int>>();
		one_d = vector<edge>();
		desired_clique_weight = -1;
		gamma = 0;
		total_edge_costs = 0;
	}
};

class node_gamma
{
public:
	vector<int> solution;//format of mix
	vector<int> candidate; //format of mix
	double weight;
	node_gamma(vector<int>& bs, vector<int>& cs, double w)
	{
		solution = bs;
		candidate = cs;
		weight = w;
	}
	node_gamma()
	{
		weight = 0;
	}
};

void a_sort_ratio(vector<pair<int, double>>& candidate);

class node
{
public:
	list<int> kv;//union of vertices in not blocking and candidate to calculate K
	vector<int> status;//0 block,1 not block,2 candidate,3 not block edge with endpoint blocked, format of mix
	list<pair<int, double>> candidate;//idm,theta
	int h;
	double solution_cost;
	bool b;//true if it's associated with block h
	int count_edge;
	double desired_edge;//number of desired edge based on updated K
	vector<int> degrees;
	int count_v;//number of vertices with degree not less than degree bound
	double degree_bound;
	int K;
	double LB;
	//Parametrized Constructor
	node(list<int>& nb, vector<int>& ub, list<pair<int, double>>& cc, int vh, double sc, bool vb, int ce, vector<int> dd, int cv, double de, double db, int k, double lb)
	{
		kv = nb;
		status = ub;
		candidate = cc;
		h = vh;
		solution_cost = sc;
		b = vb;
		count_edge = ce;
		degrees = dd;
		count_v = cv;
		desired_edge = de;
		degree_bound = db;
		K = k;
		LB = lb;
	}
	node(WKCI_instance& instance)
	{//create root node
		vector<pair<int, double>>tkv(instance.vertices_size);
		for (int i = 0; i < instance.vertices_size; i++)
		{
			tkv[i] = make_pair(i, instance.mix[i].vertex_weight);
		}
		a_sort_ratio(tkv);

		while (tkv.size() > 0)
		{
			kv.push_back(tkv.back().first);
			tkv.pop_back();
		}
		for (int i = 0; i < instance.vertices_size; i++)
		{
			double theta = instance.neighbor[i].size();
			candidate.push_back(make_pair(i, theta / instance.mix[i].blocking_cost));
		}
		for (int i = instance.vertices_size; i < instance.mix.size(); i++)
		{
			candidate.push_back(make_pair(i, 1 / instance.mix[i].blocking_cost));
		}
		status.resize(instance.mix.size(), 2);
		h = -1;
		solution_cost = 0;
		count_edge = instance.one_d.size();
		desired_edge = -1;
		degree_bound = 0;
		count_v = 0;
		degrees.resize(instance.vertices_size);
		for (int i = 0; i < instance.vertices_size; i++)
		{
			degrees[i] = instance.neighbor[i].size();
			if (instance.neighbor[i].size() > 0)
			{
				count_v++;
			}
		}
		K = 0;
	}
	node()
	{
		h = -1;
		solution_cost = 0;
		desired_edge = -1;
		degree_bound = -1;
		K = 0;
	}
};

class Hvertex
{//used in lower bound to indicate a vertex
public:
	int id;//id of mix in instance
	double score;//theta*weight
	bool isud;//true if it's in the undecided set
	//Parametrized Constructor
	Hvertex(int h, double sc, bool vb)
	{
		id = h;
		score = sc;
		isud = vb;
	}
	Hvertex()
	{
		id = -1;
		score = 0;
	}
};

class EGC
{//store an EGC
public:
	vector<int> EGC1;//id of mix in instance
	int cost;//cost of minimum cost element
	vector<bool> EGC2;//true if element is in the EGC
	int mh;//minimum cost element
	//Parametrized Constructor
	EGC(vector<int> egc1, int sc, vector<bool> egc2, int h)
	{
		mh = h;
		cost = sc;
		EGC1 = egc1;
		EGC2 = egc2;
	}
	EGC()
	{
		mh = -1;
		cost = 0;
	}
};

vector<int> find_gamma_clique(WKCI_instance& instance, list<int>& vertices, vector<int>& Ghat, int hq, clock_t start,double);

int readdata(WKCI_instance& instance, string filename);

void runmodel(string);

int c_branch(WKCI_instance& instance, int t, vector<node_gamma>& tree);

bool c_check_inf(WKCI_instance& instance, vector<int>& vertices, vector<int>& Ghat, node_gamma target);

int branch(WKCI_instance& instance, int t, vector<node>& tree, int& node_count);

int HTbranch(WKCI_instance& instance, int q, vector<node>& tree, int& node_count);
bool check_infeasibility(WKCI_instance& instance, node& tq, vector<EGC>& egclist, clock_t start,double);

int get_K(WKCI_instance& instance, list<int> weight_v, int lK, int pos);

vector<bool> get_EGC_nb(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat);

vector<bool> get_EGC_un(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat);

vector<int> get_A(WKCI_instance& instance, vector<int> vertices, vector<int> candidate);

int get_BIA_nb(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat, int& Nedge, vector<int>& a);
int get_BIA_un(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat, int& Nedge, vector<int>& a);

double rho(WKCI_instance& instance, vector<int>& vertices, vector<int>& Ghat, node_gamma& target, int D);

double a(WKCI_instance& instance, vector<int>& vertices, vector<int>& Ghat, node_gamma& target, int v, bool isU, int D);
//update edge
void update_b_edges(WKCI_instance& instance, node& tp, int v, double degree_bound);

void update_edge_theta(WKCI_instance& instance, node& tp, int e);
//remove blocked edge incident to the vertex from the block set and return blocking cost of edges
double update_beb(WKCI_instance& instance, node& tp, int v);

int get_edge_nb(WKCI_instance& instance, vector<int>& candidates, int v);

int get_edge_un(WKCI_instance& instance, vector<int>& candidates, int v);

int find_EGC_union(vector<EGC>& egclist, vector<int>& mix, int start);

bool find_EGC_nb(vector<EGC>& egclist, vector<int>& nb, int hq);

bool check_gamma_nb(WKCI_instance& instance, vector<int>& Ghat, vector<int>& vertices);

bool check_gamma_un(WKCI_instance& instance, vector<int>& Ghat, vector<int>& vertices);

double lower_bound(WKCI_instance& instance, node q, vector<EGC>& egclist, clock_t start);

int insert_EGC(WKCI_instance& instance, vector<EGC>& egclist, vector<bool>& egc);

double LB_unprocessed_nodes(WKCI_instance& instance, int cn, vector<node>& tree);

int get_max(vector<pair<int, double>> candidate);

bool d_v_w(const Hvertex& a, const Hvertex& b);

void d_sort_lb(vector<node>& candidate);

int main()
{
	ofstream fout;
	fout.open("CBB.csv");
	fout << "File Name,|V|,density,Eta,Gamma,number_of_EGC,Time Limit (s),Run Time (s),Gap,Best Sol, Best Bound, Nodes #,processed_node,fathomed by bound,fathomed by infeasibility, fathomed by feasibility,update solution,sols,removed vetices,removed edges" << endl;
	fout.close();

	ifstream fin;
	fin.open("filename.txt");
	string str3;
	int number_of_files;
	fin >> number_of_files;
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
	double total_time_limit = 3600.0;
	WKCI_instance instance;
	readdata(instance, filename);

	cout << filename << " start" << endl;
	clock_t pre_process_start = clock();
	vector<node> tree;
	tree.push_back(node(instance));

	int K_all = get_K(instance, tree[0].kv, 1, 0);
	double degree_bound = instance.gamma * K_all;
	double Desired_edge = static_cast<double>((K_all) * (K_all - 1)) / 2 * instance.gamma;
	tree[0].desired_edge = Desired_edge;
	tree[0].degree_bound = degree_bound;
	tree[0].K = K_all;
	vector<EGC> EGC_list;
	tree[0].LB = 0;
	vector<int> incumbent_solution;
	incumbent_solution.resize(instance.vertices_size, 1);
	incumbent_solution.resize(instance.mix.size(), 0);
	double incumbent_cost = instance.total_edge_costs;
	clock_t now = clock();
	cout << filename << ' ' << static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC << "s, UB done." << endl;

	int current_node = 0;
	int node_count = 1;//the root node
	int processed_node = 0;//the root node
	int fathomed_by_bound = 0;
	int fathomed_by_infeasibility = 0;
	int fathomed_by_feasibility = 0;
	int update_solution = 0;

	double rest_LB = 0;
	bool solved = true;
	bool find_s = false;
	while (current_node != -1)
	{
		//check time limit
		clock_t now = clock();
		double duration = static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC;
		if ((find_s && duration > total_time_limit / 2)||duration>total_time_limit)
		{
			cout << filename << ' ' << static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC << "s, half time, start to close gap." << endl;
			break;
		}
		bool infeasibility = false;
		if (!tree[current_node].b)
		{//associated with not blocking
			if (tree[current_node].h < instance.vertices_size)
			{
				infeasibility = check_infeasibility(instance, tree[current_node], EGC_list, pre_process_start, total_time_limit);
			}
			else
			{//edge
				edge tempe = instance.one_d[instance.mix[tree[current_node].h].e];
				if (tree[current_node].status[tempe.head] == 1 && tree[current_node].status[tempe.tail] == 1)
				{
					infeasibility = check_infeasibility(instance, tree[current_node], EGC_list, pre_process_start, total_time_limit);
				}
			}
		}
		if (infeasibility)
		{
			fathomed_by_infeasibility++;
			current_node--;
		}
		else if (tree[current_node].count_v == 0 || tree[current_node].count_edge < tree[current_node].desired_edge || tree[current_node].candidate.size() == 0)
		{//fathom by feasibility
			fathomed_by_feasibility++;
			if (tree[current_node].solution_cost < incumbent_cost)
			{
				find_s = true;
				incumbent_cost = tree[current_node].solution_cost;
				incumbent_solution = tree[current_node].status;
				update_solution++;
			}
			current_node--;
		}
		else
		{
			double newLB = lower_bound(instance, tree[current_node], EGC_list, pre_process_start);
			if (newLB < incumbent_cost)
			{
				tree[current_node].LB = newLB;
				current_node = branch(instance, current_node, tree, node_count);
			}
			else//fathom by bound
			{
				fathomed_by_bound++;
				current_node--;
			}
		}
		processed_node++;
	}
	tree.resize(current_node + 1);
	d_sort_lb(tree);

	while (current_node != -1)
	{
		//check time limit
		clock_t now = clock();
		if (static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC > total_time_limit)
		{
			cout << filename << ' ' << static_cast<double>(now - pre_process_start) / CLOCKS_PER_SEC << "s, TL, start to calculate gap." << endl;
			//calculate lower bound for all unprocessed nodes
			rest_LB = LB_unprocessed_nodes(instance, current_node + 1, tree);

			solved = false;
			break;
		}

		bool infeasibility = false;
		if (!tree[current_node].b)
		{//associated with not blocking
			if (tree[current_node].h < instance.vertices_size)
			{
				infeasibility = check_infeasibility(instance, tree[current_node], EGC_list, pre_process_start, total_time_limit);
			}
			else
			{//edge
				edge tempe = instance.one_d[instance.mix[tree[current_node].h].e];
				if (tree[current_node].status[tempe.head] == 1 && tree[current_node].status[tempe.tail] == 1)
				{
					infeasibility = check_infeasibility(instance, tree[current_node], EGC_list, pre_process_start, total_time_limit);
				}
			}
		}
		if (infeasibility)
		{
			fathomed_by_infeasibility++;
			current_node--;
		}
		else if (tree[current_node].count_v == 0 || tree[current_node].count_edge < tree[current_node].desired_edge || tree[current_node].candidate.size() == 0)
		{//fathom by feasibility
			fathomed_by_feasibility++;
			if (tree[current_node].solution_cost < incumbent_cost)
			{
				incumbent_cost = tree[current_node].solution_cost;
				incumbent_solution = tree[current_node].status;
				update_solution++;
			}
			current_node--;
		}
		else
		{
			double newLB = lower_bound(instance, tree[current_node], EGC_list, pre_process_start);
			if (newLB < incumbent_cost)
			{
				if (newLB > tree[current_node].LB)
				{
					tree[current_node].LB = newLB;
				}
				current_node = HTbranch(instance, current_node, tree, node_count);
			}
			else//fathom by bound
			{
				fathomed_by_bound++;
				current_node--;
			}
		}
		processed_node++;
	}
	clock_t end = clock();
	cout << filename << " end" << endl;
	ofstream fout;
	fout.open("CBB.csv", std::ios_base::app);
	fout << filename << ',' << instance.vertices_size;
	fout << ',' << static_cast<double>(instance.one_d.size()) / static_cast<double>((instance.vertices_size * (instance.vertices_size - 1)) / 2);
	fout << ',' << instance.desired_clique_weight << ',' << instance.gamma;
	fout << ',' << EGC_list.size();
	fout << ',' << total_time_limit << ",";

	if (solved)
	{
		fout << static_cast<double>(end - pre_process_start) / CLOCKS_PER_SEC;
		fout << ",0,";
		fout << incumbent_cost << ","<< incumbent_cost<<',';
	}
	else
	{
		fout << "TL" << ',';
		fout << static_cast<double> (incumbent_cost - rest_LB) / (incumbent_cost) << ',';
		fout << incumbent_cost << "," << rest_LB << ",";
	}
	fout << node_count << ',';
	fout << processed_node << ',';
	fout << fathomed_by_bound << ',';
	fout << fathomed_by_infeasibility << ',';
	fout << fathomed_by_feasibility << ',';
	fout << update_solution << ',';
	int countv = 0;
	int counte = 0;
	for (int i = 0; i < instance.vertices_size; i++)
	{
		if (incumbent_solution[i]==0)
		{
			fout << i + 1 << ';';
			countv++;
		}
	}
	for (int i = instance.vertices_size; i < incumbent_solution.size(); i++)
	{
		if (incumbent_solution[i] == 0)
		{
			fout << '(' << instance.one_d[instance.mix[i].e].head + 1 << ';' << instance.one_d[instance.mix[i].e].tail + 1 << ')';
			counte++;
		}
	}
	fout << ',' << countv << ',' << counte;
	fout << endl;
	fout.close();
}

vector<int> find_gamma_clique(WKCI_instance& instance, list<int>& kv, vector<int>& Ghat, int hq, clock_t start, double ttl)
{//return vertices in format of mix,return {-1} if reach time limit
	vector<node_gamma> tree;
	//initialize the root node
	vector<pair<int, double>> temp_vertices;
	vector<int> vertices;
	for (list<int>::iterator it = kv.begin(); it != kv.end(); it++)
	{
		if (Ghat[*it] == 1)
		{//vertex is in not blocking set
			vertices.push_back(*it);
			int ti = *it;
			double temp = instance.mix[ti].vertex_weight;
			temp_vertices.push_back(make_pair(ti, temp));
		}
	}
	a_sort_ratio(temp_vertices);
	tree.push_back(node_gamma());
	if (hq < instance.vertices_size)
	{
		for (int i = 0; i < temp_vertices.size(); i++)
		{//initialize candidate set of root node
			if (temp_vertices[i].first != hq)
			{
				tree[0].candidate.push_back(temp_vertices[i].first);
			}
			else
			{
				tree[0].solution.push_back(hq);
				tree[0].weight += instance.mix[hq].vertex_weight;

			}
		}
	}
	else
	{
		edge temp_e = instance.one_d[instance.mix[hq].e];
		for (int i = 0; i < temp_vertices.size(); i++)
		{//initialize candidate set of root node
			int tempv = temp_vertices[i].first;
			if (tempv != temp_e.head && tempv != temp_e.tail)
			{
				tree[0].candidate.push_back(temp_vertices[i].first);
			}
			else
			{
				tree[0].solution.push_back(temp_vertices[i].first);
				tree[0].weight += instance.mix[temp_vertices[i].first].vertex_weight;

			}
		}
	}
	//search tree
	int current_node = 0;
	while (current_node != -1)
	{
		clock_t now = clock();
		if (static_cast<double>(now - start) / CLOCKS_PER_SEC > ttl)
		{
			return { -1 };
		}
		if (tree[current_node].weight > instance.desired_clique_weight && check_gamma_nb(instance, Ghat, tree[current_node].solution))
		{
			return tree[current_node].solution;
		}
		else
		{
			if (c_check_inf(instance, vertices, Ghat, tree[current_node]) == false)
			{
				current_node = c_branch(instance, current_node, tree);
			}
			else
			{
				current_node--;
			}
		}
	}
	return vector<int>();
}

int c_branch(WKCI_instance& instance, int t, vector<node_gamma>& tree)
{
	vector<node_gamma> temp_nodes;
	while (tree[t].candidate.size() > 0)
	{
		double temp_weight = tree[t].weight + instance.mix[tree[t].candidate.back()].vertex_weight;
		vector<int> temp_solution_set = tree[t].solution;
		temp_solution_set.push_back(tree[t].candidate.back());
		tree[t].candidate.pop_back();//remove the last one
		vector<int> temp_candidate_set = tree[t].candidate;
		temp_nodes.push_back(node_gamma(temp_solution_set, temp_candidate_set, temp_weight));
	}
	if (temp_nodes.size() == 1)
	{
		tree[t] = temp_nodes[0];
		return t;
	}
	else
	{
		int i = temp_nodes.size();
		tree.resize(t + i);
		for (int j = 0; j < temp_nodes.size(); j++)
		{
			tree[t + j] = temp_nodes[temp_nodes.size() - j - 1];
		}
		return t + temp_nodes.size() - 1;
	}
}

bool check_gamma_nb(WKCI_instance& instance, vector<int>& Ghat, vector<int>& vertices)
{//vertices in format of mix in instance
	int count = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			if (instance.two_d[vertices[i]][vertices[j]] > -1)
			{
				if (Ghat[instance.one_d[instance.two_d[vertices[i]][vertices[j]]].IDm] == 1)
				{
					count++;
				}
			}
		}
	}
	double density = static_cast<double>(count) / (static_cast<double>(vertices.size() * (vertices.size() - 1)) / 2);
	if (density >= instance.gamma)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool check_gamma_un(WKCI_instance& instance, vector<int>& Ghat, vector<int>& vertices)
{//vertices in format of mix in instance
	int count = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			if (instance.two_d[vertices[i]][vertices[j]] > -1)
			{
				if (Ghat[instance.one_d[instance.two_d[vertices[i]][vertices[j]]].IDm] > 0)
				{
					count++;
				}
			}
		}
	}
	double density = static_cast<double>(count) / (static_cast<double>(vertices.size() * (vertices.size() - 1)) / 2);
	if (density >= instance.gamma)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool c_check_inf(WKCI_instance& instance, vector<int>& vertices, vector<int>& Ghat, node_gamma target)
{
	int u = target.solution.size() + target.candidate.size();
	while (u > target.solution.size() && rho(instance, vertices, Ghat, target, u) < instance.gamma)
	{
		u--;
	}
	double wbar = 0;
	for (int i = 0; i < u - target.solution.size(); i++)
	{
		wbar += instance.mix[target.candidate[target.candidate.size() - i - 1]].vertex_weight;
	}
	if (u == target.solution.size() || (target.weight + wbar) <= instance.desired_clique_weight)
	{
		return true;
	}
	else
	{
		return false;
	}
}

double rho(WKCI_instance& instance, vector<int>& vertices, vector<int>& Ghat, node_gamma& target, int D)
{
	double front = 1 / static_cast<double>(D * (D - 1));
	double back = 0;
	for (int i = 0; i < target.solution.size(); i++)
	{
		back += a(instance, vertices, Ghat, target, target.solution[i], false, D);
	}
	vector<double> temp_vertices(target.candidate.size());
	for (int i = 0; i < target.candidate.size(); i++)
	{
		temp_vertices[i] = a(instance, vertices, Ghat, target, target.candidate[i], true, D);
	}
	sort(temp_vertices.begin(), temp_vertices.end());
	for (int i = 0; i < D - target.solution.size(); i++)
	{
		back += temp_vertices[temp_vertices.size() - i - 1];
	}
	return front * back;
}

double a(WKCI_instance& instance, vector<int>& vertices, vector<int>& Ghat, node_gamma& target, int v, bool isU, int D)
{//v is in format of vertices in instance
	int NGi = 0;
	int NGiS = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		int ti = vertices[i];
		if (instance.two_d[ti][v] > -1 && Ghat[instance.one_d[instance.two_d[ti][v]].IDm] == 1)
		{
			NGi++;
			for (int j = 0; j < target.solution.size(); j++)
			{
				if (target.solution[j] == ti)
				{
					NGiS++;
				}
			}
		}

	}
	if (isU)
	{
		return NGiS + min(NGi - NGiS, static_cast<int>(D - target.solution.size() - 1));
	}
	else
	{
		return NGiS + min(NGi - NGiS, static_cast<int>(D - target.solution.size()));
	}
}

int readdata(WKCI_instance& instance, string filename)
{//idv is the same as idm
	ifstream fin;
	fin.open(filename);
	string str3;
	int number_of_points, number_of_edges;
	fin >> number_of_points;
	fin >> number_of_edges;
	fin >> instance.desired_clique_weight;
	fin >> instance.gamma;
	instance.total_edge_costs = 0;
	instance.two_d.resize(number_of_points);
	instance.neighbor.resize(number_of_points);
	instance.vertices_size = number_of_points;
	for (int i = 0; i < number_of_points; i++)
	{
		instance.two_d[i].resize(number_of_points, -1);
	}
	getline(fin, str3);//finish current line
	for (int i = 0; i < number_of_points; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		double temp_w;
		ss0 >> temp_w;
		double bc;
		ss0 >> bc;
		instance.mix.push_back(element(i, bc, temp_w));
	}
	int one_code = 0, m_code = instance.mix.size();//code ID of edges and mix
	for (int i = 0; i < number_of_edges; i++)
	{
		getline(fin, str3);
		istringstream ss0(str3);
		int head;
		ss0 >> head;
		head--;
		int tail;
		ss0 >> tail;
		tail--;
		double temp_cost;
		ss0 >> temp_cost;
		if (head!=tail&&instance.two_d[head][tail] == -1)
		{
			instance.total_edge_costs += temp_cost;
			instance.one_d.push_back(edge(head, tail, m_code));
			instance.two_d[head][tail] = one_code;
			instance.two_d[tail][head] = one_code;
			instance.mix.push_back(element(one_code, temp_cost));
			m_code++;
			instance.neighbor[head].push_back(tail);
			instance.neighbor[tail].push_back(head);
			one_code++;
		}
	}
	return number_of_edges;
}

int branch(WKCI_instance& instance, int t, vector<node>& tree, int& node_count)
{//return position of tree node
	list<pair<int, double>>::iterator ith = tree[t].candidate.begin();
	for (list<pair<int, double>>::iterator it0 = tree[t].candidate.begin(); it0 != tree[t].candidate.end(); ++it0)
	{
		if ((*it0).second > (*ith).second)
		{
			ith = it0;
		}
	}

	int h = (*ith).first;
	tree[t].candidate.erase(ith);
	//create tp associated with blocking h
	double temp_solution_cost = tree[t].solution_cost + instance.mix[h].blocking_cost;
	vector<int> temp_s = tree[t].status;
	temp_s[h] = 0;
	list<pair<int, double>> temp_c = tree[t].candidate;
	int temp_ec = tree[t].count_edge;
	vector<int> temp_dd = tree[t].degrees;
	int temp_cv = tree[t].count_v;
	list<int> temp_kv = tree[t].kv;
	double tempde = tree[t].desired_edge;
	int temp_k = tree[t].K;
	node b_node = node(temp_kv, temp_s, temp_c, h, temp_solution_cost, true, temp_ec, temp_dd, temp_cv, tempde, tree[t].degree_bound, temp_k, tree[t].LB);

	if (h < instance.vertices_size)
	{
		int pos = 0;
		for (list<int>::iterator it1 = b_node.kv.begin(); it1 != b_node.kv.end(); ++it1)
		{
			if (*it1 == h)
			{
				b_node.K = get_K(instance, b_node.kv, b_node.K, pos);
				b_node.degree_bound = instance.gamma * b_node.K;
				b_node.kv.erase(it1);
				break;
			}
			pos++;
		}

		if (temp_dd[h] >= b_node.degree_bound)
		{
			b_node.count_v--;
		}
		update_b_edges(instance, b_node, h, b_node.degree_bound);
		b_node.solution_cost -= update_beb(instance, b_node, h);
		double Desired_edge = static_cast<double>((b_node.K) * (b_node.K - 1)) / 2 * instance.gamma;
		b_node.desired_edge = Desired_edge;
	}
	else
	{
		b_node.count_edge--;
		//update theta and countv in update_edge_theta
		update_edge_theta(instance, b_node, instance.mix[h].e);
	}
	//create tq associated with not blocking h
	temp_s = tree[t].status;
	temp_s[h] = 1;
	temp_kv = tree[t].kv;
	//temp_c is also tree[t].candidate;
	temp_c = tree[t].candidate;
	temp_solution_cost = tree[t].solution_cost;
	temp_ec = tree[t].count_edge;
	temp_dd = tree[t].degrees;
	temp_cv = tree[t].count_v;
	node nb_node = node(temp_kv, temp_s, temp_c, h, temp_solution_cost, false, temp_ec, temp_dd, temp_cv, tempde, tree[t].degree_bound, temp_k, tree[t].LB);
	//node [t+1] becomes node block ,[t] is not block
	node_count += 2;//block and not block
	tree.resize(t + 2);
	tree[t] = nb_node;
	tree[t + 1] = b_node;
	return t + 1;
}
int HTbranch(WKCI_instance& instance, int q, vector<node>& tree, int& node_count)
{//return position of tree node
	list<pair<int, double>>::iterator ith = tree[q].candidate.begin();
	for (list<pair<int, double>>::iterator it0 = tree[q].candidate.begin(); it0 != tree[q].candidate.end(); ++it0)
	{
		if ((*it0).second > (*ith).second)
		{
			ith = it0;
		}
	}

	int h = (*ith).first;
	tree[q].candidate.erase(ith);
	//create tp associated with blocking h
	double temp_solution_cost = tree[q].solution_cost + instance.mix[h].blocking_cost;
	vector<int> temp_s = tree[q].status;
	temp_s[h] = 0;
	list<pair<int, double>> temp_c = tree[q].candidate;
	int temp_ec = tree[q].count_edge;
	vector<int> temp_dd = tree[q].degrees;
	int temp_cv = tree[q].count_v;
	list<int> temp_kv = tree[q].kv;
	double tempde = tree[q].desired_edge;
	int temp_k = tree[q].K;
	node b_node = node(temp_kv, temp_s, temp_c, h, temp_solution_cost, true, temp_ec, temp_dd, temp_cv, tempde, tree[q].degree_bound, temp_k, tree[q].LB);

	if (h < instance.vertices_size)
	{
		int pos = 0;
		for (list<int>::iterator it1 = b_node.kv.begin(); it1 != b_node.kv.end(); ++it1)
		{
			if (*it1 == h)
			{
				b_node.K = get_K(instance, b_node.kv, b_node.K, pos);
				b_node.degree_bound = instance.gamma * b_node.K;
				b_node.kv.erase(it1);
				break;
			}
			pos++;
		}

		if (temp_dd[h] >= b_node.degree_bound)
		{
			b_node.count_v--;
		}
		update_b_edges(instance, b_node, h, b_node.degree_bound);
		b_node.solution_cost -= update_beb(instance, b_node, h);
		double Desired_edge = static_cast<double>((b_node.K) * (b_node.K - 1)) / 2 * instance.gamma;
		b_node.desired_edge = Desired_edge;
	}
	else
	{
		b_node.count_edge--;
		//update theta and countv in update_edge_theta
		update_edge_theta(instance, b_node, instance.mix[h].e);
	}
	//create tq associated with not blocking h
	temp_s = tree[q].status;
	temp_s[h] = 1;
	temp_kv = tree[q].kv;
	//temp_c is also tree[t].candidate;
	temp_c = tree[q].candidate;
	temp_solution_cost = tree[q].solution_cost;
	temp_ec = tree[q].count_edge;
	temp_dd = tree[q].degrees;
	temp_cv = tree[q].count_v;
	node nb_node = node(temp_kv, temp_s, temp_c, h, temp_solution_cost, false, temp_ec, temp_dd, temp_cv, tempde, tree[q].degree_bound, temp_k, tree[q].LB);
	node_count += 2;//block and not block
	tree.resize(tree.size() + 1);
	for (int i = q - 1; i >= 0; i--)
	{
		if (i < 0)
		{
			tree[0] = b_node;
			tree[1] = nb_node;
			return 1;
		}
		else if (tree[i].LB < b_node.LB)
		{
			tree[i + 2] = tree[i];
		}
		else
		{
			tree[i + 1] = b_node;
			tree[i + 2] = nb_node;
			return q + 1;
		}
	}
}

bool find_EGC_nb(vector<EGC>& egclist, vector<int>& nb, int hq)
{//check if there exist EGC in S+q containing hq, return true if find EGC
	for (int i = 0; i < egclist.size(); i++)
	{
		if (egclist[i].EGC2[hq])
		{
			bool contain = true;
			for (int j = 0; j < egclist[i].EGC1.size(); j++)
			{
				if (nb[egclist[i].EGC1[j]] != 1)
				{
					contain = false;
					break;
				}
			}
			if (contain)
			{
				return true;
			}
		}
	}
	return false;
}

bool check_infeasibility(WKCI_instance& instance, node& tq, vector<EGC>& egclist, clock_t start,double ttl)
{//return true if infeasible
	if (tq.h >= instance.vertices_size)
	{
		edge uhqe = instance.one_d[instance.mix[tq.h].e];
		for (list<pair<int, double>>::iterator it0 = tq.candidate.begin(); it0 != tq.candidate.end(); ++it0)
		{
			if ((*it0).first < instance.vertices_size)
			{
				if (uhqe.head == (*it0).first || uhqe.tail == (*it0).first)
				{
					return false;
				}
			}
		}
	}

	double temp_weight = 0;
	for (list<int>::iterator it = tq.kv.begin(); it != tq.kv.end(); it++)
	{
		if (tq.status[*it] > 0)
		{//vertex is in not blocking set
			temp_weight += instance.mix[*it].vertex_weight;
			if (temp_weight > instance.desired_clique_weight)
			{
				break;
			}
		}
	}
	if (temp_weight <= instance.desired_clique_weight)
	{
		return false;
	}

	vector<pair<int, double>> candidate;
	vector<double> weight_v;
	for (list<int>::iterator it = tq.kv.begin(); it != tq.kv.end(); it++)
	{
		if (tq.status[*it] > 0)
		{//vertex is in not blocking set
			int ti = *it;
			int theta = get_edge_nb(instance, tq.status, ti);
			candidate.push_back(make_pair(ti, theta));
		}
	}
	a_sort_ratio(candidate);
	int K = get_K(instance, tq.kv, 1, 0);

	double sum_theta = 0;
	for (int i = 0; i < K; i++)
	{
		sum_theta += candidate[candidate.size() - i - 1].second;
	}
	double temp = static_cast<double>(K * (K - 1));
	if (sum_theta / temp < instance.gamma)
	{
		return false;
	}

	if (find_EGC_nb(egclist, tq.status, tq.h))
	{
		return true;
	}

	vector<int> H = find_gamma_clique(instance, tq.kv, tq.status, tq.h, start,ttl);
	if (H.size() == 0)
	{
		return false;
	}
	else if (H.size() == 1&&H[0]==-1)
	{
		return true;
	}
	else
	{
		vector<bool> EGC2 = get_EGC_nb(instance, H, tq.status);
		insert_EGC(instance, egclist, EGC2);
		return true;
	}

}

int get_K(WKCI_instance& instance, list<int> weight_v, int lK, int pos)
{
	if (pos >= lK)
	{
		return lK;
	}
	int K = 0;
	double temp_weight = 0;
	for (list<int>::iterator it = weight_v.begin(); it != weight_v.end(); it++)
	{
		temp_weight += instance.mix[*it].vertex_weight;
		K++;
		if (temp_weight > instance.desired_clique_weight)
		{
			return K;
		}
	}
}

vector<bool> get_EGC_nb(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat)
{//vertices in format of mix
	int Nedge = 0;
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			int Eid = instance.two_d[vertices[i]][vertices[j]];
			if (Eid != -1 && Ghat[instance.one_d[Eid].IDm] == 1)
			{
				Nedge++;
			}
		}
	}
	vector<int> A = get_A(instance, vertices, vertices);
	int inter_pos = get_BIA_nb(instance, vertices, Ghat, Nedge, A);
	while (inter_pos != -1)
	{
		vertices[inter_pos] = vertices.back();
		vertices.pop_back();
		A = get_A(instance, vertices, A);
		inter_pos = get_BIA_nb(instance, vertices, Ghat, Nedge, A);
	}
	vector<bool> EGC;
	EGC.resize(instance.mix.size(), false);
	double desired_edge = static_cast<double>(vertices.size() * (vertices.size() - 1)) / 2 * instance.gamma;
	int edge_number = ceil(desired_edge);//need to get ceiling
	int temp_en = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		EGC[vertices[i]] = true;
	}
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			int Eid = instance.two_d[vertices[i]][vertices[j]];
			if (Eid != -1 && Ghat[instance.one_d[Eid].IDm] == 1)
			{
				temp_en++;
				EGC[instance.one_d[Eid].IDm] = true;
				if (temp_en == edge_number)
				{
					return EGC;
				}
			}
		}
	}
}

vector<bool> get_EGC_un(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat)
{//vertices in format of mix, e is edge that egc must contain
	int Nedge = 0;
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			int Eid = instance.two_d[vertices[i]][vertices[j]];
			if (Eid != -1 && Ghat[instance.one_d[Eid].IDm] > 0)
			{
				Nedge++;
			}
		}
	}
	vector<int> A = get_A(instance, vertices, vertices);
	int inter_pos = get_BIA_un(instance, vertices, Ghat, Nedge, A);
	while (inter_pos != -1)
	{
		vertices[inter_pos] = vertices.back();
		vertices.pop_back();
		A = get_A(instance, vertices, A);
		inter_pos = get_BIA_un(instance, vertices, Ghat, Nedge, A);
	}
	vector<bool> EGC;
	EGC.resize(instance.mix.size(), false);
	double desired_edge = static_cast<double>(vertices.size() * (vertices.size() - 1)) / 2 * instance.gamma;
	int edge_number = ceil(desired_edge);//need to get ceiling
	int temp_en = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		EGC[vertices[i]] = true;
	}
	for (int i = 0; i < vertices.size() - 1; i++)
	{
		for (int j = i + 1; j < vertices.size(); j++)
		{
			int Eid = instance.two_d[vertices[i]][vertices[j]];
			if (Eid != -1 && Ghat[instance.one_d[Eid].IDm] > 0)
			{
				temp_en++;
				EGC[instance.one_d[Eid].IDm] = true;
				if (temp_en == edge_number)
				{
					return EGC;
				}
			}
		}
	}
}

vector<int> get_A(WKCI_instance& instance, vector<int> vertices, vector<int> candidate)
{
	double weight = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		weight += instance.mix[vertices[i]].vertex_weight;
	}
	vector<int> temp_a;
	for (int i = 0; i < candidate.size(); i++)
	{
		if (weight - instance.mix[candidate[i]].vertex_weight >= instance.desired_clique_weight)
		{
			temp_a.push_back(candidate[i]);
		}
	}
	return temp_a;
}

int get_BIA_nb(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat, int& Nedge, vector<int>& a)
{//if find a vertex in B intersect with a, then return position of the vertex in vertices, else return -1
	int potential_edge = 0;
	double Desired_edge = static_cast<double>((vertices.size() - 1) * (vertices.size() - 2)) / 2 * instance.gamma;
	int position = 0;//store the position of a in vertices
	for (int i = 0; i < a.size(); i++)
	{
		int degree = 0;
		for (int j = 0; j < vertices.size(); j++)
		{
			if (a[i] == vertices[j])
			{
				position = j;
				continue;
			}
			int Eid = instance.two_d[a[i]][vertices[j]];
			if (Eid != -1 && Ghat[instance.one_d[Eid].IDm] == 1)
			{
				degree++;
			}
		}
		if ((Nedge - degree) > Desired_edge)
		{
			Nedge -= degree;
			a[i] = a.back();
			a.pop_back();

			return position;//position of the interesect vertex in vertices
		}
	}
	return -1;
}
int get_BIA_un(WKCI_instance& instance, vector<int> vertices, vector<int>& Ghat, int& Nedge, vector<int>& a)
{//if find a vertex in B intersect with a, then return position of the vertex in vertices, else return -1
	int potential_edge = 0;
	double Desired_edge = static_cast<double>((vertices.size() - 1) * (vertices.size() - 2)) / 2 * instance.gamma;
	int position = 0;//store the position of a in vertices
	for (int i = 0; i < a.size(); i++)
	{
		int degree = 0;
		for (int j = 0; j < vertices.size(); j++)
		{
			if (a[i] == vertices[j])
			{
				position = j;
				continue;
			}
			int Eid = instance.two_d[a[i]][vertices[j]];
			if (Eid != -1 && Ghat[instance.one_d[Eid].IDm] > 0)
			{
				degree++;
			}
		}
		if ((Nedge - degree) > Desired_edge)
		{
			Nedge -= degree;
			a[i] = a.back();
			a.pop_back();

			return position;//position of the interesect vertex in vertices
		}
	}
	return -1;
}

int get_max(vector<pair<int, double>> candidate)
{//return index of candidate
	int max = 0;
	double mRatio = 0;
	for (int i = 0; i < candidate.size(); i++)
	{
		if (candidate[i].second > mRatio)
		{//larger ratio
			max = i;
			mRatio = candidate[i].second;
		}
	}
	return max;
}

void update_b_edges(WKCI_instance& instance, node& tp, int v, double degree_bound)
{//remove incident edges of the blocked vertex from the candidate
	for (list<pair<int, double>>::iterator it = tp.candidate.begin(); it != tp.candidate.end(); )
	{
		if ((*it).first >= instance.vertices_size)
		{
			element temp = instance.mix[(*it).first];//edge
			if (instance.one_d[temp.e].head == v)
			{//remove the edge
				tp.count_edge--;
				int other_end = instance.one_d[temp.e].tail;
				for (list<pair<int, double>>::iterator it0 = tp.candidate.begin(); it0 != tp.candidate.end(); it0++)
				{
					if ((*it0).first < instance.vertices_size)
					{
						double temp1b = instance.mix[(*it0).first].blocking_cost;
						if ((*it0).first == other_end)
						{//find the other endpoint
							double temp_theta = (*it0).second * temp1b - 1;
							(*it0).second = temp_theta / temp1b;
							if (tp.degrees[other_end] >= degree_bound)
							{
								tp.degrees[other_end]--;
								if (tp.degrees[other_end] < degree_bound)
								{
									tp.count_v--;
								}
							}
							break;
						}
					}
				}
				tp.status[(*it).first] = 3;
				it = tp.candidate.erase(it);
			}
			else if (instance.one_d[temp.e].tail == v)
			{//remove the edge
				tp.count_edge--;
				int other_end = instance.one_d[temp.e].head;
				for (list<pair<int, double>>::iterator it0 = tp.candidate.begin(); it0 != tp.candidate.end(); it0++)
				{
					if ((*it0).first < instance.vertices_size)
					{
						element temp1 = instance.mix[(*it0).first];
						if ((*it0).first == other_end)
						{//find the other endpoint
							double temp_theta = (*it0).second * temp1.blocking_cost - 1;
							(*it0).second = temp_theta / temp1.blocking_cost;
							if (tp.degrees[other_end] >= degree_bound)
							{
								tp.degrees[other_end]--;
								if (tp.degrees[other_end] < degree_bound)
								{
									tp.count_v--;
								}
							}
							break;
						}
					}
				}
				tp.status[(*it).first] = 3;
				it = tp.candidate.erase(it);
			}
			else
			{
				it++;
			}
		}
		else
		{
			it++;
		}
	}
}

void update_edge_theta(WKCI_instance& instance, node& tp, int e)
{//update theta of the endpoints of the blocked edge in the candidate
	for (list<pair<int, double>>::iterator it = tp.candidate.begin(); it != tp.candidate.end(); it++)
	{
		int count = 0;
		if ((*it).first < instance.vertices_size)
		{
			double tempb = instance.mix[(*it).first].blocking_cost;
			if (instance.one_d[e].head == (*it).first || instance.one_d[e].tail == (*it).first)
			{//update theta
				double temp_theta = (*it).second * tempb - 1;
				(*it).second = temp_theta / tempb;
				//update degree 
				if (tp.degrees[(*it).first] >= tp.degree_bound)
				{
					tp.degrees[(*it).first]--;
					if (tp.degrees[(*it).first] < tp.degree_bound)
					{
						tp.count_v--;
					}
				}
				count++;
				if (count == 2)
				{
					return;
				}
			}
		}
	}
}

double update_beb(WKCI_instance& instance, node& tp, int v)
{//remove blocked edge incident to the vertex from the block set and return blocking cost of edges
	double bc = 0;
	for (int i = instance.vertices_size; i < tp.status.size(); i++)
	{//edge
		if (tp.status[i] == 0)
		{//blocked edge
			element e = instance.mix[i];//edge
			if (instance.one_d[e.e].head == v || instance.one_d[e.e].tail == v)
			{//remove the edge, set status of edge to 3
				bc += e.blocking_cost;
				tp.status[i] = 3;
			}
		}
	}
	return bc;
}

int get_edge_nb(WKCI_instance& instance, vector<int>& candidates, int v)
{//get the number of incident edges of the vertex from the not blocking set
	//v in format of vertices in instance
	int count = 0;
	for (int i = 0; i < instance.neighbor[v].size(); i++)
	{
		if (candidates[instance.one_d[instance.two_d[v][instance.neighbor[v][i]]].IDm] == 1)
		{
			count++;
		}
	}
	return count;
}

int get_edge_un(WKCI_instance& instance, vector<int>& candidates, int v)
{//get the number of incident edges of the vertex from the union of nb and U
	int count = 0;
	for (int i = 0; i < instance.neighbor[v].size(); i++)
	{
		if (candidates[instance.one_d[instance.two_d[v][instance.neighbor[v][i]]].IDm] > 0)
		{
			count++;
		}
	}
	return count;
}

int find_EGC_union(vector<EGC>& egclist, vector<int>& mix, int start)
{//check if there exist EGC in given vertices and edges, from start in the list
	for (int i = start; i < egclist.size(); i++)
	{
		bool contain = true;
		for (int j = 0; j < egclist[i].EGC1.size(); j++)
		{
			if (mix[egclist[i].EGC1[j]] == 0)
			{//element in EGC is not contained
				contain = false;
				break;
			}
		}
		if (contain)
		{
			return i;//contain this EGC
		}
	}
	return egclist.size();
}

double lower_bound(WKCI_instance& instance, node q, vector<EGC>& egclist, clock_t start)
{
	list<Hvertex> Hv;
	vector<int> Hmix;//element in H format of mix
	vector<int> He2;//id of edge, format of mix
	vector<bool> Heud;//1 if edge in ud, format of one_d in instance

	Hmix = q.status;
	Heud.resize(instance.one_d.size(), false);
	int HUq = 0;//number of elements in H that is in undecided set
	double HVweight = 0;//total weight of vertices in H
	for (int i = 0; i < instance.vertices_size; i++)
	{//vertex
		if (q.status[i] > 0)
		{
			double theta = get_edge_un(instance, Hmix, i);
			double temp = instance.mix[i].vertex_weight;
			if (q.status[i] == 1)
			{//not blocking vertex
				Hv.push_back(Hvertex(i, theta * temp, false));
			}
			else
			{//candidate vertex
				HUq++;
				Hv.push_back(Hvertex(i, theta * temp, true));
			}
			HVweight += temp;
		}
	}
	for (int i = instance.vertices_size; i < q.status.size(); i++)
	{//edge
		if (q.status[i] > 0)
		{
			if (q.status[i] == 2)
			{//candidate edge
				Heud[instance.mix[i].e] = true;
				HUq++;
			}
			He2.push_back(i);
		}
	}
	double c_t = 0;
	Hv.sort(d_v_w);
	int EGCid = find_EGC_union(egclist, Hmix, 0);
	while (HUq > 0 && HVweight > instance.desired_clique_weight)
	{
		if (EGCid < egclist.size())
		{
			bool egchud = false;
			if (q.status[egclist[EGCid].mh] == 2)
			{//is in Ud
				c_t += egclist[EGCid].cost;
				for (list<Hvertex>::iterator it = Hv.begin(); it != Hv.end(); )
				{
					if (egclist[EGCid].EGC2[(*it).id] && (*it).isud)
					{//vertex in EGC and vertex is in ud
						HUq--;
						HVweight -= instance.mix[(*it).id].vertex_weight;
						Hmix[(*it).id] = 0;
						it = Hv.erase(it);
					}
					else
					{
						it++;
					}
				}
				for (int i = 0; i < He2.size(); i++)
				{
					int e = instance.mix[He2[i]].e;
					bool headt = egclist[EGCid].EGC2[instance.one_d[e].head] && (q.status[instance.one_d[e].head] != 1);
					bool tailt = egclist[EGCid].EGC2[instance.one_d[e].tail] && (q.status[instance.one_d[e].tail] != 1);
					bool edt = egclist[EGCid].EGC2[He2[i]] && Heud[instance.mix[He2[i]].e];
					if (edt)
					{//edge in EGC and edge is in ud
						Hmix[He2[i]] = 0;
						HUq--;
						He2[i] = He2.back();
						He2.pop_back();
						i--;
					}
					else if (headt || tailt)
					{//head or tail in EGC and v is in ud
						Hmix[He2[i]] = 0;
						HUq--;
						He2[i] = He2.back();
						He2.pop_back();
						i--;
					}
				}
			}
			else
			{
				double c_min = 100 * instance.mix.size();
				for (list<Hvertex>::iterator it = Hv.begin(); it != Hv.end(); )
				{
					if (egclist[EGCid].EGC2[(*it).id] && (*it).isud)
					{//vertex in EGC and vertex is in ud
						if (instance.mix[(*it).id].blocking_cost < c_min)
						{
							c_min = instance.mix[(*it).id].blocking_cost;
						}
						HUq--;
						HVweight -= instance.mix[(*it).id].vertex_weight;
						Hmix[(*it).id] = 0;
						it = Hv.erase(it);
					}
					else
					{
						it++;
					}
				}
				for (int i = 0; i < He2.size(); i++)
				{
					int e = instance.mix[He2[i]].e;
					bool headt = egclist[EGCid].EGC2[instance.one_d[e].head] && (q.status[instance.one_d[e].head] != 1);
					bool tailt = egclist[EGCid].EGC2[instance.one_d[e].tail] && (q.status[instance.one_d[e].tail] != 1);
					bool edt = egclist[EGCid].EGC2[He2[i]] && Heud[instance.mix[He2[i]].e];
					if (edt)
					{//edge in EGC and edge is in ud
						Hmix[He2[i]] = 0;
						HUq--;
						if (instance.mix[He2[i]].blocking_cost < c_min)
						{
							c_min = instance.mix[He2[i]].blocking_cost;
						}
						He2[i] = He2.back();
						He2.pop_back();
						i--;
					}
					else if (headt || tailt)
					{//head or tail in EGC and v is in ud
						Hmix[He2[i]] = 0;
						HUq--;
						He2[i] = He2.back();
						He2.pop_back();
						i--;
					}
				}
				c_t += c_min;
			}
			EGCid = find_EGC_union(egclist, Hmix, EGCid + 1);
		}
		else
		{
			double WAV = 0;
			bool AUq = true;//elements in A that is in U is empty
			vector<int> AV2;//store vertices in format of vertices
			vector<int> AVm;//store vertices in format of mix
			vector<int> AEU;//store edges in format of mix in A
			list<Hvertex>::iterator it = Hv.begin();
			AV2.resize(instance.vertices_size, 0);
			int fv1 = -1;
			int fv2;
			int fe = -1;
			while (WAV <= instance.desired_clique_weight || AUq)
			{
				double tempw = instance.mix[(*it).id].vertex_weight;

				AV2[(*it).id] = 1;
				AVm.push_back((*it).id);
				WAV += tempw;
				if (AUq)
				{
					if ((*it).isud)
					{//is in undecided set
						AUq = false;
						fv1 = (*it).id;
						it++;
						continue;
					}
					int v = (*it).id;//v in vertex
					for (int i = 0; i < instance.neighbor[v].size(); i++)
					{
						int vp = instance.neighbor[v][i];
						int e = instance.two_d[v][vp];
						if (Heud[e] && AV2[vp] && Hmix[instance.one_d[e].IDm] > 0)
						{
							fe = instance.one_d[e].IDm;//thie edge need to be fix in egc
							fv1 = v;
							fv2 = vp;
							AUq = false;
							break;
						}
					}
				}
				it++;
			}
			if (check_gamma_un(instance, Hmix, AVm))
			{
				EGCid++;//make sure EGCid is not less than the size of egclist
				vector<bool> EGC = get_EGC_un(instance, AVm, Hmix);
				int egcpos = insert_EGC(instance, egclist, EGC);
				if (fe != -1)
				{
					EGC[fe] = true;
					EGC[fv1] = true;
					EGC[fv2] = true;
				}
				else if (fv1 != -1)
				{
					EGC[fv1] = true;
				}
				if (q.status[egclist[egcpos].mh] == 2)
				{//is Ud
					c_t += egclist[egcpos].cost;
					for (list<Hvertex>::iterator it = Hv.begin(); it != Hv.end(); )
					{
						if (EGC[(*it).id] && (*it).isud)
						{//vertex in EGC and vertex is in ud
							HUq--;
							HVweight -= instance.mix[(*it).id].vertex_weight;
							Hmix[(*it).id] = 0;
							it = Hv.erase(it);
						}
						else
						{
							it++;
						}
					}
					for (int i = 0; i < He2.size(); i++)
					{
						int e = instance.mix[He2[i]].e;
						bool headt = EGC[instance.one_d[e].head] && (q.status[instance.one_d[e].head] != 1);
						bool tailt = EGC[instance.one_d[e].tail] && (q.status[instance.one_d[e].tail] != 1);
						if (EGC[He2[i]] && Heud[instance.mix[He2[i]].e])
						{//edge in EGC and edge is in ud
							Hmix[He2[i]] = 0;
							HUq--;
							He2[i] = He2.back();
							He2.pop_back();
							i--;
						}
						else if (headt || tailt)
						{
							Hmix[He2[i]] = 0;
							HUq--;
							He2[i] = He2.back();
							He2.pop_back();
							i--;
						}
					}
				}
				else
				{
					double c_min = 100 * instance.mix.size();
					for (list<Hvertex>::iterator it = Hv.begin(); it != Hv.end(); )
					{
						if (EGC[(*it).id] && (*it).isud)
						{//vertex in EGC and vertex is in ud
							element temp = instance.mix[(*it).id];
							HUq--;
							HVweight -= temp.vertex_weight;
							Hmix[(*it).id] = 0;
							if (temp.blocking_cost < c_min)
							{
								c_min = temp.blocking_cost;//update c_min
							}
							it = Hv.erase(it);
						}
						else
						{
							it++;
						}
					}
					for (int i = 0; i < He2.size(); i++)
					{
						int e = instance.mix[He2[i]].e;
						bool headt = EGC[instance.one_d[e].head] && (q.status[instance.one_d[e].head] != 1);
						bool tailt = EGC[instance.one_d[e].tail] && (q.status[instance.one_d[e].tail] != 1);
						if (EGC[He2[i]] && Heud[instance.mix[He2[i]].e])
						{//edge in EGC and edge is in ud
							Hmix[He2[i]] = 0;
							HUq--;
							element temp = instance.mix[He2[i]];
							if (temp.blocking_cost < c_min)
							{
								c_min = temp.blocking_cost;//update c_min
							}
							He2[i] = He2.back();
							He2.pop_back();
							i--;
						}
						else if (headt || tailt)
						{
							Hmix[He2[i]] = 0;
							HUq--;
							He2[i] = He2.back();
							He2.pop_back();
							i--;
						}
					}
					if (c_min != 100 * instance.mix.size())
					{
						c_t += c_min;
					}
				}
			}
			else
			{
				for (list<Hvertex>::iterator it = Hv.begin(); it != Hv.end(); )
				{
					if ((*it).isud && AV2[(*it).id])
					{//vertex in EGC and vertex is in ud
						HUq--;
						HVweight -= instance.mix[(*it).id].vertex_weight;
						Hmix[(*it).id] = 0;
						it = Hv.erase(it);
					}
					else
					{
						it++;
					}
				}
				for (int i = 0; i < He2.size(); i++)
				{
					int e = instance.mix[He2[i]].e;
					bool headt = AV2[instance.one_d[e].head] && (q.status[instance.one_d[e].head] != 1);
					bool tailt = AV2[instance.one_d[e].tail] && (q.status[instance.one_d[e].tail] != 1);
					bool edt = AV2[instance.one_d[e].head] && AV2[instance.one_d[e].tail] && Heud[e];
					if (headt || tailt || edt)
					{//edge in EGC and edge is in ud
						Hmix[He2[i]] = 0;
						HUq--;
						He2[i] = He2.back();
						He2.pop_back();
						i--;
					}
				}
			}
		}
	}
	return c_t + q.solution_cost;
}

int insert_EGC(WKCI_instance& instance, vector<EGC>& egclist, vector<bool>& egc)
{
	EGC egc1;
	egc1.cost = 100;
	egc1.mh = 0;
	egc1.EGC2 = egc;
	for (int i = 0; i < egc.size(); i++)
	{
		if (egc[i])
		{
			egc1.EGC1.push_back(i);
			if (instance.mix[i].blocking_cost < egc1.cost)
			{
				egc1.cost = instance.mix[i].blocking_cost;
				egc1.mh = i;
			}
		}
	}
	int i = egclist.size() - 1;
	egclist.resize(egclist.size() + 1);
	while (i >= 0)
	{
		if (egclist[i].cost < egc1.cost)
		{
			egclist[i + 1] = egclist[i];
			i--;
		}
		else
		{
			break;
		}
	}
	egclist[i + 1] = egc1;
	return i + 1;
}

double LB_unprocessed_nodes(WKCI_instance& instance, int cn, vector<node>& tree)
{//cn+1 passed
	double LB = instance.mix.size() * 10;
	for (int i = 0; i < cn; i++)
	{
		if (tree[i].LB < LB)
		{
			LB = tree[i].LB;
		}
	}
	return LB;
}

//descending the list vertex weight in lowerbound
bool d_v_w(const Hvertex& a, const Hvertex& b)
{
	return (a.score > b.score);
}
//ascending ratio
void a_sort_ratio(vector<pair<int, double>>& candidate)
{
	sort(candidate.begin(), candidate.end(), [](const pair<int, double>& f, const pair<int, double>& s) { return f.second < s.second; });
}
//dscending tre nodes based on LB
void d_sort_lb(vector<node>& candidate)
{
	sort(candidate.begin(), candidate.end(), [](const node& f, const node& s) { return f.LB > s.LB; });
}

