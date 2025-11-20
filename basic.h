#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <set>

using namespace std;

int cccccccccc = -2;

#define Pop(stack) stack[--stack##_fill_pointer]
#define Push(item, stack) stack[stack##_fill_pointer++] = item
double sum_store_ver;
//#define NDEBUG

typedef long long llong;
typedef unsigned int uint;

struct two_ele {
	int v;
	int degree;
};

struct Edge {
  int v1;
  int v2;
};

chrono::steady_clock::time_point start;

llong max_steps;
llong step;
int try_step;
uint seed;
int cutoff_time;
int mode;
int flag_1000_out = 1;
int v_num;
int e_num;

Edge *edge;
int *edge_weight;

int *score;
// int *sub_score;
llong *time_stamp;
double s_alpha;
double s_beta;
double s_gamma;

int **v_adj;
int *v_degree;
int *v_threshold;
int *v_fixed;

int *dominated;
int c_size;
int *v_in_c;
int *candidate;
int *index_in_candidate;
int candidate_size;
llong now_weight;

int best_c_size;
int *best_v_in_c;
int *best_dominated;
double best_comp_time;
llong best_step;

int *undom_stack;
int undom_stack_fill_pointer;
int *index_in_undom_stack;

int *conf_change;
int *conf_count;

int ave_weight;
int s_gamma_total_weight;
int threshold;
double p_scale;

int para1;
int para2;
int para3;
int para4;
int para;

llong check_size = 0;


int *must_in;
int *record_in;
int must_in_beg = 0;
int must_in_end = 0;
int *I;

//added data structure
int unimprove;


int* zero_stack;
int* index_in_zero_stack;
int zero_stack_fill_pointer;
int* cover_set;
int* cover_set2;
string file;

//顶点开邻居中未被覆盖的个数
int *v_degree_tmp;
//闭邻居record_in未赋值为0的个数
int *v_degree_record;
int* is_covered;
int* pos_in_must_in;

int *propagate_degree;//
int *propagate_valid;//没有放入sort_propagate_degree的点，propagate_degree较大
int *index_in_propagate_valid;
int propagate_valid_fill_pointer = 0;

int unlock_propagate_ver_num;//还剩多少顶点没有加入sort_propagate_degree
int unlock_propagate_edge_num;

int max_num = 0;
int first_flag = 1;
int fix_flag = 1;
vector<vector <int> > store;
int* store_ver;//每个点的传播度是多少

int Max_degree = -1;
int store_count = 0;
int reduction_num;

int* temp_store_reduction;
int temp_store_reduction_fill_pointer;
int* temp_flag;
int* temp_flag_stack;
int temp_flag_stack_fill_pointer;

int thre5;
int thre32;

//reduction
vector<two_ele> v_d_vec;
bool* addable;//addable[v]=0，代表v已经在redcution列表中
int* neighbor_indicator;
int* neighbor_indicator2;//u的二邻居
int* working_vertice;
int* pos_in_working_vertice;
int* next_working_vertice;
int num_working_vertice;
int num_next_working_vertice;
int* new_fix;
int* pos_in_new_fix;
int num_new_fix;
int* new_fix_neightbor;
int* pos_in_new_fix_neightbor;
int num_new_fix_neightbor;
long long fix_time_stamp;//fix时间戳
long long* time_out;//time_out[i]代表i出working_vertice数组的时间
long long* time_need_in;//需要加入next_working_vertice数组的时间
long long* fix_time;//fix_time[i]代表i的fix为1的时间，也是i的邻居邻居因为i而需要加入next_working_vertice的时间
double t_redu;
int num_redu;
int u_neightbor;//顶点u未被覆盖邻居的个数
int u_neightbor2;

//covers_set
bool cover_set_flag;//为true代表使用cover_set
vector<int> fix_must_in;
int best_candidate_record;
int* best_candidate;
vector<int> equal_Reduction;
bool* need_check_X;

//Ignore
int checkI = 0;
int flagX = -3, flagI = -3;
int flagX1 = -3, flagI1 = -3;
int countccc;
bool* need_check_I;
int num_new_X;
int* new_X;
bool* new_covered;

//delete vertex
vector<int>delete_vertex;

//reduction rule 4
int* N2;
int num_N2 = 0;
int* pos_in_N2;
bool* is_N2_neightbor;

//large graph
bool is_large = false;
bool rule4_is_large = false;
