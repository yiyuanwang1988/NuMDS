#pragma once

#include "my_heap.h"
// #define NDEBUG

double TimeElapsed();
int BuildInstance(string);
void FreeMemory();
inline void Undom(int);
inline void Dom(int);
void Add(int);
void Remove(int);
void ConstructDS();
void ConstructByDegree();
void ConstructByScore();
void ResetCandidate();
void UpdateBestSolution();
int CheckSolution();
void UpdateTargetSize();
int ChooseRemoveV(int);
int ChooseAddV();
void AddInit();
void init_increase_dominate(int, int);
void init_decrease_dominate(int);
void init_removeRedundant();
void removeRedundant();
void LocalSearch();


int t1 = 100; //25 + rand() % 10
int t2 = 50;  //1

int total_count = 0;
int total_32 = 0;
int total_5 = 0;

double TimeElapsed() {
	chrono::steady_clock::time_point finish = chrono::steady_clock::now();
	chrono::duration<double> duration = finish - start;
	return duration.count();
}

void swap(int x, int low, int high) //交换两个数的值
{
	if (low == high)
		return;
	int t = v_adj[x][low];    v_adj[x][low] = v_adj[x][high];    v_adj[x][high] = t;
}

void adjust(int x, int index, int len)
{
	int leftIndex = 2 * index + 1;
	int rightIndex = 2 * index + 2;
	int bigIndex = index;

	if (leftIndex < len && v_adj[x][bigIndex] < v_adj[x][leftIndex]) {
		bigIndex = leftIndex;
	}

	if (rightIndex < len && v_adj[x][bigIndex] < v_adj[x][rightIndex]) {
		bigIndex = rightIndex;
	}

	if (bigIndex != index) {
		swap(x, index, bigIndex);
		adjust(x, bigIndex, len);
	}
}

void heapsort(int x, int low, int high)
{
	int len = high - low + 1;
	for (int i = (len / 2) - 1; i >= 0; i--) {
		adjust(x, i, len);
	}

	for (int i = 0; i < len; i++) {
		swap(x, 0, len - 1 - i);
		adjust(x, 0, len - 1 - i);
	}
}

int partition(int x, int low, int high)  //计算基准点，分割为左右两个数组
{
	int point = v_adj[x][low];//基准点等于第一个元素
	while (low < high)
	{
		while (low < high && v_adj[x][high] >= point)//控制high指针比较并左移
		{
			high--;
		}
		swap(x, low, high);
		//}
		while (low < high && v_adj[x][low] <= point)//控制low指针比较并右移
		{
			low++;
		}
		swap(x, low, high);
	}
	return low;//返回基准点位置
}

void quicksort(int x, int low, int high)  //low:起始位置 high:末尾位置
{
	if (low < high) {
		int point = partition(x, low, high);//计算基准点
		quicksort(x, low, point - 1);  //对基准点的左边进行排序
		quicksort(x, point + 1, high);//对基准点的右边进行排序
	}
}

bool CompareBySecond2(const two_ele& lhs, const two_ele& rhs) {
	return lhs.degree > rhs.degree;
}

bool CompareBySecond(const two_ele &lhs, const two_ele &rhs) {
	return lhs.degree < rhs.degree;
}

int get_index_v_adj(int v, int u) {//v的哪个邻居是v_adj 索引从0到v_degree - 1
	int low = 0;
	int high = v_degree[v] - 1;
	int mid = (low + high) / 2;
	while (low <= high) {
		if (v_adj[v][mid] == u) {
			return mid;
		}
		else if (v_adj[v][mid] > u) {
			high = mid - 1;
		}
		else {
			low = mid + 1;
		}
		mid = (low + high) / 2;
	}
	return -1;
}

int get_index_v_adj_range(int v, int u, int lowest) {//v的哪个邻居是v_adj 索引从0到v_degree - 1
	int low = lowest;
	int high = v_degree[v] - 1;
	int mid = (low + high) / 2;
	while (low <= high) {
		if (v_adj[v][mid] == u) {
			return mid;
		}
		else if (v_adj[v][mid] > u) {
			high = mid - 1;
		}
		else {
			low = mid + 1;
		}
		mid = (low + high) / 2;
	}
	return -1;
}

int belong(int u, int v)
{
	int i, j;
	if (v_degree[u] <= v_degree[v])
	{
		for (i = 0, j = 0; i < v_degree[u] && j < v_degree[v]; )
		{
			if (is_covered[i] == 1) {
				i++;
				continue;
			}
			if (is_covered[j] == 1) {
				j++;
				continue;
			}
			if (v_adj[u][i] < v_adj[v][j])  return -2;//N[u]不属于N[v] 并且 N[v]不属于N[u]
			if (v_adj[u][i] == v_adj[v][j])  i++, j++;
			else if (v_adj[u][i] > v_adj[v][j])  j++;
		}
		if (i == v_degree[u] && j == v_degree[v] && i == j)  return 0;//N[u]==N[v]
		if (i == v_degree[u])  return 1;//N[u]属于N[v]
		return -2;//N[u]不属于N[v] 并且 N[v]不属于N[u]
	}
	else {
		for (j = 0, i = 0; j < v_degree[v] && i < v_degree[u];)//看N[u]是否为N[v]的超集//u的度为17 v的度为2
		{
			if (is_covered[i] == 1) {
				i++;
				continue;
			}
			if (is_covered[j] == 1) {
				j++;
				continue;
			}
			if (v_adj[v][j] < v_adj[u][i])  return -2;//N[u]不属于N[v] 并且 N[v]不属于N[u]
			if (v_adj[v][j] == v_adj[u][i])  i++, j++;
			else if (v_adj[v][j] > v_adj[u][i])  i++;
		}
		if (j == v_degree[v])  return -1;//N[v]属于N[u]
		return -2;//N[u]不属于N[v] 并且 N[v]不属于N[u]
	}
}

int belong_1(int u, int v)
{
	int i, j;
	int v_de_tmp_u = v_degree_tmp[u];
	int v_de_tmp_v = v_degree_tmp[v];
	v_de_tmp_u--;
	v_de_tmp_v--;
	if (is_covered[v] == 1 || I[v] == 1)
		v_de_tmp_u++;
	if (is_covered[u] == 1 || I[v] == 1)
		v_de_tmp_v++;
	if (v_de_tmp_u < v_de_tmp_v)//这种情况需要看v是否为u的超集
	{
		int lowest_111 = 0;
		for (int i = 0; i < v_degree[u]; i++) {//时间复杂度v_degree_tmp[u]*log(v_degree[v])
			int u_adjcent = v_adj[u][i];
			if (is_covered[u_adjcent] == 1 || u_adjcent == v) {
				continue;
			}
			if (v_degree[u_adjcent] > v_degree[v]) {//u_adjcent的度大于v，比如则查看v的v_adj有没有u;
				int temp111 = get_index_v_adj_range(v, u_adjcent, 0);//查看v的邻居有没有u_adjcent；
				if (temp111 == -1)
					return -2;
				if (temp111 < 0 || temp111 > v_degree[v]) {
					cout << file << " temp111 < 0 || temp111 > v_degree[v] " << endl;
				}
				if (v_adj[v][temp111] != u_adjcent) {
					cout << " v_adj[v][temp111] != u_adjcent " << endl;
				}
				lowest_111 = temp111;//这里要求v_adj的各个顶点根据度从小到大排列
			}
			else {//v的度大于u_adjcent
				int temp111 = get_index_v_adj_range(u_adjcent, v, 0);//查看v的邻居有没有u_adjcent；
				if (temp111 == -1)
					return -2;
				if (temp111 < 0 || temp111 > v_degree[u_adjcent]) {
					cout << file << " temp111 < 0 || temp111 > v_degree[u] " << endl;
				}
				if (v_adj[u_adjcent][temp111] != v) {
					cout << " v_adj[u][temp111] != v " << endl;//若u_adjcent的邻居不为v
				}
				lowest_111 = temp111;
			}
		}
		return 1;//u_adjcent是v的子集
	}
	else if (v_de_tmp_u > v_de_tmp_v) {
		int lowest_111 = 0;
		for (int i = 0; i < v_degree[v]; i++) {
			int v_adjcent = v_adj[v][i];
			if (is_covered[v_adjcent] == 1 || v_adjcent == u) {
				continue;
			}
			if (v_degree[u] > v_degree[v_adjcent]) {//u的度大于v_adjcent，比如则查看v_adjcent的v_adj有没有u；
				int temp111 = get_index_v_adj_range(v_adjcent, u, 0);//查看v_adjcent的邻居有没有u；
				if (temp111 == -1)
					return -2;
				if (temp111 < 0 || temp111 > v_degree[v_adjcent]) {
					cout << file << " temp111 < 0 || temp111 > v_degree[v_adjcent] " << endl;
				}
				if (v_adj[v_adjcent][temp111] != u) {
					cout << " v_adj[v_adjcent][temp111] != u " << endl;
				}
				lowest_111 = temp111;
			}
			else {
				int temp111 = get_index_v_adj_range(u, v_adjcent, 0);//查看v_adjcent的邻居有没有u；
				if (temp111 == -1)
					return -2;
				if (temp111 < 0 || temp111 > v_degree[u]) {
					cout << file << " temp111 < 0 || temp111 > v_degree[u] " << endl;
				}
				if (v_adj[u][temp111] != v_adjcent) {
					cout << " v_adj[u][temp111] != v_adjcent " << endl;//若u的邻居不为v_adjcent
				}
				lowest_111 = temp111;
			}
		}
		return -1;//v_adjcent是u的子集
	}
	else {
		int lowest_111 = 0;
		for (int i = 0; i < v_degree[v]; i++) {
			int v_adjcent = v_adj[v][i];
			if (is_covered[v_adjcent] == 1 || v_adjcent == u) {
				continue;
			}
			if (v_degree[u] > v_degree[v_adjcent]) {//u的度大于v_adjcent，比如则查看v_adjcent的v_adj有没有u；
				int temp111 = get_index_v_adj_range(v_adjcent, u, 0);//查看v_adjcent的邻居有没有u；
				if (temp111 == -1)
					return -2;
				if (temp111 < 0 || temp111 > v_degree[v_adjcent]) {
					cout << file << " temp111 < 0 || temp111 > v_degree[v_adjcent] " << endl;
				}
				if (v_adj[v_adjcent][temp111] != u) {
					cout << " v_adj[v_adjcent][temp111] != u " << endl;
				}
				lowest_111 = temp111;
			}
			else {
				int temp111 = get_index_v_adj_range(u, v_adjcent, 0);//查看v_adjcent的邻居有没有u；
				if (temp111 == -1)
					return -2;
				if (temp111 < 0 || temp111 > v_degree[u]) {
					cout << file << " temp111 < 0 || temp111 > v_degree[u] " << endl;
				}
				if (v_adj[u][temp111] != v_adjcent) {
					cout << " v_adj[u][temp111] != v_adjcent " << endl;//若u的邻居不为v_adjcent
				}
				lowest_111 = temp111;
			}
		}
		return 0;//u于v_adjcent等价
	}
}

void check_v_degree_tmp(int ver) {
	int count = 0;
	for (int j = 0; j < v_degree[ver]; j++) {
		int v1 = v_adj[ver][j];
		if (is_covered[v1] == 0 && I[v1] == 0)
			count++;
	}
	if (count != v_degree_tmp[ver]) {
		cout << file << " count != v_degree_tmp[ver] " << ver << endl;
	}
}

void check_v_record(int ver) {
	int count = 0;
	for (int j = 0; j < v_degree[ver]; j++) {
		int v1 = v_adj[ver][j];
		if (record_in[v1] != 0)
			count++;
	}
	if (record_in[ver] != 0)
		count++;
	if (count != v_degree_record[ver] && is_covered[ver] == 0) {
		cout << file << " count != v_degree_record[ver] " << ver << endl;
	}
}

void check_2(int a) {
	for (int i = 0; i < v_degree[a]; i++) {
		int vadj = v_adj[a][i];
		if (record_in[vadj] == -1) {
			check_v_record(vadj);
			check_v_degree_tmp(vadj);
		}
	}
	if (record_in[a] == -1) {
		check_v_record(a);
		check_v_degree_tmp(a);
	}
	if (v_degree[a] > 1)
		for (int i = 1; i <= v_num; i++) {
			if (v_degree[i] != v_degree[i])
				cout << "i" << i << endl;
		}
}

void check4633(int k) {
	for (int i = 1; i <= v_num; i++) {
		if (v_degree[i] != v_degree[i])
			cout << "i " << k << endl;
	}
}

int get_de_tmp(int v, int u) {
	int count = 0;
	for (int i = 0; i < v_degree[v]; i++) {
		if (v_adj[v][i] != u && is_covered[v_adj[v][i]] == 0)
			count++;
	}
	return count;
}

void check_domi() {
	for (int i = 1; i <= v_num; i++) {
		int domi;
		if (v_in_c[i] == 0) {
			domi = 0;
		}
		else {
			domi = 1;
		}

		for (int j = 0; j < v_degree[i]; j++) {
			int v_adj_ver = v_adj[i][j];
			if (v_in_c[v_adj_ver] == 1) {
				domi++;
			}
		}
		if (domi != dominated[i])
			cout <<"(domi != dominated[i])"<< i << dominated[i] << endl;
	}
}


void check_undomi() {
	for (int i = 1; i <= v_num; i++) {
		if (dominated[i] == 0) {
			int flag = 0;
			for (int j = 0; j < undom_stack_fill_pointer; j++) {
				if (undom_stack[j] == i) {
					flag = 1;
					break;
				}

			}
			if (flag == 0)
				cout << "not find undomi" << endl;
		}
	}
}

bool check_sort(int v) {
	for (int i = 0; i < v_degree[v] - 1; i++) {
		if (v_adj[v][i] > v_adj[v][i + 1])
			return true;
	}
	return false;
}

void check_candidate() {
	for (int i = 1; i <= v_num; i++) {
		if (v_in_c[i] == 1 && record_in[i] != 1 && index_in_candidate[i] < 0) {
			cout << " v_in_c[i] == 1 && v_fixed[i] == 0 && index_in_candidate[i] < 0 " << i << endl;
		}
		if (index_in_candidate[i] >= 0 && candidate[index_in_candidate[i]] != i)
			cout << "index_in_candidate[i] >= 0 && candidate[index_in_candidate[i]]" << i << " " << file << endl;
	}
	for (int i = 0; i < candidate_size; i++) {
		if (candidate[i] < 0)
			cout << "candidate[i] < 0 " << candidate[i] << endl;
		if (v_in_c[candidate[i]] == 0 || record_in[candidate[i]] == 1)
			cout << " v_in_c[i] == 0 || v_fixed[i] == 1 " << i << endl;
	}
}


void check_candidate1() {
	for (int i = 0; i < candidate_size; i++) {
		if (candidate[i] == 511833) {
			cout << "wrong" << endl;
		}
	}
}

void add_zero_vertex(int v) {
	zero_stack[zero_stack_fill_pointer] = v;
	index_in_zero_stack[v] = zero_stack_fill_pointer;
	zero_stack_fill_pointer++;
}

void delete_zero_vertex(int v) {
	int vpos = index_in_zero_stack[v];
	int last_v = zero_stack[zero_stack_fill_pointer - 1];
	zero_stack[vpos] = last_v;
	index_in_zero_stack[last_v] = vpos;
	index_in_zero_stack[v] = -1;
	zero_stack_fill_pointer--;
}

void push_propagate_valid(int ver) {//预处理前所有顶点都是valid的，逐渐变为非valid并向sort_propagate_degree中移入
	propagate_valid[propagate_valid_fill_pointer] = ver;
	index_in_propagate_valid[ver] = propagate_valid_fill_pointer;
	propagate_valid_fill_pointer++;
}

void add_to_new_X(int u)
{
	new_X[num_new_X++] = u;
}

void add_to_new_fix_neightbor(int u)
{
	if (pos_in_new_fix_neightbor[u] == -1)
	{
		pos_in_new_fix_neightbor[u] = num_new_fix_neightbor;
		new_fix_neightbor[num_new_fix_neightbor++] = u;
	}
}

void delete_fix_edges(int u)
{
	int i, j, k, v;
	for (i = 0; i < v_degree[u]; i++)
	{
		v = v_adj[u][i];//v是u的邻居
		if (v_degree[v] == 0)  continue;
		for (j = 0; j < v_degree[v]; j++)
		{
			if (v_adj[v][j] == u)  break;
		}
		if (j == v_degree[v])
		{
			printf("error: j == v_degree[v]\n");
			printf("v_degree[v]=%d\n", v_degree[v]);
			exit(0);
		}
		v_adj[v][j] = v_adj[v][v_degree[v]-1];
		v_degree[v]--;
	}
	v_degree[u] = 0;
	//delete v_adj[u];
}

void update_next_working_vertice()
{
	for (int ii = 0, u, v, kk; ii < num_new_X; ii++)
	{
		u = new_X[ii];//u在上一轮循环中record_in被置为0
		for (kk = 0; kk < v_degree[u]; kk++)
		{
			v = v_adj[u][kk];
			if (addable[v])
			{
				next_working_vertice[num_next_working_vertice] = v;
				num_next_working_vertice++;
				addable[v] = false;
				//need_check_I[v] = true;
			}
			need_check_I[v] = true;
		}
		//if (is_covered[u] == 1)  delete_vertex.push_back(u);
	}
	for (int ii = 0, u, v, kk; ii < num_new_fix; ii++)//???
	{
		u = new_fix[ii];//u在上一轮循环中被固定
		for (kk = 0; kk < v_degree[u]; kk++)
		{
			v = v_adj[u][kk];
			if (new_covered[v])//只有当v是上一轮刚刚被覆盖时，他的邻居才可以加入下次的检查队列
			{
				add_to_new_fix_neightbor(v);
				new_covered[v] = false;//重置
				//if (record_in[v] == 0)  delete_vertex.push_back(v);
			}
			if (new_covered[u])//如果u是上一轮新被覆盖时，他的邻居才可以加入下次的检查队列
			{
				if (addable[v])//因为被固定为1的顶点的is_covered也被置为1，所以在这里固定顶点邻居也应该加入后续的work数组中
				{
					next_working_vertice[num_next_working_vertice] = v;
					num_next_working_vertice++;
					addable[v] = false;
					//need_check_X[v] = true;
				}
				need_check_X[v] = true;
			}
		}
		new_covered[u] = false;//重置
		pos_in_new_fix[u] = -1;
	}
	for (int ii = 0, u, v, kk; ii < num_new_fix_neightbor; ii++)
	{
		u = new_fix_neightbor[ii];//u的邻居需要加入next_working_vertice中
		for (kk = 0; kk < v_degree[u]; kk++)
		{
			v = v_adj[u][kk];
			if (addable[v])
			{
				next_working_vertice[num_next_working_vertice] = v;
				num_next_working_vertice++;
				addable[v] = false;
			}
			need_check_X[v] = true;
		}
		pos_in_new_fix_neightbor[u] = -1;
	}

	for (int ii = 0, u, v, kk; ii < num_new_X; ii++)
	{
		u = new_fix[ii];//u在上一轮循环中被固定
		//delete_fix_edges(u);//v1固定后，将与v1连接的边删除
	}

	num_new_fix = 0;
	num_new_fix_neightbor = 0;
	num_new_X = 0;
	/*
	for (int ii = 0, u; ii < delete_vertex.size(); ii++)
	{
		u = delete_vertex[ii];
		if (record_in[u] != 0 || is_covered[u] != 1)
		{
			cout << file << " ";
			printf("record_in[u] != 0 || is_covered[u] != 1\n");
			exit(0);
		}
		delete_fix_edges(u);
	}
	delete_vertex.clear();
	*/
}

void add_to_must_in(int u)
{
	record_in[u] = 1;
	if (pos_in_must_in[u] == -1)
	{
		pos_in_must_in[u] = must_in_end;
		must_in[must_in_end++] = u;
		fix_must_in.push_back(u);
	}
	else {
		cout << file << " add to must in wrong" << seed << endl;
	}
}

void add_to_new_fix(int u)
{
	if (pos_in_new_fix[u] == -1)
	{
		pos_in_new_fix[u] = num_new_fix;
		new_fix[num_new_fix++] = u;
		fix_time[u] = fix_time_stamp;//记录u被固定的时间
	}
}

void update_v_degree_tmp(int v1)
{
	//printf("update degree ---------\n");
	for (int ii = 0; ii < v_degree[v1]; ii++) {
		int vv = v_adj[v1][ii];
		if (is_covered[vv] == 0 ) {
			is_covered[vv] = 1;//v的唯一邻居n被固定为1后，n的邻居一定是is cover 
			for (int j = 0; j < v_degree[vv]; j++) {//若顶点v是is cover的，则v的邻居需要覆盖的顶点减少，若v_degree_temp为0则说明该顶点可以被fix为0
				int vv1 = v_adj[vv][j];//此外v_degree_tmp大小还可以被用来估计俩点之间是否存在包含关系
				v_degree_tmp[vv1]--;
			}
			new_covered[vv] = true;
		}
	}
	if (is_covered[v1] == 0 ) {
		is_covered[v1] = 1;
		for (int j = 0; j < v_degree[v1]; j++) {//若顶点v是is cover的，则v的邻居需要覆盖的顶点减少，若v_degree_tmp为0则说明该顶点可以被fix为0
			int v11 = v_adj[v1][j];//此外v_degree_tmp大小还可以被用来估计俩点之间是否存在包含关系。
			v_degree_tmp[v11]--;//v_degree_tmp[v1]需要跟v1的邻居中is_coverd为0的顶点个数相等。
		}
		new_covered[v1] = true;
	}
}

void fix_only_neighbor1(int u)
{
	int v1 = -1;
	for (int j = 0; j < v_degree[u]; j++) {
		int verradj = v_adj[u][j];
		if (record_in[verradj] == -1) {
			v1 = verradj;
			break;
		}
	}
	if (v1 == -1)
		cout << file << " " << v1 << "v1 == -1" << endl;
	if (record_in[v1] == 0)
		cout << "record_in[v1] == 0 3 " << endl;
	add_to_must_in(v1);
	update_v_degree_tmp(v1);
	add_to_new_fix(v1);
	//delete_fix_edges(v1);//v1固定后，将与v1连接的边删除
}

void fix_only_neighbor2(int u)
{
	int v1 = -1;
	for (int j = 0; j < v_degree[u]; j++) {
		int verradj = v_adj[u][j];
		if (record_in[verradj] == -1) {
			v1 = verradj;
			break;
		}
	}
	if (v1 == -1 && record_in[u] == -1)
		v1 = u;
	else if (v1 == -1 && record_in[u] != -1)
		cout << file << " " << v1 << "v1 == -1" << endl;
	if (record_in[v1] == 0)
		cout << "record_in[v1] == 0 3 " << endl;
	add_to_must_in(v1);
	update_v_degree_tmp(v1);
	add_to_new_fix(v1);
	//delete_fix_edges(v1);//v1固定后，将与v1连接的边删除
}

int new_check_belong(int u, int v)
{
	int v_de_tmp_u = v_degree_tmp[u];
	int v_de_tmp_v = v_degree_tmp[v];
	v_de_tmp_u--;
	v_de_tmp_v--;
	if (is_covered[v] == 1 || I[v] == 1)  v_de_tmp_u++;
	if (is_covered[u] == 1 || I[u] == 1)  v_de_tmp_v++;
	//printf("111:v_de_tmp_u=%d  v_de_tmp_v=%d\n", v_de_tmp_u, v_de_tmp_v);
	//if (is_covered[u] == 0)  v_de_tmp_u++;
	//if (is_covered[v] == 0)  v_de_tmp_v++;
	if (v_de_tmp_u == v_de_tmp_v) {
		if (v_de_tmp_u == 0)
			return 0;

		bool is_dominated = false;
		for (int i = 0, j; i < v_degree[v]; i++)
		{
			j = v_adj[v][i];
			if (is_covered[j] == 1 || j == u || I[j] == 1)  continue;
			if (neighbor_indicator[j] == false)
			{
				is_dominated = true;
				break;
			}
		}
		if (is_dominated)  return -2;

		if (record_in[u] == 1) {
			equal_Reduction.push_back(v);
			return -1;//N[v]属于N[u],将record_in[v]置为0
		}
		if (record_in[v] == 1) {
			equal_Reduction.push_back(u);
			return 1;//N[u]属于N[v],将record_in[u]置为0
		}
		int minu = v_num + 1, minv = v_num + 1;
		//if (record_in[u] == -1)  minu = min(minu, v_degree_record[u]);
		//if (record_in[v] == -1)  minv = min(minv, v_degree_record[v]);
		for (int ii = 0, jj; ii < v_degree[u]; ii++)
		{
			jj = v_adj[u][ii];
			if (record_in[jj] != -1)  continue;
			//if (jj == v)  continue;
			if (minu > v_degree_record[jj])
				minu = v_degree_record[jj];
		}

		for (int ii = 0, jj; ii < v_degree[v]; ii++)
		{
			jj = v_adj[v][ii];
			if (record_in[jj] != -1)  continue;
			//if (jj == u)  continue;
			if (minv > v_degree_record[jj])
				minv = v_degree_record[jj];
		}
		if (minu < minv) {
			equal_Reduction.push_back(u);
			return 1;//N[u]属于N[v],将record_in[u]置为0
		}
		if (minu > minv) {
			equal_Reduction.push_back(v);
			return -1;//N[v]属于N[u],将record_in[v]置为0
		}
		return 0;

	}//X
	else if (v_de_tmp_u < v_de_tmp_v)//这里只能是u为v的子集
	{
		//if (record_in[u] == 0)  return -2;
		//printf("--------\n");
		//return -2;
		if (v_de_tmp_u == 0)
			return 1;

		int v_neightbor = 0;
		if (is_covered[v] == 0 && I[v] == 0)  v_neightbor++;
		for (int i = 0, j; i < v_degree[v]; i++)
		{
			j = v_adj[v][i];
			if (is_covered[j] == 1 || j == u || I[j] == 1)  continue;
			if (neighbor_indicator[j] == true)
				v_neightbor++;
		}
		//printf("flag:  v_neightbor=%d  u_neightbor=%d\n", v_neightbor, u_neightbor);
		if (v_neightbor != u_neightbor)  return -2;
		return 1;

	}
	else {
		//if (record_in[v] == 0)  return -2;
		if (v_de_tmp_v == 0)
			return -1;

		bool is_dominated = false;
		for (int i = 0, j; i < v_degree[v]; i++)
		{
			j = v_adj[v][i];
			if (is_covered[j] == 1 || j == u || I[j] == 1)  continue;
			if (neighbor_indicator[j] == false)
			{
				is_dominated = true;
				break;
			}
		}
		if (is_dominated)  return -2;
		return -1;

	}
	return 0;
}

int check_I(int u, int v)
{
	
	//N[u]\X属于N[v]\X
	bool flag = false;
	for (int i = 0, j, x, y; i < v_degree[u]; i++)
	{
		x = v_adj[u][i];
		if (record_in[x] != -1)  continue;
		if (x == v)  continue;
		for (j = 0; j < v_degree[v]; j++)
		{
			y = v_adj[v][j];
			if (x == y)  break;
		}
		if (j == v_degree[v])
		{
			flag = true;
			break;//在v的邻居中没有找到u的邻居x，所以N[u]\X不属于N[v]\X
		}
	}
	if (flag)  return false;
	return true;
}

int check_I2(int u, int v)
{
	//N[u]\X属于N[v]\X
	bool flag = false;
	for (int i = 0, j, x, y; i < v_degree[u]; i++)
	{
		x = v_adj[u][i];
		if (record_in[x] != -1)  continue;
		if (x == v)  continue;
		for (j = 0; j < v_degree[v]; j++)
		{
			y = v_adj[v][j];
			if (x == y)  break;
		}
		if (j == v_degree[v])
		{
			flag = true;
			break;//在v的邻居中没有找到u的邻居x，所以N[u]\X不属于N[v]\X
		}
	}
	if (!flag)
	{
		int x = u, j, y;
		if (record_in[x] == -1)
		{
			for (j = 0; j < v_degree[v]; j++)
			{
				y = v_adj[v][j];
				if (x == y)  break;
			}
			if (j == v_degree[v])
			{
				flag = true;
			}
		}
	}
	if (flag)  return false;
	return true;
}

int check_X(int u, int v)
{
	if (record_in[u] == 0 || record_in[v] == 0 || record_in[u] == 1)  return false;
	//N[u]\N[S]属于N[v]\N[S]
	bool flag = false;
	bool is_zero = true;
	for (int i = 0, j, x, y; i < v_degree[u]; i++)
	{
		x = v_adj[u][i];
		if (is_covered[x])  continue;
		if (x == v)  continue;
		is_zero = false;
		for (j = 0; j < v_degree[v]; j++)
		{
			y = v_adj[v][j];
			if (x == y)  break;
		}
		if (j == v_degree[v])
		{
			flag = true;
			break;//在v的邻居中没有找到u的邻居x，所以N[u]\N[S]不属于N[v]\N[S]
		}
	}
	if (is_zero)  return true;

	if (flag)  return false;
	return true;
}

int new_check_belong1(int u, int v, int round)
{
	//重置flag
	flagX = -3, flagI = -3;
	
	//flagX=-2代表u和v不会加入入X中。flagI=-2代表u和v不会加入I中
	//计算N[u]\X和N[v]\X
	int u_not_in_X = v_degree_record[u];//N[u]\X
	int v_not_in_X = v_degree_record[v];
	if (is_covered[u] || is_covered[v])  flagI = -2;//判断I时，只需要判断两个未被支配的顶点即可
	if (record_in[u] == 0 || record_in[v] == 0)  flagX = -2;
	//计算N[u]\(N[S]并I)和N[v]\(N[S]并I)
	if (v_num > 500000) {
		if (v_degree_tmp[v] > 500 && v_degree_tmp[u] > 500) {
			flagX = -2;
		}
		if (v_degree_record[v] > 500 && v_degree_record[u] > 500) {
			flagI = -2;
		}
	}
	int v_de_tmp_u = v_degree_tmp[u];//N(u)\(N[S]并I)
	int v_de_tmp_v = v_degree_tmp[v];
	v_de_tmp_u--;
	v_de_tmp_v--;
	if (is_covered[v] == 1)  v_de_tmp_u++;
	if (is_covered[u] == 1)  v_de_tmp_v++;

	if (need_check_I[u] == 0 && need_check_I[v] == 0)  flagI = -2;//两个都不用检查I关系
	if (need_check_X[u] == 0 && need_check_X[v] == 0)  flagX = -2;//两个都不用检查X关系
	/*
	if (round > 1)
	{
		if (need_check_I[u] == 0 && need_check_I[v] == 0)  flagI = -2;//两个都不用检查I关系
		else if (need_check_I[u] == 1 && need_check_I[v] == 1);//两个都需要检查，就不做任何处理
		else if (need_check_I[v] == 1)//只有v需要检查I关系
		{
			
			if (v_not_in_X > u_not_in_X)//这里v只需要在N[v]\X <= N[u]\X的情况下进行检查
			{
				flagI = -2;
			}
			
		}
		else if (need_check_I[u] == 1)//只有u需要检查I关系
		{
			if (u_not_in_X > v_not_in_X)//这里v只需要在N[u]\X <= N[v]\X的情况下进行检查
			{
				flagI = -2;
			}
		}

		if (need_check_X[u] == 0 && need_check_X[v] == 0)  flagX = -2;//两个都不用检查X关系
		else if (need_check_X[u] == 1 && need_check_X[v] == 1);//两个都需要检查，就不做任何处理
		else if (need_check_X[v] == 1)//只有v需要检查X关系
		{
			if (v_de_tmp_v > v_de_tmp_u)//这里v只需要在N[v]\N[S] <= N[u]\N[S]的情况下进行检查
			{
				flagX = -2;
			}
		}
		else if (need_check_X[u] == 1)//只有u需要检查X关系
		{
			if (v_de_tmp_u > v_de_tmp_v)//这里v只需要在N[u]\N[S] <= N[v]\N[S]的情况下进行检查
			{
				flagX = -2;
			}
		}
	}
	*/
	if ((v_de_tmp_u < v_de_tmp_v && record_in[u] == 1))// || (v_degree_tmp[u] == 0 && record_in[u] == 1)
		flagX = -2;
	if ((v_de_tmp_u > v_de_tmp_v && record_in[v] == 1))// || (v_degree_tmp[v] == 0 && record_in[v] == 1)//若v被固定为1则v_degree_tmp[v]一定为0，
		flagX = -2;

	if (is_large && flagX == -2)
	{
		return 0;
	}
	
	u_neightbor = v_degree_tmp[u];//N(u)\(N[S])
	if (is_covered[u] == 0)  u_neightbor++;//N[u]\(N[S])
	int v_neightbor = 0, v_neightbor_I = 0;
	if (is_covered[v] == 0)  v_neightbor++;//因为u的邻居中肯定有v,v_neightbor表示N[v]中u的is_covered不为0邻居的个数,但下面遍历的是N(v)
	if (record_in[v] == -1)  v_neightbor_I++;//因为u的邻居中肯定有v,v_neightbor_I表示N[v]中u的record_in不为0邻居的个数,但下面遍历的是N(v)
	for (int i = 0, j; i < v_degree[v]; i++)
	{
		j = v_adj[v][i];
		//X: C(u) = N[u] \ (N[S] 并 I) 属于 N[v], add u to X
		if (flagX == -3)
		{
			if (v_de_tmp_u < v_de_tmp_v)//u可能为v的子集
			{
				if (is_covered[j] == 0 && neighbor_indicator[j] == true)
				{
					v_neightbor++;
					if (v_neightbor == u_neightbor)
					{
						flagX = 1;//X判断时，u的邻居全是v的邻居，代表N[u]\{N[S]并I}为N[v]的子集，需要将u加入X
					}
				}
				if (v_neightbor + v_degree[v] - i - 1 < u_neightbor)
				{
					flagX = -2;
				}

			}
			else {//N[v]可能为N[u]的子集
				if (is_covered[j] == 0 && neighbor_indicator[j] == false)
				{
					//flagX == -1代表还没有判断出u是否加入X
					//j属于N[v]\{N[S]并I}
					//v的邻居j不属于u的邻居，所以这次检查关系中u不会加入X
					flagX = -2;
				}
			}
		}

		//I: D(u) = N[u] \ X 属于 N[v], add v to I
		if (flagI == -3)
		{
			if (u_not_in_X < v_not_in_X)//N[u]\X可能为N[v]\X的子集
			{
				if (record_in[j] != 0 && neighbor_indicator[j] == true)
				{
					v_neightbor_I++;
					if (v_neightbor_I == u_not_in_X)
						flagI = -1;//I判断时，u的邻居全是v的邻居，代表N[u]\X为N[v]\X的子集，需要将v加入I
				}
				if (v_neightbor_I + v_degree[v] - i - 1 < u_not_in_X)
				{
					flagI = -2;
				}
			}
			else {//N[v]\X可能为N[u]\X的子集
				if (record_in[j] != 0 && neighbor_indicator[j] == false)
				{
					//j属于N[v]\X
					//v的邻居j不属于u的邻居，所以这次检查关系中u不会加入I
					flagI = -2;
				}
			}
		}
		if (flagX != -3 && flagI != -3)
		{
			return 0;//全部判断完成
		}
	}
	//运行到这里，证明flagX和flagI中至少有一个不为-3
	if (flagX == -3)
	{
		if (v_de_tmp_u < v_de_tmp_v)//u可能为v的子集
		{
			//printf("v_neightbor=%d  u_neightbor=%d\n", v_neightbor, u_neightbor);
			if (v_neightbor == u_neightbor)
			{
				flagX = 1;//X判断时，u的邻居全是v的邻居，代表N[u]\{N[S]并I}为N[v]的子集，需要将u加入X
			}
			else {//如果v_neightbor ！= u_neightbor，说明uv没有关系，flagX=-2
				flagX = -2;
			}
		}
		else {//v可能为u的子集
			if (v_de_tmp_u == v_de_tmp_v)  flagX = 0;
			else flagX = -1;//X判断时，没有出现v的邻居不是u的邻居的情况，代表N[v]\{N[S]并I}为N[u]的子集，需要将v加入X
		}
	}
	if (flagI == -3)
	{
		if (u_not_in_X < v_not_in_X)//N[u]\X可能为N[v]\X的子集
		{
			//printf("v_neightbor_I=%d  temp=%d\n", v_neightbor_I,temp);
			if (v_neightbor_I == u_not_in_X)
				flagI = -1;//I判断时，u的邻居全是v的邻居，代表N[u]\X为N[v]\X的子集，需要将v加入I
			else {//说明uv没有关系，flagI=-2
				flagI = -2;
			}
		}
		else {//N[v]\X可能为N[u]\X的子集
				flagI = 1;//I判断时，没有出现N[v]\X不属于N[u]\X的情况，代表N[v]\X为N[u]\X的子集，需要将u加入I
		}
	}

	return 0;
}

int new_check_belong2(int u, int v, int round)
{
	//重置flag
	flagX = -3, flagI = -3;
	//flagX=-2代表u和v不会加入入X中。flagI=-2代表u和v不会加入I中
	//计算N[u]\X和N[v]\X
	if (is_covered[u] || is_covered[v])  flagI = -2;//规定uv都必须是未被覆盖的顶点
	if (record_in[u] == 0 || record_in[v] == 0)  flagX = -2;//自己规定这样不检查
	if (record_in[u] == -1 && record_in[v] == -1)  flagI = -2;//因为uv不相连，所以record_in都是-1的情况下，检查I关系必为-2
	if (is_covered[u] == 0 && is_covered[v] == 0)  flagX = -2;//此时的v不可能是u的邻居，因此只有u和v有一个被覆盖前提下，uv的X关系检查才有用


	if (v_num > 500000) {
		if (v_degree_tmp[v] > 500 && v_degree_tmp[u] > 500) {
			flagX = -2;
		}
		if (v_degree_record[v] > 500 && v_degree_record[u] > 500) {
			flagI = -2;
		}
	}

	int u_not_in_X = v_degree_record[u];//N[u]\X
	int v_not_in_X = v_degree_record[v];//N[v]\X
	
	//计算N[u]\(N[S]并I)和N[v]\(N[S]并I),判二邻居时，这里的uv不可能有边相连，因此和一邻居的计算方法不一样
	int v_de_tmp_u = v_degree_tmp[u];//N(u)\(N[S]并I)
	int v_de_tmp_v = v_degree_tmp[v];
	if (is_covered[v] == 0)  v_de_tmp_v++;
	if (is_covered[u] == 0)  v_de_tmp_u++;

	if (need_check_I[u] == 0 && need_check_I[v] == 0)  flagI = -2;//两个都不用检查I关系
	if (need_check_X[u] == 0 && need_check_X[v] == 0)  flagX = -2;//两个都不用检查X关系
	/*
	if (round > 1)
	{
		if (need_check_I[u] == 0 && need_check_I[v] == 0)  flagI = -2;//两个都不用检查I关系
		else if (need_check_I[u] == 1 && need_check_I[v] == 1);//两个都需要检查，就不做任何处理
		else if (need_check_I[v] == 1)//只有v需要检查I关系
		{
			if (v_not_in_X > u_not_in_X)//这里v只需要在N[v]\X <= N[u]\X的情况下进行检查
			{
				flagI = -2;
			}
		}
		else if (need_check_I[u] == 1)//只有u需要检查I关系
		{
			if (u_not_in_X > v_not_in_X)//这里v只需要在N[u]\X <= N[v]\X的情况下进行检查
			{
				flagI = -2;
				//printf("u=%d v=%d  u_not_in_X=%d  v_not_in_X=%d\n", u, v, u_not_in_X, v_not_in_X);
			}
		}

		if (need_check_X[u] == 0 && need_check_X[v] == 0)  flagX = -2;//两个都不用检查X关系
		else if (need_check_X[u] == 1 && need_check_X[v] == 1);//两个都需要检查，就不做任何处理
		else if (need_check_X[v] == 1)//只有v需要检查X关系
		{
			if (v_de_tmp_v > v_de_tmp_u)//这里v只需要在N[v]\N[S] <= N[u]\N[S]的情况下进行检查
			{
				flagX = -2;
			}
		}
		else if (need_check_X[u] == 1)//只有u需要检查X关系
		{
			if (v_de_tmp_u > v_de_tmp_v)//这里v只需要在N[u]\N[S] <= N[v]\N[S]的情况下进行检查
			{
				flagX = -2;
			}
		}
	}*/
	
	if (((v_de_tmp_u < v_de_tmp_v) && record_in[u] == 1))// || (v_degree_tmp[u] == 0 && record_in[u] == 1)
		flagX = -2;
	if (((v_de_tmp_v < v_de_tmp_u) && record_in[v] == 1))// || (v_degree_tmp[v] == 0 && record_in[v] == 1)//若v被固定为1则v_degree_tmp[v]一定为0，
		flagX = -2;

	if (v_de_tmp_u <= v_de_tmp_v && is_covered[u] == 0)//u不属于N[v]\N[S]
		flagX = -2;
	if (v_de_tmp_v <= v_de_tmp_u && is_covered[v] == 0)//v不属于N[u]\N[S]
		flagX = -2;

	if (u_not_in_X <= v_not_in_X && record_in[u] == -1)  flagI = -2;//u不属于N[v]\X
	if (v_not_in_X <= u_not_in_X && record_in[v] == -1)  flagI = -2;//v不属于N[u]\X

	u_neightbor = v_degree_tmp[u];//在这里u一定被覆盖，所以N[u]=N(u),因为N[v]\N[S] > N[u]\N[S]且u未被覆盖,flagX直接为-2
	int v_neightbor = 0, v_neightbor_I = 0;
	//if (is_covered[v] == 0)  v_neightbor++;//因为判二邻居关系时，u的邻居里面肯定没有v，所以不需要这样
	//if (record_in[v] != 0)  v_neightbor_I++;//因为判二邻居关系时，u的邻居里面肯定没有v，所以不需要这样
	for (int i = 0, j; i < v_degree[v]; i++)
	{
		j = v_adj[v][i];
		//X: C(u) = N[u] \ (N[S] 并 I) 属于 N[v], add u to X
		if (flagX == -3)
		{
			if (v_de_tmp_u < v_de_tmp_v)//u可能为v的子集
			{
				if (is_covered[j] == 0 && neighbor_indicator[j] == true)
				{
					v_neightbor++;
					if (v_neightbor == u_neightbor)
						flagX = 1;//X判断时，u的邻居全是v的邻居，代表N[u]\{N[S]并I}为N[v]的子集，需要将u加入X
				}
				if (v_neightbor + v_degree[v] - i - 1 < u_neightbor)
				{
					flagX = -2;
				}
			}
			else {//N[v]可能为N[u]的子集
				if (is_covered[j] == 0 && neighbor_indicator[j] == false)
				{
					//flagX == -1代表还没有判断出u是否加入X
					//j属于N[v]\{N[S]并I}
					//v的邻居j不属于u的邻居，所以这次检查关系中u不会加入X
					flagX = -2;
				}
			}
		}

		//I: D(u) = N[u] \ X 属于 N[v], add v to I
		if (flagI == -3)
		{
			if (u_not_in_X < v_not_in_X)//N[u]\X可能为N[v]\X的子集
			{
				if (record_in[j] != 0 && neighbor_indicator[j] == true)
				{
					v_neightbor_I++;
					if (v_neightbor_I == u_not_in_X)
						flagI = -1;//I判断时，u的邻居全是v的邻居，代表N[u]\X为N[v]\X的子集，需要将v加入I
				}
				if (v_neightbor_I + v_degree[v] - i - 1 < u_not_in_X)
				{
					flagI = -2;
				}
			}
			else {//N[v]\X可能为N[u]\X的子集
				if (record_in[j] != 0 && neighbor_indicator[j] == false)
				{
					//j属于N[v]\X
					//v的邻居j不属于u的邻居，所以这次检查关系中u不会加入I
					flagI = -2;
				}
			}
		}
		if (flagX != -3 && flagI != -3)
		{
			//printf("return flagX = %d 111\n",flagX);
			return 0;//全部判断完成
		}
	}
	//运行到这里，证明flagX和flagI中至少有一个不为-3
	if (flagX == -3)
	{
		if (v_de_tmp_u < v_de_tmp_v)//u可能为v的子集
		{
			if (v_neightbor == u_neightbor)
				flagX = 1;//X判断时，u的邻居全是v的邻居，代表N[u]\{N[S]并I}为N[v]的子集，需要将u加入X
			else {//如果v_neightbor ！= u_neightbor，说明uv没有关系，flagX=-2
				flagX = -2;
			}
		}
		else {//v可能为u的子集
			//if (flagX == -3)
			if (v_de_tmp_u == v_de_tmp_v)  flagX = 0;
			else flagX = -1;//X判断时，没有出现v的邻居不是u的邻居的情况，代表N[v]\{N[S]并I}为N[u]的子集，需要将v加入X
		}
	}
	//printf("u_not_in_X=%d  v_not_in_X=%d\n", u_not_in_X, v_not_in_X);
	if (flagI == -3)
	{
		if (u_not_in_X < v_not_in_X)//N[u]\X可能为N[v]\X的子集
		{
			if (v_neightbor_I == u_not_in_X)
				flagI = -1;//I判断时，u的邻居全是v的邻居，代表N[u]\X为N[v]\X的子集，需要将v加入I
			else {//说明uv没有关系，flagI=-2
				flagI = -2;
			}
		}
		else {//N[v]\X可能为N[u]\X的子集
			//if (flagI == -3)
			flagI = 1;//I判断时，没有出现N[v]\X不属于N[u]\X的情况，代表N[v]\X为N[u]\X的子集，需要将u加入I
		}
	}
	return 0;
}

void add_to_N2(int u)
{
	if (pos_in_N2[u] == -1)
	{
		N2[num_N2] = u;
		pos_in_N2[u] = num_N2;
		num_N2++;
	}
}

void find_N2_rule5(int u, int v)
{
	int i, j, k, x;
	for (i = 0; i < v_degree[u]; i++)
	{
		k = v_adj[u][i];
		//if (record_in[k] != -1)  continue;
		if (k == v || is_covered[k])  continue;//v3和v4要求未被覆盖
		for (j = 0; j < v_degree[k]; j++)
		{
			x = v_adj[k][j];
			if (neighbor_indicator[x] == 0 && record_in[x] == -1)  break;//因为v3和v4不会被覆盖，因此record_in[x]不会等于1
		}
		if (j == v_degree[k])
		{
			add_to_N2(k);
		}
	}

	for (i = 0; i < v_degree[v]; i++)
	{
		k = v_adj[v][i];
		//if (record_in[k] != -1 || N2[k])  continue;
		if (pos_in_N2[k] != -1 || k == u || is_covered[k])  continue;
		for (j = 0; j < v_degree[k]; j++)
		{
			x = v_adj[k][j];
			if (neighbor_indicator[x] == 0 && record_in[x] == -1)  break;
		}
		if (j == v_degree[k])
		{
			add_to_N2(k);
		}
	}
}

void find_N2(int u, int v)
{
	int i, j, k, x;
	for (i = 0; i < v_degree[u]; i++)
	{
		k = v_adj[u][i];
		//if (record_in[k] != -1)  continue;
		if (k == v)  continue;
		for (j = 0; j < v_degree[k]; j++)
		{
			x = v_adj[k][j];
			if (neighbor_indicator[x] == 0 && is_covered[x] == 0)  break;
		}
		if (j == v_degree[k])
		{
			add_to_N2(k);
		}
	}

	for (i = 0; i < v_degree[v]; i++)
	{
		k = v_adj[v][i];
		//if (record_in[k] != -1 || N2[k])  continue;
		if (pos_in_N2[k] != -1 || k == u)  continue;
		for (j = 0; j < v_degree[k]; j++)
		{
			x = v_adj[k][j];
			if (neighbor_indicator[x] == 0 && is_covered[x] == 0)  break;
		}
		if (j == v_degree[k])
		{
			add_to_N2(k);
		}
	}
}

void new_reduction()
{
	//printf("new reduction\n");
	int k, i, j, flag;
	int cccc = 0;
	int degree_bigger_than_2_index = 0;
	bool useI = false;
	useI = true;
	bool userTN = false;
	userTN = true;
	bool use_rule4 = false;
	use_rule4 = true;
	bool only_round1 = false;
	//only_round1 = true;

	if (e_num > 1888888)
	{
		is_large = true;
		userTN = false;
	}
	if (v_num > 1000)
	{
		rule4_is_large = true;
	}
	use_rule4 = false;
	countccc = 0;
	int countttt = 0;
	bool backflag = true;
	int round = 1;
	double time_check = 0, time_update = 0;
	double t_strat = TimeElapsed();
	while (num_working_vertice > 0)
	{
		int l1fix1 = 0, l2fix1 = 0;
		int l1fix0 = 0, l2fix0 = 0;
		int l1fixI = 0, l2fixI = 0;
		for (int kk = 0; kk < num_working_vertice; kk++)
		{
			fix_time_stamp++;
			int u = working_vertice[kk];
			int v;

			if (v_degree[u] == 0)//只有自己
			{
				if (round > 1)  continue;
				add_to_must_in(u);
				degree_bigger_than_2_index++;
			}
			else if (v_degree[u] == 1)
			{
				if (round > 1)  continue;
				int n = v_adj[u][0];
				if (v_degree[n] == 1) {
					if (v_adj[n][0] != u)
						cout << "v_adj[n][0] != u" << endl;
					if (u > n) {
						if (record_in[n] == 0)
							cout << "record_in[v1] == 0 1 " << endl;
						record_in[u] = 0;
						add_to_must_in(n);
						add_to_new_fix(n);
					}
				}
				else {
					record_in[u] = 0;//若u的度为1，则u固定为0
					if (record_in[n] != 1) {
						if (record_in[n] == 0)
						{
							cout << "record_in[n] == 0 2 " << file << endl;
							exit(0);
						}
						add_to_must_in(n);
						update_v_degree_tmp(n);
						add_to_new_fix(n);
					}
				}
				degree_bigger_than_2_index++;
				time_out[u] = fix_time_stamp;
			}
			else {

				if (is_large && record_in[u] == 0)//大图专用
				{
					continue;
				}

				//标记邻居***
				u_neightbor = 0;
				neighbor_indicator[u] = true;
				for (i = 0; i < v_degree[u]; i++)
				{
					j = v_adj[u][i];
					neighbor_indicator[j] = true;
				}

				bool ff = false;
				for (i = 0; i < v_degree[u]; i++)
				{
					v = v_adj[u][i];//u的度很小 v的度很大
					
					//if (!addable[v] && round == 1)
					if (!addable[v])
					{
						if (v_degree[u] < v_degree[v])
							continue;
						else if (v_degree[u] == v_degree[v] && u > v)
							continue;
					}
					if (record_in[u] == 1 && record_in[v] == 1)  continue;

					if (is_large)//大图专用
					{
						if (record_in[u] == 0)  break;
						if (record_in[v] == 0)  continue;
						if (need_check_X[u] == 0 && need_check_X[v] == 0)  continue;
						int v_de_tmp_u = v_degree_tmp[u];//N(u)\(N[S]并I)
						int v_de_tmp_v = v_degree_tmp[v];
						v_de_tmp_u--;
						v_de_tmp_v--;
						if (is_covered[v] == 1)  v_de_tmp_u++;
						if (is_covered[u] == 1)  v_de_tmp_v++;
						if ((v_de_tmp_u < v_de_tmp_v && record_in[u] == 1))// || (v_degree_tmp[u] == 0 && record_in[u] == 1)
							continue;
						if ((v_de_tmp_u > v_de_tmp_v && record_in[v] == 1))// || (v_degree_tmp[v] == 0 && record_in[v] == 1)//若v被固定为1则v_degree_tmp[v]一定为0，
							continue;
					}

					neighbor_indicator2[v] = true;//后面二邻居不需要再判uv关系
					//t_strat = TimeElapsed();
					new_check_belong1(u, v, round);
					//time_check += TimeElapsed() - t_strat;
					
					if (flagX == 1)
					{
						if (record_in[u] == 1)
							cout << file << " record_in[u] == 1 " << u << endl;                                                                                                                                                                                                                            
						if (record_in[u] == 0)  continue;
						record_in[u] = 0;
						if (!is_large)  add_to_new_X(u);
						l1fix0++;
						v_degree_record[u]--;
						if (v_degree_record[u] == 1 && is_covered[u] == 0) {
							fix_only_neighbor1(u);
							l1fix1++;
						}
						for (k = 0; k < v_degree[u]; k++)
						{
							int verr = v_adj[u][k];
							v_degree_record[verr]--;
							if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {
								fix_only_neighbor2(verr);
								l1fix1++;
							}
						}
					}
					else if (flagX == -1)
					{//需要考虑一种情况，若二者的v_degree_tmp = 0且有一个固定为1了另一个固定为0
						if (record_in[v] == 1)
							cout << file << " record_in[v] == 1 " << v << endl;
						if (record_in[v] == 0)  continue;
						record_in[v] = 0;
						if (!is_large)  add_to_new_X(v);
						l1fix0++;
						v_degree_record[v]--;
						if (v_degree_record[v] == 1 && is_covered[v] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
							fix_only_neighbor1(v);
							l1fix1++;
						}
						for (k = 0; k < v_degree[v]; k++)
						{
							int verr = v_adj[v][k];
							v_degree_record[verr]--;
							if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
								fix_only_neighbor2(verr);
								l1fix1++;
							}
						}
					}
					else if (flagX == 0) {//需要考虑一种情况，若二者的v_degree_tmp = 0且有一个固定为1了另一个固定为0
						if (record_in[v] == 1) {
							v = u;
						}
						if (record_in[v] == 0)  continue;
						record_in[v] = 0;
						if (!is_large)  add_to_new_X(v);
						l1fix0++;
						equal_Reduction.push_back(v);
						v_degree_record[v]--;
						if (v_degree_record[v] == 1 && is_covered[v] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
							fix_only_neighbor1(v);
							l1fix1++;
						}
						for (k = 0; k < v_degree[v]; k++)
						{
							int verr = v_adj[v][k];
							v_degree_record[verr]--;
							if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
								fix_only_neighbor2(verr);
								l1fix1++;
							}
						}
					}
					
					//Rule 3 ：Ignorable Vertices
					if (flagI == 1)//N[v]\X为N[u]\X的子集，将u加入I
					{
						if (is_covered[u])//在上面将v加入X中，可能会导致is_covered[u]=1
							continue;
						//因为u不再需要被覆盖，所以要更新u邻居的v_degree_tmp
						for (int j = 0; j < v_degree[u]; j++) {
							int v11 = v_adj[u][j];
							v_degree_tmp[v11]--;//v_degree_tmp[v1]需要跟v1的邻居中is_coverd为0的顶点个数相等。
						}
						add_to_new_fix_neightbor(u);
						is_covered[u] = 1;
						new_covered[u] = true;
						l1fixI++;
					}
					else if (flagI == -1)////N[u]\X为N[v]\X的子集，将v加入I
					{
						if (is_covered[v])//在上面将u加入X中，可能会导致is_covered[v]=1
							continue;
						for (int j = 0; j < v_degree[v]; j++) {
							int v11 = v_adj[v][j];
							v_degree_tmp[v11]--;//v_degree_tmp[v1]需要跟v1的邻居中is_coverd为0的顶点个数相等。
						}
						add_to_new_fix_neightbor(v);
						is_covered[v] = 1;
						new_covered[v] = true;
						l1fixI++;
					}
				}
			
				//-----------------------------双邻居---------------------------
				if (userTN && !(is_covered[u] == 0 && record_in[u] == -1))
				{
					int w = -1;
					for (int _i = 0; _i < v_degree[u]; _i++)
					{
						w = v_adj[u][_i];
						for (int i = 0; i < v_degree[w]; i++)
						{
							v = v_adj[w][i];
							if (neighbor_indicator2[v] == true)//之前已经判断过uv的关系了
								continue;
							if (v == u)  continue;

							//if (!addable[v] && round == 1)
							/*
							if (!addable[v])
							{
								if (v_degree[u] < v_degree[v])
									continue;
								else if (v_degree[u] == v_degree[v] && u > v)
									continue;
							}*/
							if (record_in[u] == 1 && record_in[v] == 1)
								continue;
							neighbor_indicator2[v] = true;
							//因为上面u和u的邻居已经判断过关系了，所以这里的二邻居v不可能是u的邻居
							new_check_belong2(u, v, round);

							if (flagX == 1)
							{
								if (record_in[u] == 1)
									cout << file << " record_in[u] == 1 " << u << endl;
								if (record_in[u] == 0)  continue;
								record_in[u] = 0;
								add_to_new_X(u);
								l2fix0++;
								v_degree_record[u]--;
								if (v_degree_record[u] == 1 && is_covered[u] == 0) {
									fix_only_neighbor1(u);
									l2fix1++;
								}
								for (k = 0; k < v_degree[u]; k++)
								{
									int verr = v_adj[u][k];
									v_degree_record[verr]--;
									if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {
										fix_only_neighbor2(verr);
										l2fix1++;
									}
								}
							}
							else if (flagX == -1)
							{//需要考虑一种情况，若二者的v_degree_tmp = 0且有一个固定为1了另一个固定为0
								if (record_in[v] == 1)
									cout << file << " record_in[v] == 1 " << v << endl;
								if (record_in[v] == 0)  continue;
								record_in[v] = 0;
								add_to_new_X(v);
								l2fix0++;
								v_degree_record[v]--;
								if (v_degree_record[v] == 1 && is_covered[v] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
									fix_only_neighbor1(v);
									l2fix1++;
								}
								for (k = 0; k < v_degree[v]; k++)
								{
									int verr = v_adj[v][k];
									v_degree_record[verr]--;
									if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
										fix_only_neighbor2(verr);
										l2fix1++;
									}
								}
							}
							else if (flagX == 0) {//需要考虑一种情况，若二者的v_degree_tmp = 0且有一个固定为1了另一个固定为0
								if (record_in[v] == 1) {
									v = u;
								}
								if (record_in[v] == 0)  continue;
								record_in[v] = 0;
								add_to_new_X(v);
								l2fix0++;
								equal_Reduction.push_back(v);
								v_degree_record[v]--;
								if (v_degree_record[v] == 1 && is_covered[v] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
									fix_only_neighbor1(v);
									l2fix1++;
								}
								for (k = 0; k < v_degree[v]; k++)
								{
									int verr = v_adj[v][k];
									v_degree_record[verr]--;
									if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
										fix_only_neighbor2(verr);
										l2fix1++;
									}
								}
							}
							//Rule 3 ：Ignorable Vertices
							if (flagI == 1)//N[v]\X为N[u]\X的子集，将u加入I
							{
								if (is_covered[u])//在上面将v加入X中，可能会导致is_covered[u]=1
									continue;
								//因为u不再需要被覆盖，所以要更新u邻居的v_degree_tmp
								for (int j = 0; j < v_degree[u]; j++) {
									int v11 = v_adj[u][j];
									v_degree_tmp[v11]--;//v_degree_tmp[v1]需要跟v1的邻居中is_coverd为0的顶点个数相等。
								}
								add_to_new_fix_neightbor(u);
								is_covered[u] = 1;
								new_covered[u] = true;
								l2fixI++;
							}
							else if (flagI == -1)////N[u]\X为N[v]\X的子集，将v加入I
							{
								if (is_covered[v])//在上面将u加入X中，可能会导致is_covered[v]=1
									continue;
								//因为u不再需要被覆盖，所以要更新u邻居的v_degree_tmp
								for (int j = 0; j < v_degree[v]; j++) {
									int v11 = v_adj[v][j];
									v_degree_tmp[v11]--;//v_degree_tmp[v1]需要跟v1的邻居中is_coverd为0的顶点个数相等。
								}
								add_to_new_fix_neightbor(v);
								is_covered[v] = 1;
								new_covered[v] = true;
								l2fixI++;
							}
						}
					}
					for (int _i = 0, w; _i < v_degree[u]; _i++)//清除双邻居搜索标记
					{
						w = v_adj[u][_i];
						for (int i = 0; i < v_degree[w]; i++)
						{
							v = v_adj[w][i];
							neighbor_indicator2[v] = false;
						}
					}
				}

				//清除邻居标记***
				neighbor_indicator[u] = false;
				for (i = 0; i < v_degree[u]; i++)
				{
					j = v_adj[u][i];
					neighbor_indicator[j] = false;
					neighbor_indicator2[j] = false;
				}
			}
			
			//need_check_I[u] = false;
			//need_check_X[u] = false;
			//addable[u] = true;
		}

		if(only_round1)  return ;//只看第一轮结果

		for (int kk = 0, x; kk < num_working_vertice; kk++)
		{
			x = working_vertice[kk];
			addable[x] = true;
			need_check_I[x] = false;
			need_check_X[x] = false;
		}
		//t_strat = TimeElapsed();
		update_next_working_vertice();
		//time_update += TimeElapsed() - t_strat;

		round++;//check轮数

		countccc++;

		num_working_vertice = num_next_working_vertice;
		for (int kk = 0; kk < num_working_vertice; kk++)
		{
			working_vertice[kk] = next_working_vertice[kk];
			//addable[working_vertice[kk]] = false;
			//pos_in_working_vertice[working_vertice[kk]] = kk;
		}
		num_next_working_vertice = 0;
	}

	//cout << time_check <<" "<<time_update << "\n";

	//rule5
	int coffe = v_num / 10000;


	int max_count = 10000, temp_count = 0;
	if (v_num > 10000)
		max_count = max_count / coffe;

	temp_count = 0;
	if (use_rule4)
	{
		for (int v = 0; v < v_num; v++)
		{
			working_vertice[v] = v_d_vec[v].v;
		}
		num_working_vertice = v_num;
		for (int kk = 0; kk < num_working_vertice; kk++)//只遍历一遍
		{
			int u = working_vertice[kk];
			int v, flag;
			if (record_in[u] == 0)  continue;//u不会被加入解集，因此略过
			//标记邻居***
			neighbor_indicator[u] = 1;
			for (i = 0; i < v_degree[u]; i++)
			{
				j = v_adj[u][i];
				neighbor_indicator[j] = 1;
			}
			for (i = 0; i < v_degree[u]; i++)
			{
				if (rule4_is_large && temp_count > max_count)
				//if (temp_count > max_count)
				{
					break;
				}
				v = v_adj[u][i];//u的度很小 v的度很大
				if (record_in[v] == 0 || record_in[u] == 0)  continue;
				temp_count++;//检查对数加一
				
				neighbor_indicator[v] += 1;
				for (j = 0; j < v_degree[v]; j++)
				{
					k = v_adj[v][j];
					neighbor_indicator[k] += 1;
				}

				find_N2(u, v);//寻找uv的属于N2和N3,且record_in=-1的邻居

				if (num_N2 > 1)
				{
					//接下来就是num_N2>1的情况了
					for (j = 0; j < num_N2; j++)
					{
						k = N2[j];
						if (is_covered[k])  continue;

						flag = true;
						for (int ii = 0, _kk; ii < v_degree[k]; ii++)
						{
							_kk = v_adj[k][ii];
							if (record_in[_kk] != 0 && pos_in_N2[_kk] == -1 && _kk != u && _kk != v)
							{
								flag = false;
								break;
							}
						}

						if (flag)
						{
							is_N2_neightbor[k] = true;
							for (int ii = 0, _kk; ii < v_degree[k]; ii++)
							{
								_kk = v_adj[k][ii];
								is_N2_neightbor[_kk] = true;
							}

							//
							for (int ii = 0, _kk; ii < num_N2; ii++)
							{
								_kk = N2[ii];
								if (record_in[_kk] == -1 && is_N2_neightbor[_kk] == false)
								{
									record_in[_kk] = 0;
									v_degree_record[_kk]--;
									if (v_degree_record[_kk] == 1 && is_covered[_kk] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
										fix_only_neighbor1(_kk);
									}
									for (int temp = 0; temp < v_degree[_kk]; temp++)
									{
										int verr = v_adj[_kk][temp];
										v_degree_record[verr]--;
										if (v_degree_record[verr] == 1 && is_covered[verr] == 0) {//若verr闭邻居只有一个顶点v'未被record为0。则v'固定为1
											fix_only_neighbor2(verr);
										}
									}
								}
							}

							is_N2_neightbor[k] = false;
							for (int ii = 0, _kk; ii < v_degree[k]; ii++)
							{
								_kk = v_adj[k][ii];
								is_N2_neightbor[_kk] = false;
							}
						}
					}
				}

				//清除标记
				neighbor_indicator[v] -= 1;
				for (j = 0; j < v_degree[v]; j++)
				{
					k = v_adj[v][j];
					neighbor_indicator[k] -= 1;
				}

				for (int ii = 0, _kk; ii < num_N2; ii++)
				{
					_kk = N2[ii];
					pos_in_N2[_kk] = -1;
				}
				num_N2 = 0;
			}
			neighbor_indicator[u] -= 1;
			for (i = 0; i < v_degree[u]; i++)
			{
				j = v_adj[u][i];
				neighbor_indicator[j] -= 1;
			}

			if (rule4_is_large && temp_count > max_count)
			//if (temp_count > max_count)
			{
				break;
			}
		}
	}

	
	delete[] working_vertice;
	delete[] next_working_vertice;
	delete[] addable;
	delete[] neighbor_indicator;
	delete[] new_fix;
	delete[] pos_in_new_fix;
	delete[] new_fix_neightbor;
	delete[] pos_in_new_fix_neightbor;
	delete[] time_out;
	delete[] time_need_in;
	delete[] fix_time;
}

int BuildInstance(string filename) {
  string tmp;
  string sign;
  int v, e;
  int v1, v2;

  ifstream infile(filename);
  if (!infile) {
    return 1;
  }

  while (infile.peek() != 'p') {
    getline(infile, tmp);
  }

  infile >> sign >> tmp >> v_num >> e_num;

  edge = new Edge[e_num];
  score = new int[v_num + 1];
  // sub_score = new int[v_num + 1];
  time_stamp = new llong[v_num + 1];
  v_adj = new int *[v_num + 1];
  v_degree = new int[v_num + 1];
  v_in_c = new int[v_num + 1];
  candidate = new int[v_num + 1];
  index_in_candidate = new int[v_num + 1];
  //conf_change = new int[v_num + 1];
  //conf_count = new int[v_num + 1];
  dominated = new int[v_num + 1];
  undom_stack = new int[v_num + 1];
  index_in_undom_stack = new int[v_num + 1];
  fill_n(index_in_undom_stack, v_num + 1, -1);

  v_fixed = new int[v_num + 1];
  temp_store_reduction = new int[v_num + 1];


  my_heap = new int[v_num + 1];
  pos_in_my_heap = new int[v_num + 1];
  fill_n(pos_in_my_heap, v_num + 1, -2);
  //cover_set = new int[v_num + 1];
  //cover_set2 = new int[v_num + 1];
  my_heap_count = 0;
  store_ver = new int[v_num + 1];

  unlock_propagate_ver_num = v_num;
  unlock_propagate_edge_num = e_num;
  propagate_degree = new int[v_num + 1];
  propagate_valid = new int[v_num + 1];
  index_in_propagate_valid = new int[v_num + 1];

  fill_n(v_degree, v_num + 1, 0);
  fill_n(v_in_c, v_num + 1, 0);
  fill_n(score, v_num + 1, 0);
  // fill_n(sub_score, v_num + 1, 0);
  //fill_n(conf_change, v_num + 1, 0);
  fill_n(time_stamp, v_num + 1, 0);
  fill_n(dominated, v_num + 1, 0);
  fill_n(v_fixed, v_num + 1, 0);
  
  //inference
  must_in = new int[v_num + 1];
  pos_in_must_in = new int[v_num + 5];
  fill_n(pos_in_must_in, v_num + 1, -1);
  record_in = new int[v_num + 1];
  fill_n(record_in, v_num + 1, 0);
  I = new int[v_num + 1];
  fill_n(I, v_num + 1, 0);

  temp_flag = new int[v_num + 1];
  fill_n(temp_flag, v_num + 1, 0);
  temp_flag_stack = new int[v_num + 1];

  index_in_zero_stack = new int[v_num + 1];
  zero_stack = new int[v_num + 1];
  memset(index_in_zero_stack, -1, sizeof(int) * (v_num + 1));
  memset(zero_stack, 0, sizeof(int) * (v_num + 1));
  zero_stack_fill_pointer = 0;

  //reduction
  addable = new bool[v_num + 1];
  fill_n(addable, v_num + 1, 0);
  neighbor_indicator = new int[v_num + 1];
  fill_n(neighbor_indicator, v_num + 1, 0);
  neighbor_indicator2= new int[v_num + 1];
  fill_n(neighbor_indicator2, v_num + 1, 0);
  new_fix = new int[v_num + 1];
  new_fix_neightbor = new int[v_num + 1];
  fill_n(new_fix, v_num + 1, 0);
  fill_n(new_fix_neightbor, v_num + 1, 0);
  pos_in_new_fix = new int[v_num + 1];
  pos_in_new_fix_neightbor = new int[v_num + 1];
  fill_n(pos_in_new_fix, v_num + 1, -1);
  fill_n(pos_in_new_fix_neightbor, v_num + 1, -1);
  num_new_fix = 0;
  num_new_fix_neightbor = 0;
  fix_time_stamp = 0;//fix时间戳置0
  fix_time = new long long[v_num + 1];
  time_out = new long long[v_num + 1];
  time_need_in = new long long[v_num + 1];
  fill_n(fix_time, v_num + 1, 0);
  fill_n(time_out, v_num + 1, 0);
  fill_n(time_need_in, v_num + 1, 0);
  working_vertice = new int[v_num + 1];
  pos_in_working_vertice = new int[v_num + 1];
  next_working_vertice = new int[v_num + 1];
  num_working_vertice = 0;
  num_next_working_vertice = 0;
  need_check_I = new bool[v_num + 1];
  need_check_X = new bool[v_num + 1];
  num_new_X = 0;
  new_X = new int[v_num + 1];//存储新加入X的顶点的数组
  new_covered = new bool[v_num + 1];//新被covered的顶点
  fill_n(new_covered, v_num + 1, 0);

  //reduction rule 4
  N2 = new int[v_num + 1];
  pos_in_N2 = new int[v_num + 1];
  fill_n(N2, v_num + 1, 0);
  fill_n(pos_in_N2, v_num + 1, -1);
  is_N2_neightbor = new bool[v_num + 1];
  fill_n(is_N2_neightbor, v_num + 1, 0);

  while (infile.peek() != 'e') {
	  getline(infile, tmp);
  }

  for (e = 0; e < e_num; e++) {
    infile >> tmp >> v1 >> v2;
    v_degree[v1]++;
    v_degree[v2]++;

    edge[e].v1 = v1;
    edge[e].v2 = v2;
  }
  infile.close();

  v_adj[0] = 0;
  for (v = 1; v < v_num + 1; v++) {
    v_adj[v] = new int[v_degree[v]];
	push_propagate_valid(v);
	propagate_degree[v] = v_degree[v];
  }

  v_degree_tmp = new int[v_num + 5];//顶点v的邻居中有多少未被fix_rule覆盖的顶点
  fill_n(v_degree_tmp, v_num + 5, 0);
  v_degree_record = new int[v_num + 5];//闭邻居有多少未赋值0的
  fill_n(v_degree_record, v_num + 5, 0);

  for (e = 0; e < e_num; e++) {
    v1 = edge[e].v1;
    v2 = edge[e].v2;

    v_adj[v1][v_degree_tmp[v1]] = v2;
    v_adj[v2][v_degree_tmp[v2]] = v1;

    v_degree_tmp[v1]++;
    v_degree_tmp[v2]++;
  }
  

  is_covered = new int[v_num + 1];


  for (v = 1; v < v_num + 1; v++) {

	  if (v_degree[v] > Max_degree) {
		  Max_degree = v_degree[v];
	  }
  }

  //vector<two_ele> v_d_vec;
  for (v = 1; v < v_num + 1; v++) {
	  v_degree_record[v] = v_degree[v] + 1;
	  is_covered[v] = 0;//所有顶点初始都是uncover的
	  two_ele te;
	  te.v = v;
	  te.degree = v_degree[v];
	  v_d_vec.push_back(te);
  }
  //sort(v_d_vec.begin(), v_d_vec.end(), CompareBySecond);
  
  sort(v_d_vec.begin(), v_d_vec.end(), CompareBySecond2);
  
  for (int i = 0; i < v_d_vec.size() - 1; i++) {
	  if (v_d_vec[i].degree < v_d_vec[i + 1].degree) {
		  cout << " (v_degree[i] < v_degree[i + 1]) " << file << " " << v_degree[v_d_vec[i].v] << " " << v_degree[v_d_vec[i + 1].v] << endl;
	  }
  }

  for (int v = 0; v < v_num; v++)
  {
	  working_vertice[v] = v_d_vec[v].v;
	  pos_in_working_vertice[v_d_vec[v].v] = v;
	  need_check_I[v_d_vec[v].v] = true;
	  need_check_X[v_d_vec[v].v] = true;
  }
  num_working_vertice = v_num;

  delete[] edge;
  int u, k, i, j, flag;
  for (i = 1; i <= v_num; i++)  record_in[i] = -1;//fix初始化
  int cccc = 0;

  start = chrono::steady_clock::now();
  int degree_bigger_than_2_index = 0;

  new_reduction();
  
  int xx;
  int countttt = 0, num_fix0 = 0, num_fix1 = 0, I_num = 0;
  for (int i = 1; i <= v_num; i++) {
	  if (record_in[i] == 0 || record_in[i] == 1)
		  countttt++;
	  if (record_in[i] == 0)
	  {
		  num_fix0++;
	  }
	  if (record_in[i] == 1)
	  {
		  num_fix1++;
	  }
  }
  double t_redu = TimeElapsed();
//  cout << filename << " " << num_fix1 <<" "<< num_fix0 <<" " << num_fix0+num_fix1 <<" "<<countccc << "  " << t_redu << endl;
//  exit(0);

  delete[] v_degree_tmp;
  delete[] v_degree_record;
  delete[] is_covered;
  return 0;
}

void FreeMemory() {
  int v;
  for (v = 0; v < v_num + 1; v++) {
    delete[] v_adj[v];
  }

  delete[] dominated;
  delete[] undom_stack;
  delete[] index_in_undom_stack;
  delete[] index_in_candidate;
  delete[] candidate;
  delete[] v_in_c;
  delete[] v_degree;
  delete[] v_threshold;
  delete[] v_adj;
  delete[] time_stamp;
  delete[] score;
  delete[] my_heap;
  delete[] pos_in_my_heap;
  delete[] record_in;

  if (cover_set_flag)
  {
	  delete[] cover_set;
	  delete[] cover_set2;
  }
}

inline void Undom(int v) {
  index_in_undom_stack[v] = undom_stack_fill_pointer;
  Push(v, undom_stack);
}

inline void Dom(int v) {
  int index, last_undom_vertex;

  last_undom_vertex = Pop(undom_stack);
  index = index_in_undom_stack[v];
  undom_stack[index] = last_undom_vertex;
  index_in_undom_stack[last_undom_vertex] = index;
  index_in_undom_stack[v] = -1;
}

void init_increase_dominate(int v, int added_vertex) {
  if (dominated[v] == 0) {
    --score[v];
    int pos = pos_in_my_heap[v];
	if (pos != -2) {
		my_heap_remove(pos);
		my_heap_insert(v);
	}
    // all make but one break
    for (int i = 0; i < v_degree[v]; ++i) {
      int u = v_adj[v][i];
      --score[u];
      int pos = pos_in_my_heap[u];
	  if (pos != -2) {
		  my_heap_remove(pos);
		  my_heap_insert(u);
	  }
    }
	//cover_set***
	if (cover_set_flag)
	{
		cover_set[v] = added_vertex;
	}
    Dom(v);
  } else if (dominated[v] == 1) {
    // only two vertex need to update break;
	  //cover_set***
	  if (cover_set_flag)
	  {
		  int cur_v = cover_set[v];
		  ++score[cur_v];
		  cover_set2[v] = added_vertex;
	  }
	  else {
		  if (v_in_c[v] == 1) {
			  ++score[v];
		  }
		  // TODO: speed up
		  for (int i = 0; i < v_degree[v]; ++i) {
			  int u = v_adj[v][i];
			  if (v_in_c[u]) {
				  ++score[u];
			  }
		  }
	  }
	
  }
  ++dominated[v];
}

void init_decrease_dominate(int v, int remove_vertex) {
	if (dominated[v] == 1) {
		// all score of neighbours are make
		++score[v];
		for (int i = 0; i < v_degree[v]; ++i) {
			int u = v_adj[v][i];
			++score[u];
		}
		Undom(v);
	}
	else if (dominated[v] == 2) {
		// only need to update one break
		//cover_set***
		if (cover_set_flag)
		{
			if (cover_set[v] == remove_vertex) {
				cover_set[v] = cover_set2[v];
			}
			score[cover_set[v]]--;
		}
		else {
			if (v_in_c[v]) {
				--score[v];
			}
			else {
				// TODO: speed up this procedure
				for (int i = 0; i < v_degree[v]; ++i) {
					int u = v_adj[v][i];
					if (v_in_c[u]) {
						--score[u];
						//break;//*****
					}
				}
			}
		}
	}
	else if (dominated[v] == 3 && cover_set_flag) {
		int flag = 0;
		if (v_in_c[v] == 1) {
			cover_set[v] = v;
			flag = 1;
		}
		for (int i = 0; i < v_degree[v]; ++i) {
			int u = v_adj[v][i];
			if (v_in_c[u] == 1) {
				if (flag == 0)
				{
					cover_set[v] = u;
					flag = 1;
				}
				else
				{
					cover_set2[v] = u;
					break;
				}
			}
		}
	}
	--dominated[v];
}


void RemoveInit(int v) {
	v_in_c[v] = 0;
	c_size--;
	int last_candidate_v = candidate[--candidate_size];
	int index = index_in_candidate[v];
	candidate[index] = last_candidate_v;
	index_in_candidate[last_candidate_v] = index;
	index_in_candidate[v] = -1;
	score[v] = 0;
	init_decrease_dominate(v, v);
	for (int i = 0; i < v_degree[v]; ++i) {
		int u = v_adj[v][i];
		init_decrease_dominate(u, v);
	}
	if (zero_stack_fill_pointer < 0)
		cout << "zero_stack_fill_pointer < 0 " << file << " " << seed << endl;
}

void removeRedundant() {
  int v;
  while (zero_stack_fill_pointer > 0) {
	  int zero_index = rand() % zero_stack_fill_pointer;
	  v = zero_stack[zero_index];
	  Remove(v);
  }
}


bool test_score() {
	vector<int> correct_score(v_num + 1, 0);
	vector<int> correct_dominate(v_num + 1, 0);
	for (int v = 1; v < v_num + 1; ++v) {
		if (v_in_c[v]) {
			correct_dominate[v]++;
			for (int i = 0; i < v_degree[v]; ++i) {
				int u = v_adj[v][i];
				correct_dominate[u]++;
			}
		}
	}
	for (int v = 1; v < v_num + 1; ++v) {
		if (correct_dominate[v] == 0) {
			correct_score[v]++;
			for (int i = 0; i < v_degree[v]; ++i) {
				int u = v_adj[v][i];
				correct_score[u]++;
			}
		}
		else if (correct_dominate[v] == 1) {
			if (v_in_c[v]) {
				correct_score[v]--;
			}
			else {
				for (int i = 0; i < v_degree[v]; ++i) {
					int u = v_adj[v][i];
					if (v_in_c[u]) {
						correct_score[u]--;
					}
				}
			}
		}
	}
	bool right = true;
	for (int v = 1; v < v_num + 1; ++v) {
		if (correct_dominate[v] != dominated[v]) {
			std::cout << "dominated wrong: " << v << " " << file << std::endl;
			std::cout << "correct: " << correct_dominate[v]
				<< ", wrong: " << dominated[v] << std::endl;
			right = false;
		}
		if (correct_score[v] != score[v]) {
			std::cout << "score wrong: " << v << std::endl;
			std::cout << "correct: " << correct_score[v] << ", wrong: " << score[v]
				<< std::endl;
			std::cout << v_in_c[v] << std::endl;
			right = false;
		}
	}
	if (!right)
		cout << endl;

	right = true;
	int undom_num = 0;
	for (int v = 1; v < v_num + 1; ++v) {
		if (dominated[v] == 0) {
			undom_num++;
			bool okay = false;
			for (int i = 0; i < undom_stack_fill_pointer; ++i) {
				if (undom_stack[i] == v) {
					okay = true;
					break;
				}
			}
			if (!okay) {
				std::cout << "not in undom_stack: " << v << std::endl;
				right = false;
			}
		}
	}
	if (undom_num != undom_stack_fill_pointer) {
		std::cout << "undom_num wrong!" << std::endl;
	}
	if (!right)
		abort();

	return true;
}

void init_removeRedundant() {
	int v;
	for (int i = 0; i < candidate_size; ++i) {
		v = candidate[i];
		if (score[v] == 0) {
			RemoveInit(v);
			i--;
		}
	}
}

void increase_dominate(int v, int added_vertex) {
  if (dominated[v] == 0) {//开始的状态是必须要其一个邻居v才能保持支配的。后来的状态是不一定需要这个邻居v就可以保持支配
    --score[v];//加入一个顶点使v由未支配变为支配。若v或v的邻居u不在当前解中，加入u新增的覆盖点减少；若v或v的邻居u在当前解中，与支配v前相比移除u后objectivefunction减少量加1
	for (int i = 0; i < v_degree[v]; ++i) {
      int u = v_adj[v][i];
      --score[u];//若v_in_c[u] == 1则score[u]一定不肯能等于0
    }
    Dom(v);
	//cover_set***
	if (cover_set_flag)
	{
		cover_set[v] = added_vertex;
	}
  } else if (dominated[v] == 1) {
	  if (cover_set_flag)
	  {
		  int cur_v = cover_set[v];
		  ++score[cur_v];
		  if (score[cur_v] == 0 && index_in_zero_stack[cur_v] == -1 && record_in[cur_v] != 1) {
			  add_zero_vertex(cur_v);
		  }
		  cover_set2[v] = added_vertex;
	  }
	  else
	  {
		  int cur_v;
		  if (v_in_c[v] == 1) {
			  ++score[v];
			  cur_v = v;
		  }
		  else {
			  // TODO: speed up
			  for (int i = 0; i < v_degree[v]; ++i) {
				  int u = v_adj[v][i];
				  if (v_in_c[u]) {
					  ++score[u];
					  cur_v = u;
					  break;
				  }
			  }
		  }
		  if (score[cur_v] == 0 && index_in_zero_stack[cur_v] == -1 && record_in[cur_v] != 1) {
			  add_zero_vertex(cur_v);
		  }
	  }
  }
  ++dominated[v];
}

int is_in_score0_stack(int v) {
	for (int i = 0; i < zero_stack_fill_pointer; i++) {
		if (v == zero_stack[i])
			return 1;
	}
	return 0;
}

void check_score0_stack() {
	for (int i = 1; i <= v_num; i++) {
		if (v_in_c[i] == 1 && score[i] == 0 && index_in_zero_stack[i] < 0 && v_fixed[i] == 0) {
			cout << " (v_in_c[i] == 1 && score[i] == 0 && index_in_zero_stack[i] < 0) " << i << " " << file << endl;
		}
		if (index_in_zero_stack[i] >= 0 && score[i] != 0) {
			cout << "index_in_zero_stack[i] > 0 && score[i] != 0 " << i << endl;
		}
		if (index_in_zero_stack[i] >= 0 && v_in_c[i] != 1) {
			cout << "index_in_zero_stack[i] >= 0 && v_in_c[i] != 1 " << i << " " << file << endl;
		}
		if (index_in_zero_stack[i] < 0 && is_in_score0_stack(i) == 1) {
			cout << " index_in_zero_stack[i] < 0 && is_in_score0_stack(i) == 1 " << i << endl;
		}
		if (v_in_c[i] == 1 && score[i] == 0 && is_in_score0_stack(i) == 0) {
			cout << " v_in_c[i] == 1 && score[i] == 0 && is_in_score0_stack(i) == 0 " << i << endl;
		}
		if (index_in_zero_stack[i] >= 0 && zero_stack[index_in_zero_stack[i]] != i)
			cout<<"index_in_zero_stack[i] >= 0 && zero_stack[index_in_zero_stack[i]] != i" << i << " " << file <<endl;
	}
	for (int i = 0; i < zero_stack_fill_pointer; i++) {
		if (score[zero_stack[i]] != 0 || v_in_c[zero_stack[i]] == 0) {
			cout << " score[zero_stack[i]] != 0 || v_in_c[zero_stack[i]] == 0 " << i << " " << score[zero_stack[i]] << " " << v_in_c[zero_stack[i]] << endl;
		}
		if (index_in_zero_stack[zero_stack[i]] != i) {
			cout << "index_in_zero_stack[zero_stack[i]] != i" << i << endl;
		}
	}

}

bool test_score();
void Add(int v) {
  int new_score = -score[v];
  increase_dominate(v, v);//v是否由被支配1次变为被支配两次或由被支配0次变为被支配1次
  for (int i = 0; i < v_degree[v]; ++i) {
    int u = v_adj[v][i];
    increase_dominate(u, v);
  }
  score[v] = new_score;
  if (score[v] == 0 && index_in_zero_stack[v] == -1 && record_in[v] != 1) {
	  add_zero_vertex(v);
  }
  if (v_in_c[v] == 1) {
	  cout<<" v_in_c[v] == 1 "<<endl;
  }
  v_in_c[v] = 1;
  c_size++;
  candidate[candidate_size] = v;
  index_in_candidate[v] = candidate_size++;
  // sub_score[v] = new_subscore;
  if (zero_stack_fill_pointer < 0)
	  cout << "zero_stack_fill_pointer < 0 " << file << " " << seed << endl;
  //check_candidate1();
  //check_candidate();
}

void decrease_dominate(int v, int remove_vertex) {
  if (dominated[v] == 1) {//v被支配了一次然后马上变为未被支配，v一定不在解中其邻居也一定不在解中
    // all score of neighbours are make
    ++score[v];
    for (int i = 0; i < v_degree[v]; ++i) {
      int u = v_adj[v][i];
      ++score[u];
    }
    Undom(v);
  } 
  else if (dominated[v] == 2) {
    // only need to update one break
	//cover_set***
	  if (cover_set_flag)
	  {
		  if (cover_set[v] == remove_vertex) {
			  cover_set[v] = cover_set2[v];
		  }
		  score[cover_set[v]]--;
		  if (index_in_zero_stack[cover_set[v]] != -1) {
			  if (score[cover_set[v]] == 0)
				  cout << "score[v] == 0" << endl;
			  delete_zero_vertex(cover_set[v]);//从zero_stack中删除
		  }
	  }
	  else {
		  int cur_v;
		  if (v_in_c[v]) {
			  --score[v];
			  cur_v = v;
		  }
		  else {
			  // TODO: speed up this procedure
			  for (int i = 0; i < v_degree[v]; ++i) {
				  int u = v_adj[v][i];
				  if (v_in_c[u]) {
					  --score[u];
					  cur_v = u;
					  break;
				  }
			  }
		  }
		  if (index_in_zero_stack[cur_v] != -1) {
			  if (score[cur_v] == 0)
				  cout << "score[v] == 0" << endl;
			  delete_zero_vertex(cur_v);//从zero_stack中删除
		  }
	  }
  }
  else if (dominated[v] == 3 && cover_set_flag) {
	  int flag = 0;
	  if (v_in_c[v] == 1) {
		  cover_set[v] = v;
		  flag = 1;
	  }
	  for (int i = 0; i < v_degree[v]; ++i) {
		  int u = v_adj[v][i];
		  if (v_in_c[u] == 1) {
			  if (flag == 0)
			  {
				  cover_set[v] = u;
				  flag = 1;
			  }
			  else
			  {
				  cover_set2[v] = u;
				  break;
			  }
		  }
	  }
  }
  --dominated[v];
}

void Remove(int v) {
  if (index_in_zero_stack[v] != -1) {
	delete_zero_vertex(v);
  }
  int last_candidate_v = candidate[--candidate_size];
  int index = index_in_candidate[v];
  candidate[index] = last_candidate_v;
  index_in_candidate[last_candidate_v] = index;
  index_in_candidate[v] = -1;
  if (v_in_c[v] != 1) {
	  cout << " v_in_c[v] != 1 " << file << " " << seed << " " << v << endl;
  }
  v_in_c[v] = 0;
  c_size--;

  int new_score = -score[v];
  // int new_subscore = -sub_score[v];
  decrease_dominate(v, v);
  int neighbours = v_degree[v];
  for (int i = 0; i < neighbours; ++i) {
    int u = v_adj[v][i];
    decrease_dominate(u, v);
  }
  score[v] = new_score;
  if (zero_stack_fill_pointer < 0)
	  cout << "zero_stack_fill_pointer < 0 " << file << " " << seed << endl;
  //check_candidate1();
  //check_candidate();
}

void ResetCandidate() {
  int v;
  int j = 0;

  for (v = 1; v < v_num + 1; v++) {
    if (v_in_c[v] == 1 && v_fixed[v] == 0) {
      candidate[j] = v;
      index_in_candidate[v] = j;
      j++;
    } else {
      index_in_candidate[v] = -1;
    }
  }
  candidate_size = j;
}

void UpdateBestSolution() {
  int v;
  unimprove = 0;
  if (c_size < best_c_size) {
    for (int i1 = 0; i1 < candidate_size; i1++) {
	   best_candidate[i1] = candidate[i1];
	}
	best_candidate_record = candidate_size;
    best_c_size = c_size;
    best_comp_time = TimeElapsed();
    best_step = step;
	if (c_size != candidate_size + fix_must_in.size()) {
		cout << file << " " << seed << " " << c_size << " aaa " << candidate_size + fix_must_in.size() << endl;
	}
  }
  else {
	  cout << file << " " << seed << " bestsol wrong " << endl;
  }

}

bool CompareBySecond(const pair<int, int> &lhs, const pair<int, int> &rhs) {
  return lhs.second > rhs.second;
}

void AddInit(int v) {
  int pos;

  v_in_c[v] = 1;

  int new_score = -score[v];
  init_increase_dominate(v, v);
  for (int i = 0; i < v_degree[v]; ++i) {
    int u = v_adj[v][i];
    init_increase_dominate(u, v);
  }
  score[v] = new_score;
  pos = pos_in_my_heap[v];
  my_heap_remove(pos);
  if (zero_stack_fill_pointer < 0)
	  cout << "zero_stack_fill_pointer < 0 " << file << " " << seed << endl;
}

bool CompareBySecond1(const pair<int, int>& lhs, const pair<int, int>& rhs) {
	return lhs.second < rhs.second;
}

int GetMaxScore(int v)
{
	int pos = v, u;
	for (int i = 0; i < v_degree[v]; i++) {
		u = v_adj[v][i];
		if (score[u] > score[pos])
		{
			pos = u;
		}
	}
	return pos;
}

void init_increase_dominate1(int v, int added_vertex) {
	if (dominated[v] == 0) {
		--score[v];
		// all make but one break
		for (int i = 0; i < v_degree[v]; ++i) {
			int u = v_adj[v][i];
			--score[u];
		}
		//cover_set***
		if (cover_set_flag)
		{
			cover_set[v] = added_vertex;
		}
		Dom(v);
	}
	else if (dominated[v] == 1) {
		// only two vertex need to update break;
		//cover_set***
		if (cover_set_flag)
		{
			int cur_v = cover_set[v];
			++score[cur_v];
			cover_set2[v] = added_vertex;
		}
		else {
			if (v_in_c[v] == 1) {
				++score[v];
			}
			// TODO: speed up
			for (int i = 0; i < v_degree[v]; ++i) {
				int u = v_adj[v][i];
				if (v_in_c[u]) {
					++score[u];
					//break;//*****
				}
			}
		}
	}
	++dominated[v];
}

void AddInit1(int v) {
	int pos;

	v_in_c[v] = 1;

	int new_score = -score[v];
	init_increase_dominate1(v, v);
	for (int i = 0; i < v_degree[v]; ++i) {
		int u = v_adj[v][i];
		init_increase_dominate1(u, v);
	}
	score[v] = new_score;
}

void ConstructByOurAndScore(int parameter = 5) {
	fill_n(v_in_c, v_num + 1, 0);
	fill_n(score, v_num + 1, 0);
	fill_n(dominated, v_num + 1, 0);
	fill_n(v_fixed, v_num + 1, 0);//cai
	int v;
	undom_stack_fill_pointer = 0;
	// step = 0;
	int ds_size = 0;

	for (v = 1; v < v_num + 1; ++v) {
		score[v] = v_degree[v] + 1;
		time_stamp[v] = 0;
		Undom(v);
	}

	must_in_beg = 0;
	while (must_in_beg < must_in_end)
	{
		int u = must_in[must_in_beg];
		must_in_beg++;
		if (v_in_c[u] == 1) {
			cout << "v_in_c == 1" << endl;
		}
		ds_size++;

		AddInit1(u);
		v_in_c[u] = 1;
		time_stamp[u] = (llong)(~0ULL >> 1);
		v_fixed[u] = 1;

	}

	int res_num = undom_stack_fill_pointer;

	// TODO: speed up by using bin sort
	vector<pair<int, int>> v_d_vec(res_num + 1);
	res_num = 0;
	for (v = 1; v < v_num + 1; v++) {
		if (dominated[v] == 0)  v_d_vec[++res_num] = pair<int, int>(v, store_ver[v]);
	}
	sort(v_d_vec.begin() + 1, v_d_vec.end(), CompareBySecond1);


	int best_score = (int)(~0U >> 1);
	int best_neighbor;
	for (v = 1; v <= res_num; v++) {
		int v0 = v_d_vec[v].first;
		if (dominated[v0] == 0) {
			best_neighbor = v0;
			best_score = score[v0];
			for (int i = 0; i < v_degree[v0]; i++) {
				if (score[v_adj[v0][i]] > best_score)
				{
					best_score = score[v_adj[v0][i]];
					best_neighbor = v_adj[v0][i];
				}
			}
			AddInit1(best_neighbor);
			ds_size++;
		}
		if (store_ver[v0] == 6)
			break;
	}

	for (v = 1; v < v_num + 1; ++v) {
		if (!v_in_c[v] && score[v] != 0)  my_heap_insert(v);//把未覆盖的顶点加入大根堆
	}

	while (undom_stack_fill_pointer > 0) {
		//printf("undom_stack_fill_pointer = %d\n", undom_stack_fill_pointer);
		int best_v = my_heap[0];
		if (score[best_v] > 0) {
			AddInit(best_v);
			ds_size++;
		}
		//else  printf("best_v=%d  score=%d\n", best_v, score[best_v]);
	}

	// score
	c_size = ds_size;
	ResetCandidate();
	init_removeRedundant();

	best_candidate_record = candidate_size;
	best_candidate = new int[candidate_size];
	best_c_size = v_num;
}

void ConstructDS() {
	ConstructByOurAndScore(6);
	check_candidate();
	delete[] must_in;
	delete[] pos_in_must_in;
	delete[] v_fixed;
}

int CheckSolution() {
  int v;
  int *check = new int[v_num + 1];
  fill_n(check, v_num + 1, 0);
  for (int j = 0; j < fix_must_in.size(); j++) {
	  int ver = fix_must_in[j];
	  check[ver] = 1;
	  for (int k = 0; k < v_degree[ver]; k++) {
		  int ver_adj = v_adj[ver][k];
		  check[ver_adj] = 1;
	  }
  }
  for (int j = 0; j < best_candidate_record; j++) {
	  int ver = best_candidate[j];
	  check[ver] = 1;
	  for (int k = 0; k < v_degree[ver]; k++) {
		  int ver_adj = v_adj[ver][k];
		  check[ver_adj] = 1;
	  }
	  //if (v_in_c[ver] != 1)
		 // cout<< "(v_in_c[ver] != 1)" <<endl;
  }
 // int total_num1 = best_candidate_record + fix_must_in.size();
 // int total_Num2 = 0;
 // for (int i = 1; i <= v_num; i++) {
	//  if (v_in_c[i] == 1) {
	//	  total_Num2++;
	//}
 // }
 // if (total_Num2 != total_num1) {
	//  cout << "(total_Num1 != total_num1)" << endl;
 // }
  for (int i = 1; i <= v_num; i++) {
	  if (check[i] == 0) {
		  cout << file << " " << seed << " " << i << "(check[i] == 0)" << endl;
	  }
  }
  return 1;
}

void UpdateTargetSize() {
  int i, v;
  int best_score;
  int best_remove_v;

  best_remove_v = ChooseRemoveV(100);


  Remove(best_remove_v);
  time_stamp[best_remove_v] = step;
}

class Sub_score {
public:
  int operator[](int v) {
    int ss = 0;
    if (!v_in_c[v]) {
      if (dominated[v] == 1) {
        ss++;
      }
      for (int i = 0; i < v_degree[v]; ++i) {
        int u = v_adj[v][i];
        int d = dominated[u];
        if (d == 1) {
          ss++;
        }
      }
    } else {
      if (dominated[v] == 2) {
        ss--;
      }
      for (int i = 0; i < v_degree[v]; ++i) {
        int u = v_adj[v][i];
        int d = dominated[u];
        if (d == 2) {
          ss--;
        }
      }
    }
    return ss;
  }
};
Sub_score sub_score;

int ChooseRemoveV(int count) {
	int i, v;

	if (v_num <= 10000) {
		count = 15 + rand() % 10;
	}

	int best_v = candidate[rand() % candidate_size];
	for (i = 0; i < count; i++) {
		v = candidate[rand() % candidate_size];
		if (score[v] > score[best_v]) {
			best_v = v;
		}
		else if (score[v] == score[best_v] && time_stamp[v] < time_stamp[best_v]) {
			best_v = v;
		}
	}
	if (record_in[best_v] == 1)
		cout << file << " " << best_v << endl;
	return best_v;
}


int ChooseAddV() {
  // TODO: proof there is at least on avaible add_v
  int base_v, add_v;
  int best_score = -1;
  int best_add_v = -1;

  int aaa = rand() % 1000;
  if (v_num <= 10000 && aaa < 50) {
	  return undom_stack[rand() % undom_stack_fill_pointer];
  }


  for (int i = 0; i < undom_stack_fill_pointer; ++i) {
    base_v = undom_stack[i];
    for (int j = 0; j < v_degree[base_v]; ++j) {
      add_v = v_adj[base_v][j];
	  if (record_in[add_v] == 0) {
		  continue;
	  }
      /*if (conf_change[add_v] <= v_threshold[add_v]) {
        continue;
      }*/
      if (score[add_v] > best_score) {
        best_add_v = add_v;
        best_score = score[add_v];
      } else if (score[add_v] == best_score) {
        if (time_stamp[add_v] < time_stamp[best_add_v]) {
          best_add_v = add_v;
          best_score = score[add_v];
        }
      }
    }
    //if (conf_change[base_v] > v_threshold[base_v]) {
      if (score[base_v] > best_score && record_in[base_v] != 0) {
        best_add_v = base_v;
        best_score = score[base_v];
      } else if (score[base_v] == best_score) {
        if (time_stamp[base_v] < time_stamp[best_add_v]) {
          best_add_v = base_v;
          best_score = score[base_v];
        }
      }
    //}
  }
  if (best_add_v == -1) {
	  cout << file << " " << best_add_v << " " << seed << endl;
  }
  return best_add_v;
}

void delete_propagate_valid(int ver) {
	int lastver = propagate_valid[--propagate_valid_fill_pointer];//这里应该是--放前面
	int index = index_in_propagate_valid[ver];
	if (index == -1)
		cout << " wrong 1" << endl;
	index_in_propagate_valid[ver] = -1;
	propagate_valid[index] = lastver;
	if (propagate_valid_fill_pointer == index)
		index = -1;
	index_in_propagate_valid[lastver] = index;
}

void propagate_del_degree(int ver, int thre1) {//删除propagate_degree小于等于propagate_degree_thre的点
	for (int j = 0; j < v_degree[ver]; j++) {
		int ver_adj = v_adj[ver][j];
		if (index_in_propagate_valid[ver_adj] != -1) {//若删除的点的邻居还在propagate_valid中
			propagate_degree[ver_adj]--;
			unlock_propagate_edge_num--;
			if (propagate_degree[ver_adj] <= thre1 && temp_flag[ver_adj] == 0) {//若其邻居传播度小于等于thre1且没被删除则删除
				temp_store_reduction[temp_store_reduction_fill_pointer] = ver_adj;
				temp_store_reduction_fill_pointer++;
				
				temp_flag[ver_adj] = 1;
				temp_flag_stack[temp_flag_stack_fill_pointer] = ver_adj;
				temp_flag_stack_fill_pointer++;
			}
		}
	}
	delete_propagate_valid(ver);
	unlock_propagate_ver_num--;
}

int reduction_degree_thre(int propagate_degree_thre) {
	temp_store_reduction_fill_pointer = 0;
	temp_flag_stack_fill_pointer = 0;
	for (int i = 0; i < propagate_valid_fill_pointer; i++) {
		int ver = propagate_valid[i];
		if (propagate_degree[ver] <= propagate_degree_thre) {//有可能已经在这里的点又删除一遍
			temp_store_reduction[temp_store_reduction_fill_pointer] = ver;
			temp_store_reduction_fill_pointer++;

			temp_flag[ver] = 1;
			temp_flag_stack[temp_flag_stack_fill_pointer] = ver;
			temp_flag_stack_fill_pointer++;
		}
	}
	int index = 0;
	while (index != temp_store_reduction_fill_pointer) {
		propagate_del_degree(temp_store_reduction[index], propagate_degree_thre);
		store_ver[temp_store_reduction[index]] = propagate_degree_thre;
		index++;
		store_ver[temp_store_reduction[index]] = propagate_degree_thre;
		if (store_ver[temp_store_reduction[index]] < 5 && record_in[temp_store_reduction[index]] == -1) {
			total_5++;
		}
		if (store_ver[temp_store_reduction[index]] > 32 && record_in[temp_store_reduction[index]] == -1) {
			total_32++;
		}
		if (record_in[temp_store_reduction[index]] == -1) {
			total_count++;
		}
	}
	for (int i = 0; i < temp_flag_stack_fill_pointer; i++) {
		temp_flag[temp_flag_stack[i]] = 0;
	}

	temp_flag_stack_fill_pointer = 0;
	store_count++;
	return temp_store_reduction_fill_pointer;
}

void reduction() {
	int count = 1;
	while (1) {
		int ccc = reduction_degree_thre(count);
		if (ccc == 0) {
			if (unlock_propagate_edge_num != 0 && unlock_propagate_ver_num == 0) {
				cout << " reduction wrong " << file << endl;
				break;
			}
		}
		if (unlock_propagate_ver_num == 0)
			break;
		count++;
	}
	reduction_num = count;
	delete[] temp_store_reduction;
	delete[] temp_flag;
	delete[] temp_flag_stack;
	delete[] propagate_degree;
	delete[] propagate_valid;
	delete[] index_in_propagate_valid;

	//cover_set
	sum_store_ver = 0;
	int temp = 0;
	for (int v = 1; v < v_num + 1; v++)
	{
		if (record_in[v] == -1)
		{
			sum_store_ver += store_ver[v];
			temp++;
		}
	}
	sum_store_ver /= temp;
	if (sum_store_ver < 40 && v_num > 10000)
	{
		cover_set_flag = true;
		cover_set = new int[v_num + 1];
		cover_set2 = new int[v_num + 1];
		fill_n(cover_set, v_num + 1, 0);
		fill_n(cover_set2, v_num + 1, 0);
	}
	else  cover_set_flag = false;
}

int ChooseRemoveV1() {
	int i, v;
	int best_v = candidate[rand() % candidate_size];
	for (i = 0; i < 10; i++) {
		v = candidate[rand() % candidate_size];
		if (score[v] > score[best_v]) {
			best_v = v;
		}
		else if (score[v] == score[best_v] && time_stamp[v] < time_stamp[best_v]) {
			best_v = v;
		}
	}
	if (record_in[best_v] == 1)
		cout << file << " " << best_v << endl;
	return best_v;
}

int ChooseRemoveV2() {
	int i, v;
	int best_v = candidate[rand() % candidate_size];
	int best_score = 0 - (score[best_v] - 1);
	int best_age = step - time_stamp[best_v];


	if (reduction_num  < 40 || v_num < 10000) {
		for (i = 0; i < 10; i++) {
			v = candidate[rand() % candidate_size];

			int score_v = 0 - (score[v] - 1);
			int age_v = step - time_stamp[v];

			if (age_v * best_score > best_age * score_v) {
				best_v = v;
				best_age = age_v;
				best_score = score_v;
			}


		}
		if (record_in[best_v] == 1)
			cout << file << " " << best_v << endl;
	}
	else {
		for (i = 0; i < 10; i++) {
			v = candidate[rand() % candidate_size];
			if (time_stamp[v] < time_stamp[best_v]) {
				best_v = v;
			}
		}
	}
	
	return best_v;
}

void unlock() {
	fix_flag = 0;
	for (int i = 0; i < equal_Reduction.size(); i++) {
		record_in[equal_Reduction[i]] = -1;
	}
}

void LocalSearch() {
  int remove_v, add_v, remove_v2;
  step = 1LL;
  try_step = 100;
  int flag = 0;
  int flag1 = 0;
  int count_update = 0;
  while (true) {
	  //
	t1 = 45;
	t2 = 0;
    if (undom_stack_fill_pointer == 0) {
      UpdateBestSolution();
	  for (int i = 1; i <= v_num; i++) {
		  if (index_in_undom_stack[i] > 0)
			  cout << " (index_in_undom_stack[i] > 0) " << i << endl;
	  }

      UpdateTargetSize();//移除一个顶点
	  flag = 0;
      continue;
    }
    if (step % try_step == 0) {
	  double time__ = TimeElapsed();
      if (time__ >= cutoff_time) {
        return;
      }
	  if (time__ > 1000 && flag_1000_out == 1) {
		  cout << 1000 << file << " " << seed << " " << best_c_size << " " << best_comp_time << " " << sum_store_ver << " " <<  cover_set_flag << " " << step << endl;
		  flag_1000_out = 0;
	  }
	  if (unimprove > best_c_size && fix_flag == 1) {
		  unlock();
	  }
    }
	remove_v = ChooseRemoveV(45);
	Remove(remove_v);
	time_stamp[remove_v] = step;

	if (candidate_size > 0)
	{
		if (unimprove < 2500) {
			remove_v = ChooseRemoveV1();
		}
		else if (unimprove >= 2500 && unimprove < 5000) {
			remove_v = candidate[rand() % candidate_size];
		}
		else {
			remove_v = ChooseRemoveV2();
		}

		Remove(remove_v);
		time_stamp[remove_v] = step;

		if (undom_stack_fill_pointer != 0) {
			add_v = ChooseAddV();
			Add(add_v);
			time_stamp[add_v] = step;
		}
	}

	if (undom_stack_fill_pointer != 0) {
		add_v = ChooseAddV();
		Add(add_v);
		time_stamp[add_v] = step;
	}

	unimprove++;
    step++;
    if (undom_stack_fill_pointer == 0) {
      removeRedundant();
    }
  }
}
