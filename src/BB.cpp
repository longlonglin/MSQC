#include <sys/time.h>
#include <sys/types.h>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

using namespace std;

const char* preOutPutFileName;

void countSort(vector<int>& arr)
{
	int len = arr.size();
	if (len == 0)
		return;
	vector<int> tempArr(arr.begin(), arr.end());
	int min = tempArr[0], max = min;
	for (int i = 1; i < len; ++i)
	{
		if (min > tempArr[i])
			min = tempArr[i];
		if (max < tempArr[i])
			max = tempArr[i];
	}
	const int k = max - min + 1;
	int count[k];
	for (int i = 0; i < k; i++)
	{
		count[i] = 0;
	}
	for (int i = 0; i < len; ++i)
		++count[tempArr[i] - min];
	for (int i = 1; i < k; ++i)
		count[i] += count[i - 1];
	for (int i = len - 1; i >= 0; --i)
		arr[--count[tempArr[i] - min]] = tempArr[i];

}
void delRepeat(vector<int>& arr)
{
	int len = arr.size();
	if (len == 0)
		return;
	int index = 0;
	for (int i = 0; i < len; i++)
	{
		if (arr[index] != arr[i])
		{
			arr[++index] = arr[i];
		}
	}
	arr.resize(index + 1);
}
void countSort(vector<pair<int, int> >& arr)
{
	int len = arr.size();
	if (len == 0)
		return;
	vector<pair<int, int> > tempArr(arr.begin(), arr.end());
	int min = tempArr[0].first, max = min;
	for (int i = 1; i < len; ++i)
	{
		if (min > tempArr[i].first)
			min = tempArr[i].first;
		if (max < tempArr[i].first)
			max = tempArr[i].first;
	}
	const int k = max - min + 1;
	int count[k];
	for (int i = 0; i < k; i++)
	{
		count[i] = 0;
	}
	for (int i = 0; i < len; ++i)
		++count[tempArr[i].first - min];
	for (int i = 1; i < k; ++i)
		count[i] += count[i - 1];
	for (int i = len - 1; i >= 0; --i)
	{
		arr[--count[tempArr[i].first - min]] = tempArr[i];
	}
}
void countSort(vector<vector<int> >& arr)
{
	int len = arr.size();
	if (len == 0)
		return;
	for (int i = 0; i < len; i++)
	{
		countSort(arr[i]);
	}
	vector<vector<int> > tempArr(arr.begin(), arr.end());
	int min = tempArr[0][0], max = min;
	for (int i = 1; i < len; ++i)
	{
		if (min > tempArr[i][0])
			min = tempArr[i][0];
		if (max < tempArr[i][0])
			max = tempArr[i][0];
	}
	const int k = max - min + 1;
	int count[k];
	for (int i = 0; i < k; i++)
	{
		count[i] = 0;
	}
	for (int i = 0; i < len; ++i)
		++count[tempArr[i][0] - min];
	for (int i = 1; i < k; ++i)
		count[i] += count[i - 1];
	for (int i = len - 1; i >= 0; --i)
	{
		int tmp = --count[tempArr[i][0] - min];
		arr[tmp].clear();
		for (int j = 0; j < tempArr[i].size(); j++)
			arr[tmp].push_back(tempArr[i][j]);
	}
}
bool cmp1(pair<int, int> a, pair<int, int> b)
{
	if (a.first != b.first)
		return a.first < b.first;
	else
		return a.second < b.second;
}
vector<int> Union(vector<int> a, vector<int> b)
{
	countSort(a);
	countSort(b);
	if (a.size() == 0)
		return b;
	else if (b.size() == 0)
		return a;
	vector<int> result;
	set_union(a.begin(), a.end(), b.begin(), b.end(), inserter(result, result.begin()));
	return result;
}
vector<int> Difference(vector<int> a, vector<int> b)
{
	countSort(a);
	countSort(b);
	if (b.size() == 0)
		return a;
	vector<int> result;
	set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(result, result.begin()));
	return result;
}
vector<int> Intersection(vector<int> a, vector<int> b)
{
	countSort(a);
	countSort(b);
	vector<int> result;
	set_intersection(a.begin(), a.end(), b.begin(), b.end(), inserter(result, result.begin()));
	return result;
}
vector<vector<int> > Union(vector<vector<int> > a, vector<vector<int> > b)
{
	countSort(a);
	countSort(b);
	if (a.size() == 0)
		return b;
	else if (b.size() == 0)
		return a;
	vector<vector<int> > result;
	int i, j;
	for (i = 0, j = 0; i < a.size() && j < b.size();)
	{
		if (a[i][0] < b[j][0])
		{
			result.push_back(a[i]);
			i++;
		}
		else if (a[i][0] > b[j][0])
		{
			result.push_back(b[j]);
			j++;
		}
		else
		{
			bool flag = true;
			if (a[i].size() != b[j].size())
				flag = false;
			else
			{
				for (int k = 0; k < a[i].size(); k++)
				{
					if (a[i][k] != b[j][k])
					{
						flag = false;
						break;
					}
				}
			}
			if (flag == true)
			{
				result.push_back(a[i]);
			}
			else
			{
				result.push_back(a[i]);
				result.push_back(b[j]);
			}
			i++, j++;
		}
	}
	if (i == a.size())
	{
		for (; j < b.size(); j++)
			result.push_back(b[j]);
	}
	if (j == b.size())
	{
		for (; i < a.size(); i++)
			result.push_back(a[i]);
	}
	return result;
}
vector<pair<int, int> > Union(vector<pair<int, int> > a, vector<pair<int, int> > b)
{
	sort(a.begin(), a.end(), cmp1);
	sort(b.begin(), b.end(), cmp1);
	if (a.size() == 0)
		return b;
	else if (b.size() == 0)
		return a;
	vector<pair<int, int> > result;
	int i, j;
	for (i = 0, j = 0; i < a.size() && j < b.size();)
	{
		if (a[i].first < b[j].first)
		{
			result.push_back(a[i]);
			i++;
		}
		else if (a[i].first > b[j].first)
		{
			result.push_back(b[j]);
			j++;
		}
		else if (a[i].first == b[j].first && a[i].second < b[j].second)
		{
			result.push_back(a[i]);
			i++;
		}
		else if (a[i].first == b[j].first && a[i].second > b[j].second)
		{
			result.push_back(b[j]);
			j++;
		}
		else
		{
			result.push_back(a[i]);
			i++, j++;
		}
	}
	if (i == a.size())
	{
		for (; j < b.size(); j++)
			result.push_back(b[j]);
	}
	if (j == b.size())
	{
		for (; i < a.size(); i++)
			result.push_back(a[i]);
	}
	return result;
}
vector<pair<int, int> > Difference(vector<pair<int, int> > a, vector<pair<int, int> > b)
{
	sort(a.begin(), a.end(), cmp1);
	sort(b.begin(), b.end(), cmp1);
	if (b.size() == 0)
		return a;
	vector<pair<int, int> > result;
	int i, j;
	for (i = 0, j = 0; i < a.size() && j < b.size();)
	{
		if (a[i].first < b[j].first)
		{
			result.push_back(a[i]);
			i++;
		}
		else if (a[i].first > b[j].first)
		{
			j++;
		}
		else if (a[i].first == b[j].first && a[i].second < b[j].second)
		{
			result.push_back(a[i]);
			i++;
		}
		else if (a[i].first == b[j].first && a[i].second > b[j].second)
		{
			j++;
		}
		else
		{
			i++, j++;
		}
	}
	if (j == b.size())
	{
		for (; i < a.size(); i++)
			result.push_back(a[i]);
	}
	return result;
}
pair<int, int> Intersection(pair<int, int> a, pair<int, int> b)
{
	pair<int, int> result;
	if (a.first > b.second || a.second < b.first)
		result = make_pair(0, 0);
	else
		result = make_pair(max(a.first, b.first), min(a.second, b.second));
	return result;
}
vector<pair<int, int> > Intersection(vector<pair<int, int> > a, vector<pair<int, int> > b)
{
	sort(a.begin(), a.end(), cmp1);
	sort(b.begin(), b.end(), cmp1);
	vector<pair<int, int> > result;
	int i, j;
	for (i = 0, j = 0; i < a.size() && j < b.size();)
	{
		if (a[i].first < b[j].first)
		{
			i++;
		}
		else if (a[i].first > b[j].first)
		{
			j++;
		}
		else if (a[i].first == b[j].first && a[i].second < b[j].second)
		{
			i++;
		}
		else if (a[i].first == b[j].first && a[i].second > b[j].second)
		{
			j++;
		}
		else
		{
			result.push_back(a[i]);
			i++, j++;
		}
	}
	return result;
}

struct BBPre {

	vector<vector<pair<pair<int, int>, pair<int, int> > > > preneighbors;
	vector<vector<pair<pair<int, int>, pair<int, int> > > > neighbors;
	vector<vector<pair<int, int> > > InputGraph;
	int delta, Count = 0;
	double gama, ro, nowlb = 0;
	int V, E, c;
	double M;
	double tauV;
	string Str;
	vector<vector<int> > ds;
	timeval start_at, end_at;

	
	int allTimestamps(vector<pair<int, int> > a)
	{
		if (a.empty())
			return 0;
		countSort(a);
		int sum = 0;
		sum += a[0].second - a[0].first + 1;
		int index = a[0].second;
		for (int i = 1; i < a.size(); i++)
		{
			if (a[i].second > index)
			{
				if (a[i].first <= index)
				{
					sum += a[i].second - index;
				}
				else
				{
					sum += a[i].second - a[i].first + 1;
				}
				index = a[i].second;
			}
		}
		return sum;
	}

	void init(istream& in, double m, int n, double p, int q, double d)
	{
		gama = m;
		delta = n;
		ro = p;
		c = q;
		M = d;
		in >> V >> E;
		InputGraph.resize(V);
		neighbors.resize(V);
		preneighbors.resize(V);
		for (int i = 0; i < E; i++)
		{
			int u, v, t, w;
			in >> u >> v >> t;
			if (u != v)
			{
				InputGraph[u].push_back(make_pair(v, t));
				InputGraph[v].push_back(make_pair(u, t));
			}
		}
		for (int i = 0; i < V; i++)
		{
			if (InputGraph[i].empty())
				continue;
			if (InputGraph[i].size() == 1)
			{
				preneighbors[i].push_back(make_pair(make_pair(InputGraph[i][0].first, 0), make_pair(InputGraph[i][0].second, InputGraph[i][0].second)));
				continue;
			}
			int a = 0, b = 0;
			while (a < InputGraph[i].size())
			{
				if (InputGraph[i].size() == 1)
				{
					preneighbors[i].push_back(make_pair(make_pair(InputGraph[i][a].first, 1), make_pair(InputGraph[i][a].second, InputGraph[i][a].second)));
					break;
				}
				while (InputGraph[i][a].first == InputGraph[i][a + 1].first)
				{
					a++;
					if (a == InputGraph[i].size() - 1)
						break;
				}
				vector<int> tmp;
				for (int j = b; j <= a; j++)
				{
					tmp.push_back(InputGraph[i][j].second);
				}
				int u = 0, v = 0, num = 0;
				while (v < tmp.size())
				{
					if (tmp.size() == 1)
					{
						preneighbors[i].push_back(make_pair(make_pair(InputGraph[i][a].first, num++), make_pair(tmp[0], tmp[0])));
						break;
					}
					while (tmp[v] == tmp[v + 1] - 1)
					{
						v++;
						if (v == tmp.size() - 1)
							break;
					}
					preneighbors[i].push_back(make_pair(make_pair(InputGraph[i][a].first, num++), make_pair(tmp[u], tmp[v])));
					v++;
					if (v == tmp.size() - 1)
					{
						preneighbors[i].push_back(make_pair(make_pair(InputGraph[i][a].first, num++), make_pair(tmp[v], tmp[v])));
						break;
					}
					u = v;
				}
				a++;
				if (a == InputGraph[i].size() - 1)
				{
					preneighbors[i].push_back(make_pair(make_pair(InputGraph[i][a].first, num++), make_pair(InputGraph[i][a].second, InputGraph[i][a].second)));
					break;
				}
				b = a;
			}
		}
	}
	int DUT(int u, const vector<int>& Gu, pair<int, int> T, int t)
	{
		int sum = 0;
		for (int i = 0; i < preneighbors[u].size(); i++)
		{
			if (t<preneighbors[u][i].second.first || t>preneighbors[u][i].second.second)
				continue;
			else
			{
				sum += 1;
			}
		}
		return sum;
	}
	bool cmp2(pair<double, int> a, pair<double, int> b)
	{
		return a.first < b.first;
	}
	vector<pair<int, int> > CDI(int u, const vector<int>& U, vector<pair<int, int> >& tau, vector<int> Sel)
	{
		vector<pair<int, int> > tmptau;
		for (int i = 0; i < tau.size(); i++)
		{
			if (tau[i].second - tau[i].first + 1 >= delta)
				tmptau.push_back(tau[i]);
		}
		tau.assign(tmptau.begin(), tmptau.end());
		vector<pair<int, int> > result;
		for (int i = 0; i < tau.size(); i++)
		{
			vector<double> dutU;
			vector<pair<double, int> > S;
			for (int j = 0; j < tau[i].second - tau[i].first + 1; j++)
			{
				dutU.push_back((double)DUT(u, U, tau[i], tau[i].first + j) - gama * (max(max(M, (double)Sel.size()), nowlb) - 1));
				if (j == 0)
				{
					S.push_back(make_pair(dutU[0], 0));
				}
				else
					S.push_back(make_pair(S[j - 1].first + dutU[j], j));
			}
			sort(S.begin(), S.end(), cmp2);
			int min_index = numeric_limits<int>::max();
			vector<pair<int, int> > Iu;
			for (int j = 1; j < tau[i].second - tau[i].first + 1; j++)
			{
				if (S[j - 1].second < min_index)
					min_index = S[j - 1].second;
				if (S[j].second - min_index + 1 >= delta)
				{
					vector<pair<int, int> > I;
					I.push_back(make_pair(tau[i].first + min_index, tau[i].first + S[j].second));
					vector<pair<int, int> > tmpIu = Union(I, Iu);
					Iu.assign(tmpIu.begin(), tmpIu.end());
				}
			}
			countSort(Iu);
			vector<pair<int, int> > Iut;
			if (Iu.size() == 0)
				continue;
			Iut.push_back(Iu[0]);
			for (int j = 1; j < Iu.size(); j++)
			{
				if (Iu[j].first == Iut.back().first && Iu[j].second > Iut.back().second)
				{
					Iut.pop_back();
					Iut.push_back(Iu[j]);
				}
				else if (Iu[j].first > Iut.back().first&& Iu[j].second > Iut.back().second)
				{
					Iut.push_back(Iu[j]);
				}
			}
			vector<pair<int, int> > tmpresult = Union(result, Iut);
			result.assign(tmpresult.begin(), tmpresult.end());
		}
		return result;
	}
	vector<pair<int, int> > timeZon(vector<int> U)
	{
		int len = U.size();
		vector<bool> flag(V, false);
		for (int i = 0; i < U.size(); i++)
		{
			flag[U[i]] = true;
		}
		vector<pair<int, int> > tmptau;
		vector<pair<int, int> > tau;
		for (int i = 0; i < len; i++)
		{
			for (int j = 0; j < preneighbors[U[i]].size(); j++)
			{
				if (flag[preneighbors[U[i]][j].first.first] == true)
				{
					tmptau.push_back(preneighbors[U[i]][j].second);
				}
			}
		}
		if (tmptau.size() == 1 || tmptau.empty())
			return tmptau;
		sort(tmptau.begin(), tmptau.end(), cmp1);
		tau.push_back(tmptau[0]);
		for (int i = 1; i < tmptau.size(); i++)
		{
			if (tmptau[i].first > tau.back().second)
			{
				tau.push_back(tmptau[i]);
			}
			else if (tmptau[i].first <= tau.back().second && tmptau[i].second > tau.back().second)
			{
				tau.back().second = tmptau[i].second;
			}
		}
		return tau;
	}
	void prePrune()
	{
		gettimeofday(&start_at, 0);
		vector<int> U;
		vector<int> Sel;
		for (int i = 0; i < V; i++)
		{
			if (!preneighbors[i].empty())
			{
				U.push_back(i);
			}
		}
		int presize = U.size();
		vector<pair<int, int> > tauS = timeZon(Sel);
		vector<int> mark(V, 2);
		vector<pair<int, int> > timeU = timeZon(U);
		tauV = allTimestamps(timeU);
		for (int i = 0; i < U.size(); i++)
		{
			mark[U[i]] = 1;
			vector<pair<int, int> > Iu = CDI(U[i], U, timeU, Sel);
			if (allTimestamps(Iu) < ro * tauV)
			{
				mark[U[i]] = 2;
			}
		}
		gettimeofday(&end_at, 0);
		int Nodenum = 0;
		E = 0;
		for (int i = 0; i < V; i++)
		{
			if (preneighbors[i].empty() || mark[i] != 1)
				continue;
			for (int j = 0; j < preneighbors[i].size(); j++)
			{
				if (mark[preneighbors[i][j].first.first] != 1)
					continue;
				neighbors[i].push_back(preneighbors[i][j]);
				E++;
				if (i > Nodenum)
					Nodenum = i;
			}
		}
		V = Nodenum + 1;
		int num = 0;
		for (int i = 0; i < neighbors.size(); i++)
		{
			if (!neighbors[i].empty())
				num++;
		}
		ofstream out("TGRA.txt", ios::app);
		out << endl;
		out << endl;
		out << Str << " game:" << gama << " rho:" << ro << " theta:" << M << endl;
		out << "reading time: "
			<< (end_at.tv_sec - start_at.tv_sec) * 1000 +
			(end_at.tv_usec - start_at.tv_usec) / 1000
			<< "ms" << endl;
		out << "Prenum:" << presize << "  Nodenum:" << num;
		out << "  rate:" << setprecision(5) << (double)num / (double)presize * 100 << "%" << endl;
		out.close();
	}
	void outPut()
	{
		string filename;
		filename = "new." + Str;
		preOutPutFileName = filename.c_str();
		ofstream out(filename.c_str());
		out << V << " " << E << endl;
		for (int i = 0; i < neighbors.size(); i++)
		{
			if (neighbors[i].empty())
				continue;
			for (int j = 0; j < neighbors[i].size(); j++)
			{
				out << i << " " << neighbors[i][j].first.first << " " << neighbors[i][j].second.first << " " << neighbors[i][j].second.second << endl;
			}
		}
	}
};

struct BBProcess {

vector<vector<pair<pair<int, int>, pair<int, int> > > > neighbors;
vector<vector<pair<int, int> > > StoreCDI;
vector<pair<int, int> > currentTimezone;
vector<bool> currentGraph;
int delta, Count = 0;
double gama, ro, nowlb;
int V, E, mode;
double M;
double tauV;
string Str;
int timetau;
vector<vector<int> > ds;
timeval start_at, end_at;

bool is_linked(int a, int b)
{
	if (neighbors[a].empty() || neighbors[b].empty())
		return false;
	for (int i = 0; i < neighbors[a].size(); i++)
	{
		if (neighbors[a][i].first.first == b)
			return true;
	}
	return false;
}
vector<pair<int, int> > timeZon(vector<int> U)
{
	int len = U.size();
	vector<bool> flag(V, false);
	for (int i = 0; i < U.size(); i++)
	{
		flag[U[i]] = true;
	}
	vector<pair<int, int> > tmptau;
	vector<pair<int, int> > tau;
	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < neighbors[U[i]].size(); j++)
		{
			if (flag[neighbors[U[i]][j].first.first] == true)
			{
				tmptau.push_back(neighbors[U[i]][j].second);
			}
		}
	}
	if (tmptau.size() == 1 || tmptau.empty())
		return tmptau;
	sort(tmptau.begin(), tmptau.end(), cmp1);
	tau.push_back(tmptau[0]);
	for (int i = 1; i < tmptau.size(); i++)
	{
		if (tmptau[i].first > tau.back().second)
		{
			tau.push_back(tmptau[i]);
		}
		else if (tmptau[i].first <= tau.back().second && tmptau[i].second > tau.back().second)
		{
			tau.back().second = tmptau[i].second;
		}
	}
	return tau;
}

void init(istream& in, double m, double p, double d)
{
	gama = m;
	delta = 1;
	ro = p;
	M = d;
	in >> V >> E;
	neighbors.resize(V);
	StoreCDI.resize(V);
	currentGraph.resize(V);
	for (int i = 0; i < E; i++)
	{
		int u, v, t, w;
		in >> u >> v >> w >> t;
		neighbors[u].push_back(make_pair(make_pair(v, 1), make_pair(w, t)));
	}
}

int DUT(int u, const vector<int>& Gu, pair<int, int> T, int t)
{
	if (Count == 1)
	{
		int sum = 0;
		for (int i = 0; i < neighbors[u].size(); i++)
		{
			if (t<neighbors[u][i].second.first || t>neighbors[u][i].second.second)
				continue;
			else
			{
				sum += 1;
			}
		}
		return sum;
	}
	vector<bool> flag(V, false);
	for (int i = 0; i < Gu.size(); i++)
	{
		flag[Gu[i]] = true;
	}
	int sum = 0;
	for (int i = 0; i < neighbors[u].size(); i++)
	{
		if (flag[neighbors[u][i].first.first] == false)
			continue;
		if (t<neighbors[u][i].second.first || t>neighbors[u][i].second.second)
			continue;
		else
		{
			sum += 1;
		}
	}
	return sum;
}
int DUT2(int u, const vector<int>& Gu, pair<int, int> T, int t)
{
	int sum = 0;
	for (int i = 0; i < neighbors[u].size(); i++)
	{
		if (currentGraph[neighbors[u][i].first.first] == false)
			continue;
		if (t<neighbors[u][i].second.first || t>neighbors[u][i].second.second)
			continue;
		else
		{
			sum += 1;
		}
	}
	return sum;
}
int DUT(int u, const vector<int>& Gu, pair<int, int> T)
{
	vector<bool> flag(V, false);
	for (int i = 0; i < Gu.size(); i++)
	{
		flag[Gu[i]] = true;
	}
	int sum = 0;
	for (int i = 0; i < neighbors[u].size(); i++)
	{
		if (flag[neighbors[u][i].first.first] == false)
			continue;
		pair<int, int> neighT = Intersection(T, neighbors[u][i].second);
		if (neighT.first != 0)
		{
			sum += neighT.second - neighT.first + 1;
		}
	}
	return sum;
}
int DUT2(int u, const vector<int>& Gu, pair<int, int> T)
{
	int sum = 0;
	for (int i = 0; i < neighbors[u].size(); i++)
	{
		if (currentGraph[neighbors[u][i].first.first] == false)
			continue;
		pair<int, int> neighT = Intersection(T, neighbors[u][i].second);
		if (neighT.first != 0)
		{
			sum += neighT.second - neighT.first + 1;
		}
	}
	return sum;
}
double MinAveSubArray(int u, const vector<int>& S, pair<int, int> T)
{
	vector<int> dut, sum, p;
	dut.push_back(0);
	sum.push_back(0);
	int n = T.second - T.first + 1;
	p.resize(n + 1);
	for (int i = T.first; i <= T.second; i++)
	{
		dut.push_back(DUT(u, S, T, i));
		if (i == T.first)
		{
			sum.push_back(dut[1]);
		}
		else
		{
			sum.push_back(sum[i - T.first] + dut[i - T.first + 1]);
		}
	}
	int ansL = 1, ansR = delta, i = 0, j = 0;
	for (int t = delta; t <= n; t++)
	{
		while (i < j - 1 && ((sum[p[j - 1]] - sum[p[j - 2] - 1]) * (t - delta - p[j - 1] + 1) - (sum[t - delta] - sum[p[j - 1] - 1]) * (p[j - 1] - p[j - 2] + 1))>0)
			j--;
		p[j++] = t - delta + 1;
		while (i < j - 1 && ((sum[t] - sum[p[i + 1] - 1]) * (t - p[i] + 1) - (sum[t] - sum[p[i] - 1]) * (t - p[i + 1] + 1)) >= 0)
			i++;
		int c = (sum[t] - sum[p[i] - 1]) * (ansR - ansL + 1) - (sum[ansR] - sum[ansL - 1]) * (t - p[i] + 1);
		if (c > 0 || (c == 0 && t - p[i] < ansR - ansL))
		{
			ansL = p[i];
			ansR = t;
		}
	}
	int result = 0;
	for (int i = ansL; i <= ansR; i++)
	{
		result += dut[i];
	}
	return (double)result / (double)(ansR - ansL + 1);
}
double AUT_min(const vector<int>& time, const vector<int>& H, const vector<int>& S)
{
	double result = numeric_limits<double>::min();
	vector<pair<int, int> > tauH = timeZon(time);
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < tauH.size(); j++)
		{
			if (tauH[j].second - tauH[j].first + 1 < delta)
				continue;
			if (MinAveSubArray(S[i], H, tauH[j]) > result)
			{
				result = MinAveSubArray(S[i], H, tauH[j]);
			}
		}
	}
	return result;
}
/**
double lb(const vector<int> &Gh,const vector<int> &S)
{
	if(S.empty())
		return 0;
	return (double)(S.size()-gama-AUT_min(Gh,S,S))/(1-gama);
}
double ub(const vector<int> &Gh,const vector<int> &S)
{
	if(S.empty())
		return (double)Gh.size();
	return AUT_min(Gh,Gh,S)/gama+1;
}**/
bool cmp3(double a, double b)
{
	return a > b;
}
double aut_max(const vector<int>& time, const vector<int>& H, const vector<int>& S, int k)
{
	vector<pair<int, int> > tauH = timeZon(time);
	vector<double> topk;
	for (int i = 0; i < S.size(); i++)
	{
		double tmpresult = 0;
		for (int j = 0; j < tauH.size(); j++)
		{
			if (tauH[j].second - tauH[j].first + 1 < delta)
				continue;
			if (MinAveSubArray(S[i], H, tauH[j]) > tmpresult)
			{
				tmpresult = MinAveSubArray(S[i], H, tauH[j]);
			}
		}
		topk.push_back(tmpresult);
	}
	sort(topk.begin(), topk.end(), cmp3);
	vector<double> qianzhihe;
	double x = AUT_min(time, H, H);
	int len = topk.size();
	for (int i = 0; i < len; i++)
	{
		if (i == 0)
			qianzhihe.push_back(topk[0]);
		else
		{
			qianzhihe.push_back(qianzhihe[i - 1] + topk[i]);
		}
	}
	double Hsize = H.size();
	if (k == 1)
	{
		for (int i = len - 1; i >= 0; i--)
		{
			if (x + qianzhihe[i] >= Hsize * gama * (Hsize + i))
				return i + 1 + Hsize;
		}
		return Hsize;
	}
	else
	{
		double firstlb;
		int i;
		for (i = 0; i <= len - 1; i++)
		{
			if (x + qianzhihe[i] >= Hsize * gama * (Hsize + i))
			{
				firstlb = i + 1 + Hsize;
				break;
			}
		}
		if (i == len)
			firstlb = Hsize;
		return firstlb;
	}
}
double ub(const vector<int>& Gh, const vector<int>& S)
{
	if (S.empty())
		return (double)Gh.size();
	return aut_max(Gh, S, Difference(Gh, S), 1);
}
double lb(const vector<int>& Gh, const vector<int>& S)
{
	if (S.empty())
		return M;
	return aut_max(Gh, S, Difference(Gh, S), 2);
}
bool cmp2(pair<double, int> a, pair<double, int> b)
{
	return a.first < b.first;
}
vector<pair<int, int> > CDI(int u, const vector<int>& U, vector<pair<int, int> >& tau, vector<int> Sel)
{
	vector<pair<int, int> > result;
	for (int i = 0; i < tau.size(); i++)
	{
		vector<double> dutU;
		vector<pair<double, int> > S;
		for (int j = 0; j < tau[i].second - tau[i].first + 1; j++)
		{
			dutU.push_back((double)DUT(u, U, tau[i], tau[i].first + j) - gama * (max(max(M, (double)Sel.size()), nowlb) - 1));
			if (j == 0)
			{
				S.push_back(make_pair(dutU[0], 0));
			}
			else
				S.push_back(make_pair(S[j - 1].first + dutU[j], j));
		}
		stable_sort(S.begin(), S.end(), cmp2);
		int min_index = numeric_limits<int>::max();
		int index = numeric_limits<int>::max();
		vector<pair<int, int> > Iu;
		for (int j = 1; j < tau[i].second - tau[i].first + 1; j++)
		{
			if (S[j - 1].second < min_index)
			{
				min_index = S[j - 1].second;
				index = j - 1;
			}
			if (S[j].second - min_index + 1 >= delta)
			{
				if (S[j].first - S[index].first + dutU[min_index] >= 0)
				{
					vector<pair<int, int> > I;
					I.push_back(make_pair(tau[i].first + min_index, tau[i].first + S[j].second));
					vector<pair<int, int> > tmpIu = Union(I, Iu);
					Iu.assign(tmpIu.begin(), tmpIu.end());
				}
				else
				{
					vector<pair<int, int> > I;
					if (S[j].second - min_index == delta)
						continue;
					I.push_back(make_pair(tau[i].first + min_index + 1, tau[i].first + S[j].second));
					vector<pair<int, int> > tmpIu = Union(I, Iu);
					Iu.assign(tmpIu.begin(), tmpIu.end());
				}
			}
		}
		countSort(Iu);
		vector<pair<int, int> > Iut;
		if (Iu.size() == 0)
			continue;
		Iut.push_back(Iu[0]);
		for (int j = 1; j < Iu.size(); j++)
		{
			if (Iu[j].first == Iut.back().first && Iu[j].second > Iut.back().second)
			{
				Iut.pop_back();
				Iut.push_back(Iu[j]);
			}
			else if (Iu[j].first > Iut.back().first&& Iu[j].second > Iut.back().second)
			{
				Iut.push_back(Iu[j]);
			}
		}
		vector<pair<int, int> > tmpresult = Union(result, Iut);
		result.assign(tmpresult.begin(), tmpresult.end());
	}
	return result;
}
vector<pair<int, int> > CDI2(int u, const vector<int>& U, vector<pair<int, int> >& tau, vector<int> Sel)
{
	vector<pair<int, int> > result;
	for (int i = 0; i < tau.size(); i++)
	{
		vector<double> dutU;
		vector<pair<double, int> > S;
		for (int j = 0; j < tau[i].second - tau[i].first + 1; j++)
		{
			dutU.push_back((double)DUT2(u, U, tau[i], tau[i].first + j) - gama * (max(max(M, (double)Sel.size()), nowlb) - 1));
			if (j == 0)
			{
				S.push_back(make_pair(dutU[0], 0));
			}
			else
				S.push_back(make_pair(S[j - 1].first + dutU[j], j));
		}
		stable_sort(S.begin(), S.end(), cmp2);
		int min_index = numeric_limits<int>::max();
		int index = numeric_limits<int>::max();
		vector<pair<int, int> > Iu;
		for (int j = 1; j < tau[i].second - tau[i].first + 1; j++)
		{
			if (S[j - 1].second < min_index)
			{
				min_index = S[j - 1].second;
				index = j - 1;
			}
			if (S[j].second - min_index + 1 >= delta)
			{
				if (S[j].first - S[index].first + dutU[min_index] >= 0)
				{
					vector<pair<int, int> > I;
					I.push_back(make_pair(tau[i].first + min_index, tau[i].first + S[j].second));
					vector<pair<int, int> > tmpIu = Union(I, Iu);
					Iu.assign(tmpIu.begin(), tmpIu.end());
				}
				else
				{
					vector<pair<int, int> > I;
					if (S[j].second - min_index == delta)
						continue;
					I.push_back(make_pair(tau[i].first + min_index + 1, tau[i].first + S[j].second));
					vector<pair<int, int> > tmpIu = Union(I, Iu);
					Iu.assign(tmpIu.begin(), tmpIu.end());
				}
			}
		}
		countSort(Iu);
		vector<pair<int, int> > Iut;
		if (Iu.size() == 0)
			continue;
		Iut.push_back(Iu[0]);
		for (int j = 1; j < Iu.size(); j++)
		{
			if (Iu[j].first == Iut.back().first && Iu[j].second > Iut.back().second)
			{
				Iut.pop_back();
				Iut.push_back(Iu[j]);
			}
			else if (Iu[j].first > Iut.back().first&& Iu[j].second > Iut.back().second)
			{
				Iut.push_back(Iu[j]);
			}
		}
		vector<pair<int, int> > tmpresult = Union(result, Iut);
		result.assign(tmpresult.begin(), tmpresult.end());
	}
	return result;
}
int allTimestamps(vector<pair<int, int> > a)
{
	if (a.empty())
		return 0;
	countSort(a);
	int sum = 0;
	sum += a[0].second - a[0].first + 1;
	int index = a[0].second;
	for (int i = 1; i < a.size(); i++)
	{
		if (a[i].second > index)
		{
			if (a[i].first <= index)
			{
				sum += a[i].second - index;
			}
			else
			{
				sum += a[i].second - a[i].first + 1;
			}
			index = a[i].second;
		}
	}
	return sum;
}
int quality(int v, const vector<int>& C, const vector<int>& S)
{
	int thegema = 0;
	for (int i = 0; i < currentTimezone.size(); i++)
	{
		if (mode == 4) {
			thegema += (DUT2(v, C, currentTimezone[i]) * allTimestamps(StoreCDI[v]));
		}
		if (mode == 5) {
			thegema += (DUT2(v, C, currentTimezone[i]) - allTimestamps(StoreCDI[v]));
		}
	}
	return thegema;
}
int qualitySort(vector<int> U, const vector<int>& C, const vector<int>& S)
{
	if (mode < 4)
	{
		//random
		int ranNum = rand() % U.size();
		return U[ranNum];
	}
	else
	{
		//greedy
		int len = U.size();
		if (len == 0)
			return 0;
		int result = 0;
		int qua = numeric_limits<int>::min();
		for (int i = 0; i < len; i++)
		{
			if (quality(U[i], C, S) > qua)
			{
				qua = quality(U[i], C, S);
				result = U[i];
			}
		}
		return result;
	}
}
vector<pair<int, int> > CirMul(vector<pair<int, int> > a, vector<pair<int, int> > b)
{
	vector<pair<int, int> > result;
	int index = 0;
	for (int i = 0; i < b.size(); i++)
	{
		int j = index;
		while (j < a.size() && a[j].second <= b[i].first)
		{
			j++;
		}
		index = j;
		while (j < a.size() && a[j].first <= b[i].second)
		{
			int l = max(a[j].first, b[i].first), m = min(a[j].second, b[i].second);
			if (m - l + 1 >= delta)
			{
				pair<int, int> tmp = make_pair(l, m);
				vector<pair<int, int> > temp;
				temp.push_back(tmp);
				vector<pair<int, int> > tmpresult = Union(result, temp);
				result.assign(tmpresult.begin(), tmpresult.end());
			}
			j++;
		}
	}
	return result;
}
double stability(const vector<int>& a)
{
	if (a.size() > 50)
		return 0;
	vector<pair<int, int> > MIF, T1, T2;
	T1 = timeZon(a);
	vector<pair<int, int> > T = T1;
	while (!T1.empty())
	{
		T2 = CDI(a[0], a, T1, a);
		for (int i = 1; i < a.size(); i++)
		{
			vector<pair<int, int> > T3 = CirMul(T2, CDI(a[i], a, T1, a));
			T2.assign(T3.begin(), T3.end());
			if (T3.size() == 0)
				break;
		}
		vector<pair<int, int> > tmpMIF = Union(MIF, Intersection(T1, T2));
		MIF.assign(tmpMIF.begin(), tmpMIF.end());
		vector<pair<int, int> > T3 = Difference(T2, T1);
		T1.assign(T3.begin(), T3.end());
	}
	return (double)allTimestamps(MIF) / tauV;
}
/**
vector<vector<int> > allCTS(const vector<int>& Gu, const vector<int>& S)
{
	vector<vector<int> > result;
	vector<bool> flag(V, false);
	if (S.empty())
	{
		vector<bool> inGu(V, false);
		for (int i = 0; i < Gu.size(); i++)
		{
			inGu[Gu[i]] = true;
		}
		for (int i = 0; i < Gu.size(); i++)
		{
			if (flag[Gu[i]] == false)
			{
				queue<int> Q;
				vector<int> tmpResult;
				Q.push(Gu[i]);
				flag[Gu[i]] = true;
				while (!Q.empty())
				{
					int now = Q.front();
					Q.pop();
					tmpResult.push_back(now);
					for (int j = 0; j < neighbors[now].size(); j++)
					{
						int nowneighbors = neighbors[now][j].first.first;
						if (inGu[nowneighbors] == true && flag[nowneighbors] == false)
						{
							flag[nowneighbors] = true;
							Q.push(nowneighbors);
						}
					}
				}
				result.push_back(tmpResult);
			}
		}
	}
	else
	{
		vector<bool> mark(V, false);
		vector<int> tmpresult;
		int num = S.size();
		for (int i = 0; i < S.size(); i++)
		{
			mark[S[i]] = true;
		}
		queue<int> Q;
		Q.push(S[0]);
		flag[S[0]] = true;
		mark[S[0]] = false;
		num--;
		while (!Q.empty())
		{
			int now = Q.front();
			Q.pop();
			tmpresult.push_back(now);
			if (mark[now] == true)
			{
				mark[now] = false;
				num--;
			}
			for (int i = 0; i < Gu.size(); i++)
			{
				if (flag[Gu[i]] == false && is_linked(now, Gu[i]) == true)
				{
					flag[Gu[i]] = true;
					Q.push(Gu[i]);
				}
			}
		}
		if (num == 0)
			result.push_back(tmpresult);
	}
	return result;
}
**/
void Pruning(vector<int> U, vector<int> Sel)
{
	cout<<U.size()<<" "<<Sel.size()<<endl;
	gettimeofday(&end_at, 0);
	long long int timeSpent = (end_at.tv_sec - start_at.tv_sec) * 1000 + (end_at.tv_usec - start_at.tv_usec) / 1000;
	if (timeSpent > 86400000) return;
	Count++;
	if (!Sel.empty())
	{
		vector<bool> flag(V, false);
		for (int i = 0; i < U.size(); i++)
		{
			flag[U[i]] = true;
		}
		double dht = (double)(U.size() - 2) / (double)(U.size() - 1);
		if (dht < gama)
		{
			U.clear();
			for (int i = 0; i < Sel.size(); i++)
			{
				vector<int> tmp;
				tmp.push_back(Sel[i]);
				for (int j = 0; j < neighbors[Sel[i]].size(); j++)
				{
					if (flag[neighbors[Sel[i]][j].first.first] == true)
					{
						tmp.push_back(neighbors[Sel[i]][j].first.first);
					}
				}
				countSort(tmp);
				delRepeat(tmp);
				if (i == 0)
					U.assign(tmp.begin(), tmp.end());
				vector<int> tmpU = Intersection(tmp, U);
				U.assign(tmpU.begin(), tmpU.end());
			}
		}
		else
		{
			U.clear();
			for (int i = 0; i < Sel.size(); i++)
			{
				vector<int> tmp;
				tmp.push_back(Sel[i]);
				for (int j = 0; j < neighbors[Sel[i]].size(); j++)
				{
					if (flag[neighbors[Sel[i]][j].first.first] == true)
					{
						tmp.push_back(neighbors[Sel[i]][j].first.first);
					}
				}
				countSort(tmp);
				delRepeat(tmp);
				vector<int> tmptmp;
				for (int j = 0; j < tmp.size(); j++)
				{
					for (int k = 0; k < neighbors[tmp[j]].size(); k++)
					{
						if (flag[neighbors[tmp[j]][k].first.first] == true)
						{
							tmptmp.push_back(neighbors[tmp[j]][k].first.first);
						}
					}
				}
				tmp.insert(tmp.end(), tmptmp.begin(), tmptmp.end());
				countSort(tmp);
				delRepeat(tmp);
				if (i == 0)
					U.assign(tmp.begin(), tmp.end());
				vector<int> tmpU = Intersection(tmp, U);
				U.assign(tmpU.begin(), tmpU.end());
			}
		}
	}
	vector<pair<int, int> > tauS = timeZon(Sel);
	double nowub = ub(U, Sel);
	nowlb = lb(U, Sel);
	if ((double)U.size() < M)
	{
		U.clear();
		return;
	}
	if (nowub < nowlb || nowub < M || (M <= nowub && nowub < (double)Sel.size()))
	{
		U.clear();
		return;
	}
	if (M <= nowub && nowub == (double)Sel.size())
	{
		U.clear();
		U.assign(Sel.begin(), Sel.end());
		return;
	}

	if (ub(U, Sel) < M)
		return;
	if (!(Difference(Sel, U).empty()))
		return;
	if (stability(U) >= ro)
	{
		vector<vector<int> > tmpCTS;
		tmpCTS.push_back(U);
		vector<vector<int> > tmpds = Union(ds, tmpCTS);
		ds.assign(tmpds.begin(), tmpds.end());
	}
	else
	{
		//写入时间
		int tmpGama = 100 * gama, tmpRo = 100 * ro;
		string tmpFile = "tmp." + Str + to_string(tmpGama) + "_" + to_string(tmpRo) + "_" + to_string((int)M);
		ofstream out(tmpFile, ios::app);
		gettimeofday(&end_at, 0);
		long long int curTimeSpent = (end_at.tv_sec - start_at.tv_sec) * 1000 + (end_at.tv_usec - start_at.tv_usec) / 1000;
		out << "Current Time Spent:" << curTimeSpent << "ms" << endl;
		out.close();


		vector<pair<int, int>> timeZoC = timeZon(U);
		currentTimezone.assign(timeZoC.begin(), timeZoC.end());
		fill(currentGraph.begin(), currentGraph.end(), false);
		for (int j = 0; j < U.size(); j++)
		{
			currentGraph[U[j]] = true;
			vector<pair<int, int> > Iu = CDI(U[j], U, timeZoC, Sel);
			StoreCDI[U[j]].assign(Iu.begin(), Iu.end());
		}
		vector<int> tmpCTS = Difference(U, Sel);
		int node = qualitySort(tmpCTS, U, Sel);
		vector<int> D;
		bool flagB = true;
		vector<int> tmpNodeV{ node };
		vector<int> CdifV = Difference(U, tmpNodeV);
		D.push_back(node);
		vector<pair<int, int>> timeZo = timeZon(CdifV);
		vector<int> mark(V, 2);
		vector<bool> flag(V, false);
		for (int j = 0; j < CdifV.size(); ++j) {
			mark[CdifV[j]] = 1;
			//cout<<CdifV[j]<<" ";
		}
		//cout<<endl;
		for (int j = 0; j < Sel.size(); j++)
		{
			flag[Sel[j]] = true;
		}
		for (int j = 0; j < CdifV.size(); ++j)
		{
			vector<pair<int, int> > Iu = CDI(CdifV[j], CdifV, timeZo, Sel);
			StoreCDI[CdifV[j]].assign(Iu.begin(), Iu.end());
			if (allTimestamps(Iu) < ro * tauV)
			{
				if (flag[CdifV[j]] == true)
				{
					flagB = false;
					break;
				}
				D.push_back(CdifV[j]);
			}
		}
		if (mode != 2) {
			Pruning(CdifV, Sel);
		}else if (flagB == true && mode ==2)
			Pruning(Difference(CdifV, D), Sel);
		else
			Pruning(CdifV, Sel);
		vector<int> A, SUnionv = Union(Sel, tmpNodeV);
		if (mode > 2)
		{
			for (int j = 0; j < SUnionv.size(); ++j)
			{
				vector<pair<int, int>> tmpTimeZo = CirMul(timeZo, CDI(SUnionv[j], U, timeZoC, SUnionv));
				timeZo.assign(tmpTimeZo.begin(), tmpTimeZo.end());
				if (timeZo.empty())
					break;
			}
			if (allTimestamps(timeZo) < ro * tauV)
				return;
			vector<int> CDifSUniv = Difference(U, SUnionv);
			for (int j = 0; j < CDifSUniv.size(); ++j)
			{
				vector<int> tmpSUniu{ CDifSUniv[j] };
				if ((double)allTimestamps(CirMul(timeZo, CDI(CDifSUniv[j], U, timeZoC, Union(SUnionv, tmpSUniu)))) < ro * tauV)
					A.push_back(CDifSUniv[j]);
			}
		}
		Pruning(Difference(CdifV, A), SUnionv);
		
	}


}


void RTGM()
{
	string filename = "result." + Str;
	gettimeofday(&start_at, 0);
	vector<int> lastV;
	vector<int> Sel;
	for (int i = 0; i < V; i++)
	{
		if (!neighbors[i].empty())
		{
			lastV.push_back(i);
		}
	}
	Pruning(lastV, Sel);
	gettimeofday(&end_at, 0);
	ofstream out(filename.c_str(), ios::app);
	out << endl;
	out << endl;
	out << Str << " gama:" << gama << " rho:" << ro << " theta:" << M << " mode:" << mode << endl;
	out << "reading time: "
		<< (end_at.tv_sec - start_at.tv_sec) * 1000 +
		(end_at.tv_usec - start_at.tv_usec) / 1000
		<< "ms" << endl;
	if (ds.empty())
		out << "not exist!!!" << endl;
	else
	{
		out << "Count:" << ds.size() << endl;
		for (int i = 0; i < ds.size(); i++)
		{
			for (int j = 0; j < ds[i].size(); j++)
			{
				out << ds[i][j] << " ";
			}
			out << endl;
		}
	}
	out.close();
}

};

int main(int argc, char* argv[])
{
	BBPre pre = new BBPre();
	pre.Str = argv[1];
	ifstream in1(argv[1]);
	double gama = atof(argv[2]);
	int delta = 1;
	double ro = atof(argv[4]);
	int c = 2;                              //按月采样
	double M = atof(argv[5]);
	pre.init(in1, gama, delta, ro, c, M);
	pre.prePrune();
	pre.outPut();

	BBProcess process = new BBProcess();
	process.Str = argv[1];
	ifstream in2(preOutPutFileName);
	process.mode = atoi(argv[3]);
	process.tauV = pre.tauV;
	process.init(in2, gama, ro, M);
	process.RTGM();

	return 0;
}
