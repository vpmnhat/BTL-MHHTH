#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

#define ll          long long
#define ull         unsigned long long
#define for_(n)     for(ll i=0; i<n; i++)
#define for__(a,b)  for(ll i=a; i<b; i++)
#define _for(i,a,b) for(ll i=a; i<b; i++)
#define mp          make_pair
#define fi          first
#define se          second
#define pb          push_back
#define pii         pair<int, int>
#define pll         pair<long long, long long>
#define el          "\n"
#define debug(x)    cerr << "[debug] " << #x << " : " << x << endl;
#define TPPL        mbdrn(ansLB, allUB);

using namespace __gnu_pbds;
using namespace __gnu_cxx;
using namespace std;

typedef tree<ll,null_type,less<ll>,rb_tree_tag,tree_order_statistics_node_update> ordered_set;
const long long MOD=1000000007;

template<class T, class P>
ostream& operator<<(ostream& os, const pair<T,P> &a) { os << "{" << a.first << ";" << a.second << "}"; return os; }
template<class T>
ostream& operator<<(ostream& os, const vector<T> &a) { ; for(auto it: a) os << it << " "; ; return os; }
template<class T>
ostream& operator<<(ostream& os, const deque<T> &a) { ; for(auto it: a) os << it << " "; ; return os; }
template<class T>
ostream& operator<<(ostream& os, const set<T> &a) { ; for(auto it: a) os << it << " "; ; return os; }
template<class T>
ostream& operator<<(ostream& os, const multiset<T> &a) { ; for(auto it: a) os << it << " "; ; return os; }

ll gcd(ll a, ll b) { return b==0? a : gcd(b,a%b); }
ll lcm(ll a, ll b) { return a/(gcd(a,b))*b; }

ll pow_mod(ll x, ll y, ll mod) { //mod<3.10^9
    ll ans=1;
    while(y>0) {
        if(y%2==1) ans=ans*x%mod;
        y=y>>1;
        x=x*x%mod;
    }
    return ans%mod;
}

void mbdrn(float &a, float &b) { if(b<a) b = max(b, a+(float)0.03 * a); }

struct Edge
{
    Edge(int _a, int _b, int _c, int _f, ll _w) {
        a = _a; b = _b; c = _c; f = _f; w = _w;
    }
    ~Edge() { };
    int a; //from
    int b; //to
    int c; //capacity
    int f; //flow
    ll w; //weight
    Edge* r;
};

struct EdgeT
{
    EdgeT(int a_, int b_): a(a_), b(b_) { f.assign(151,0); }
    ~EdgeT() { };
    int a; //from
    int b; //to
    int c[151]; //capacity
    vector<int> f; //flow
    ll w[151]; //weight
    EdgeT* r;
};


const int MAX_NODES = 2000;
const int MAX_DIST = 2000000; //do not choose INT_MAX because you are adding weights to it
vector<Edge*> adj[MAX_NODES];
vector<EdgeT*> adjT[MAX_NODES];
ll distances[MAX_NODES];
ll LG[MAX_NODES];
Edge* parents[MAX_NODES];
EdgeT* prT[MAX_NODES];
float pen[51][51];
float gx[51][51];
int mu[51][51], mc[51][51], x[51][51];
int u[51][151][51][101], y[51][151][51][51];
float alp[51][151][51][51], c[51][151][51][51], g[51][151][51][51];

int nodecount, nS;
int nK;
int source, sink, T, TT;
float step = 0.001, ansLB = 0;

void iniAlpha(int T) {
    for(int i=0; i<nodecount; i++) {
        for(int j=0; j<nodecount; j++) {
            for(int s = 0; s < nS; s++) {
                for(int t=0; t<T; t++) alp[s][t][i][j]=0.04;
            }
        }
    }
}

bool find_path(int from, int to, vector<Edge*>& output)
{
    fill(distances, distances+nodecount, 1e12);
    fill(parents, parents+nodecount, (Edge*)0);
    distances[from] = 0;

    bool updated = true;
    int nt = nodecount;
    while(updated && nt>0)
    {
        nt--;
        updated = false;
        for(int j = 0; j < nodecount; ++j)
            for(int k = 0; k < (int)adj[j].size(); ++k){
                Edge* e = adj[j][k];
                if( e->f >= e->c ) continue;
                if( distances[e->b] > distances[e->a] + e->w )
                {
                    distances[e->b] = distances[e->a] + e->w;
                    parents[e->b] = e;
                    updated = true;
                }
            }
    }


    output.clear();
    if(distances[to] == (ll)1e12) return false;
    int cur = to;

    while(parents[cur])
    {
        output.push_back(parents[cur]);
        cur = parents[cur]->a;
    }
    return true;
}


float min_cost_max_flow(int source, int sink, int K)
{
    ll total_cost = 0;
    int tt_flow = 0;
    vector<Edge*> p;
    while(find_path(source, sink, p))
    {
        int flow = INT_MAX;
        for(int i = 0; i < p.size(); ++i)
            if(p[i]->c - p[i]->f < flow) flow = p[i]->c - p[i]->f;

        if( tt_flow + flow > K) flow = K - tt_flow;
        ll cost = 0;
        for(int i = 0; i < p.size(); ++i) {
            cost += p[i]->w;
            p[i]->f += flow;
            p[i]->r->f -= flow;
        }
        cost = cost * flow; //cost per flow
        total_cost += cost;
        tt_flow += flow;
        if(tt_flow >= K) break;
    }
    return (float)total_cost/1000.0;
}

bool find_pathT(int from, int to, vector<EdgeT*>& output)
{
    ll mxD = 200000000;
    fill(LG, LG+nodecount, mxD);
    fill(prT, prT+nodecount, (EdgeT*)0);
    LG[from] = 0;
    bool updated = true;
    while(updated)
    {
        updated = false;
        for(int j = 0; j < nodecount; ++j)
            for(int k = 0; k < (int)adjT[j].size(); ++k){
                EdgeT* e = adjT[j][k];
                if( LG[j] >= T || e->f[LG[j]] >= e->c[LG[j]] ) continue;
                if( LG[e->b] > LG[e->a] + e->w[LG[j]] )
                {
                    LG[e->b] = LG[e->a] + e->w[LG[e->a]];
                    prT[e->b] = e;
                    updated = true;
                }
            }
    }
    output.clear();
    if(LG[to] == mxD) return false;
    int cur = to;

    int it = 0;
    while(prT[cur])
    {
        if(it++ > nodecount) return false;
        output.push_back(prT[cur]);
        cur = prT[cur]->a;
    }
    reverse(output.begin(), output.end());
    return true;
}


float min_cost_max_flow_T(int source, int sink, int K)
{
    ll total_cost = 0;
    int tt_flow = 0;

    vector<EdgeT*> p;
    while(find_pathT(source, sink, p))
    {
        int flow = INT_MAX;
        ll t = 0;
        for(int i = 0; i < p.size(); t+= p[i]->w[t], ++i )
            if(p[i]->c[t] - p[i]->f[t] < flow) flow = p[i]->c[t] - p[i]->f[t];  ////

        if( tt_flow + flow > K) flow = K - tt_flow;
        ll cost = 0;
        t = 0;
        for(int i = 0; i < p.size();  t+= p[i]->w[t], ++i) {
            cost += p[i]->w[t];
            p[i]->f[t] += flow;
            p[i]->r->f[t] -= flow;
        }
        cost *= flow; //cost per flow
        total_cost += cost;
        tt_flow += flow;

        if(tt_flow >= K) break;
    }
    return (float)total_cost;
}

void add_edge_T(int a, int b, int s, vector<float> &prob ) {
    EdgeT * e = new EdgeT(a,b);
    EdgeT * re = new EdgeT(b,a);
    /*
    float sum_g =0 ;
    for(int t=0; t<T; t++) {
        sum_g += prob[s] * c[s][t][a][b];
        if(t<=TT) sum_g += alp[s][t][a][b];
    }
    */

    for(int t=0; t<T; t++) {
        if(t<=TT) {
            e ->w[t] = (int) ( ( prob[s] * c[s][t][a][b] + alp[s][t][a][b] ) ) ;
            re -> w[t] = - e ->w[t];
        }
        else {
            e ->w[t] = (int) (  prob[s] * c[s][t][a][b] );
            re -> w[t] = - e ->w[t];
        }
    }

    for(int t=0; t<T; t++) {
        e -> c[t] = u[s][t][a][b];
        re -> c[t] = 0;
    }

    e -> r = re;
    re -> r = e;

    adjT[a].push_back(e);
    adjT[b].push_back(re);
}


void add_edge(int a, int b, int c, float w)
{
    Edge* e = new Edge(a, b, c, 0, (int)(w*1000) );
    Edge* re = new Edge(b, a, 0, 0, -(int)(w*1000));
    e->r = re;
    re->r = e;
    adj[a].push_back(e);
    adj[b].push_back(re);
}

void fetch_y(int s) {
    for(int i=0; i<nodecount; i++) {
        for(auto e : adjT[i]) {
            if(e->a > e->b) continue;
            for(int t=0; t<T; t++) {
                y[s][t][e->a][e->b] = e->f[t];
            }
        }
    }
}

void fetch_x() {
    for(int i=0; i<nodecount; i++) {
        for(auto e: adj[i]) {
            if(e->a > e->b) continue;
            x[i][e->b] = e->f;
        }
    }
}

void updateAlp() {
    for(int i=0; i<nodecount; i++) {
        for(int j=0; j<nodecount; j++) {
            for(int s = 0; s < nS; s++) {
                for(int t=0; t<=TT; t++) alp[s][t][i][j] =  alp[s][t][i][j] + step * ( y[s][t][i][j] - x[i][j] ) ;
            }
        }
    }
}

void dbe() {
    for(int i=0; i<nodecount; i++) {
        cout << i << " :";
        for(auto e: adj[i]) {
            cout << e->b << "_" << e->w <<  " ";
        }
        cout << endl;
    }
}

void solve()
{

    int nIter = 1;
    cin >> nodecount;
    cin >> source >> sink >> nK;
    cin >> T >> TT;
    cin >> nS;
    vector<float> prob(nS);
    for(int i=0; i<nS; i++) {
        cin >> prob[i];
    }
    cin >> nIter;

    iniAlpha(T);

    for(int i=0; i<nodecount; i++) {
        if(i%5!=4) {
            int j=i+1;
            float sumU = 0, sumC = 0;
            for(int s=0; s < nS; s++) {
                int val = rand()%20+80;
                int cot = rand()%3+1;
                for(int t=0; t<T; t++) {
                    u[s][t][i][j]= val;
                    c[s][t][i][j]= cot;
                }
                sumU += (float)val * prob[s];
                sumC += (float)cot * prob[s];
            }
            mu[i][j] = (int)sumU;
            mc[i][j] = (int)sumC;
            pen[i][j] = 3;
        }
        if(i<nodecount-5) {
            int j=i+5;
            float sumU = 0, sumC = 0;
            for(int s=0; s < nS; s++) {
                int val = rand()%20+80;
                int cot = rand()%5+1;
                for(int t=0; t<T; t++) {
                    u[s][t][i][j]= val;
                    c[s][t][i][j]= cot;
                }
                sumU += (float)val * prob[s];
                sumC += (float)cot * prob[s];
            }
            mu[i][j] = (int)sumU;
            mc[i][j] = (int)sumC;
            pen[i][j] = 3;
        }
    }

    float allUB = 10000000;
    while(nIter--) {
        //cout << "#######" << endl;
        for(int i=0; i<nodecount; i++) {
            if(i%5!=4) {
                int j = i+1;
                float ss = 0;
                for(int s=0; s<nS; s++) {
                    for(int t=0; t<=TT; t++) ss += alp[s][t][i][j];
                }

                add_edge(i,j,mu[i][j], pen[i][j] - ss );
            }
            if(i<nodecount-5) {
                int j=i+5;

                float ss = 0;
                for(int s=0; s<nS; s++) {
                    for(int t=0; t<=TT; t++) ss += alp[s][t][i][j];
                }

                add_edge(i,j,mu[i][j], pen[i][j] - ss );
            }
        }


        //dbe();

        //cout <<"X" << endl;
        float ansSP1 = min_cost_max_flow(source, sink, nK);
        float ansSP2 = 0;



        fetch_x();

        for(int s=0; s<nS; s++) {
            float ta = 0;
            for(int i=0; i<nodecount; i++) {
                if(i%5!=4) {
                    int j=i+1;
                    add_edge_T(i,j,s,prob);
                }
                if(i<nodecount-5) {
                    int j=i+5;
                    add_edge_T(i,j,s,prob);
                }
            }
            //cout << "DD" << endl;
            ta = min_cost_max_flow_T(source,sink,nK);
            //debug(ta);
            ansSP2+= ta;
            fetch_y(s);
            for(int i=0; i<nodecount; i++) adjT[i].clear();
        }

        //debug(ansSP1);
        //debug(ansSP2);

        ansLB = max(ansLB, ansSP1 + ansSP2);

        //debug(ansLB);


        for(int i = 0; i < nodecount; ++i){
            for(unsigned int j = 0; j < adj[i].size(); ++j)
                delete adj[i][j];
            adj[i].clear();
        }

        updateAlp();

        float ub = 0;
        for(int i=0; i<nodecount; i++) {
            if(i%5!=4) {
                int j = i+1;
                for(int s=0; s<nS; s++) {
                    for(int t=0; t<TT; t++) ub += prob[s] * c[s][t][i][j] * y[s][t][i][j];
                }
                ub += pen[i][j] * x[i][j];
            }
            if(i<nodecount-5) {
                int j=i+5;

                for(int s=0; s<nS; s++) {
                    for(int t=0; t<TT; t++) ub += prob[s] * c[s][t][i][j] * y[s][t][i][j];
                }
                ub += (pen[i][j] ) * x[i][j];
            }
        }

        allUB = min(allUB, ub);

    }

    TPPL;
    cout << "Time Thresholds: " << TT << endl;
    cout << "Lower Bound: " << ansLB << el;
    cout << "Upper Bound: "  <<allUB << el;


    cout << "Prior plan: " << endl;
    for(int i=0; i<nodecount; i++) {
            if(i%5!=4) {
                int j = i+1;
                cout << i << " -> " << j << " : " << x[i][j] << el;
            }
            if(i<nodecount-5) {
                int j=i+5;
                cout << i << " -> " << j << " : " << x[i][j] << el;
            }
        }

}

int main() {
    //clock_t begin=clock();
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    freopen("INP.inp","r",stdin);
    //freopen( "OUT.out", "w", stdout );

    int t=1;
    //cin >> t;
    while(t--)
        solve();

    //cout << "TIME : " << (double)(clock()-begin)/CLOCKS_PER_SEC << "s.";
    return 0;
}


