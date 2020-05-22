#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <random>
#include <ctime>
#include <fstream>
#include "IntervalTree.h"
#include "Matrix.h"

using namespace std;

#ifndef GRAPH_H
#define GRAPH_H

struct Rect{
    int i1, i2, k1, k2, h1, h2;
    Rect(int i1, int i2, int k1, int k2, int h1, int h2): i1(i1), i2(i2), k1(k1), k2(k2), h1(h1), h2(h2){}
};


class Graph{

public:

    struct Vertex
    {
        int x, y;
        bool d, u;
        bool r, l;
        Vertex(int x, int y, bool d, bool u, bool r, bool l): x(x), y(y), d(d), u(u), r(r), l(l){}
        Vertex(const Vertex& v){
            x=v.x;y=v.y;d=v.d;u=v.u;r=v.r;l=v.l;
        }
        bool operator ==(Vertex const& other) {
            return (x == other.x && y == other.y);
        }

        bool isConcave(){
            if (d + u + r + l != 2) return false;
            if (d + u == 2) return false;
            if (r + l == 2) return false;
            return true;
        }
    };
    
    struct Edge
    {
        //int u, v; // vertex index
        int uid, vid;
        Vertex* u;
        Vertex* v;
        Edge(Vertex* u, Vertex* v, int uid, int vid): u(u), v(v), uid(uid), vid(vid){}
        Edge(const Edge& e){
            u=e.u;v=e.v;uid=e.uid;vid=e.vid;
        }
        Vertex* other(Vertex* o){
            if (o == v) return u;
            return v;
        }

        int other_id(Vertex* o){
            if(o == v) return uid;
            return vid;
        }
        
        bool operator<(const Edge& other) const {
            if (u -> y == other.u -> y){
                return u -> x < other.u -> x;   
            }
            return u -> y < other.u -> y;
        }

        bool has(Vertex* o){
            return (o == u || o == v);
        }

        int getEnd(Vertex* o){
            // cout <<"getEnd: " << uid << " " << vid << endl;
            if(o -> x == u -> x) return uid;
            if (o -> x == v -> x) return vid;
            return -1;
        }


        bool isDown(Vertex * a){
            return (a == u && u -> x < v -> x || a == v && u -> x > v -> x);
        }
        
        bool isUp(Vertex * a){
            return (a == u && u -> x > v -> x || a == v && u -> x < v -> x);
        }

        bool isLeft(Vertex * a){
            return (a == u && u -> y > v -> y || a == v && u -> y < v -> y);
        }
        
        bool isRight(Vertex * a){
            return (a == u && u -> y < v -> y || a == v && u -> y > v -> y);
        }


        void print(){
            // cout <<"<(" << u -> x << ", " << u -> y << ") - (" << v -> x << ", " << v -> y << ")>" << endl;
        }
        
    };

    vector<Vertex*> vertexes;
    vector<Vertex*> true_vertexes;
    vector<vector<int>> graph;
    vector<Edge*> edges;
    vector<vector<int>> c_v;
    vector<vector<int>> r_v;
    vector<vector<int>> true_c_v;
    vector<vector<int>> true_r_v;
    Matrix a;
    vector<Rect*> rects;

    void update(vector<Vertex*>& vertexes, vector<vector<int>>& c_v, vector<vector<int>>& r_v, int i, int j, bool d, bool u, bool r, bool l){
        int id = vertexes.size();
        vertexes.push_back((new Vertex(i, j, d, u, r, l)));
        c_v[j].push_back(id);
        r_v[i].push_back(id);
    }
    //int** g;

    bool compareEdge(const Edge& u, const Edge& v){
        return u.u -> y < v.u -> y;
    }

    void add_edge(Vertex* u, Vertex* v, int uid, int vid){
        int id = edges.size();
        edges.push_back((new Edge(u, v, uid, vid)));
        // cout <<"push edge: " << uid << " " << vid << endl;
        graph[uid].push_back(id);
        graph[vid].push_back(id);
    }

    void generate_graph(vector<Vertex*>& vertexes, vector<vector<int>>& c_v, vector<vector<int>>& r_v){
        for (int i = 0; i < r_v.size(); i++) {
            int s = r_v[i].size();
            //cout <<s - 1 << endl;
            for(int t = 0; t < s-1; t++){
                int u = r_v[i][t];
                //cout <<i << " " << t << endl;
                int v = r_v[i][t + 1];
                //cout <<i << " " << t << endl;
                if (vertexes[u] -> r == 1){
                    if(vertexes[v] -> l == 1) {
                        add_edge(vertexes[u], vertexes[v], u, v);
                    } else {
                        cout <<"Something Wrong" << endl;
                    }
                }
                //cout <<i << " " << t << endl;

            }
        }
        //cout <<"hello?" << endl;
        for (int i = 0; i < c_v.size(); i++) {
            int s = c_v[i].size();
            //cout <<s - 1 << endl;
            for(int t = 0; t < s - 1; t++){
                int u = c_v[i][t];
                int v = c_v[i][t + 1];
                //cout <<i << "~" << t << endl;
                if (vertexes[u] -> d == 1){
                    if(vertexes[v] -> u == 1) {
                        add_edge(vertexes[u], vertexes[v], u, v);
                    } else {
                        cout <<"Something Wrong" << endl;
                    }
                }
                //cout <<i << "~" << t << endl;
            }
        }
    }

    Edge* find_dir(int i, int dir){ // d, u, r, l
        int size = graph[i].size();
        Vertex* v = vertexes[i];
        for(int j = 0; j < size; j++){
            Edge* e = edges[graph[i][j]];
            if (dir == 0 && e -> isDown(v)) return e;
            if (dir == 1 && e -> isUp(v)) return e;
            if (dir == 2 && e -> isRight(v)) return e;
            if (dir == 3 && e -> isLeft(v)) return e;
        }
        return NULL;
    }


    Graph(const Matrix& m){
        //g = (int **)calloc(m.rows + 2, sizeof(int*));
        //for(int i = 0; i < m.rows + 2; i++) g[i] = (int *)calloc(m.cols + 2, sizeof(int));

        //for(int i = 0; i <= m.rows; i++) g[i][0] = g[i][m.cols] = 1;
        //for(int i = 0; i <= m.cols; i++) g[0][i] = g[m.rows][i] = 1;
        //g[m.rows][0] = g[0][m.cols] = g[m.rows][m.cols] = g[0][0] = 1;
        //a = m;
        
        c_v.reserve(m.cols + 1);
        r_v.reserve(m.rows + 1);
        for (int i = 0; i <= m.cols; i++) c_v.push_back(vector<int>());
        for (int i = 0; i <= m.rows; i++) r_v.push_back(vector<int>());
        true_c_v.reserve(m.cols + 1);
        true_r_v.reserve(m.rows + 1);
        for (int i = 0; i <= m.cols; i++) true_c_v.push_back(vector<int>());
        for (int i = 0; i <= m.rows; i++) true_r_v.push_back(vector<int>());
        
        
        update(vertexes, c_v, r_v, 0, 0, 1, 0, 1, 0);

        for(int j = 0; j < m.cols - 1; j++) { // firts row vertexes
            if (m.data[0][j] != m.data[0][j + 1]) {
                update(vertexes, c_v, r_v, 0, j + 1, 1, 0, 1, 1);
            }
        }
        //cout <<"hello?" << endl;
        update(vertexes, c_v, r_v, 0, m.cols, 1, 0, 0, 1);
        // cout <<"hello?" << endl;

        for(int i = 0; i < m.rows - 1; i++) {
            if (m.data[i][0] != m.data[i + 1][0]) {
                update(vertexes, c_v, r_v, i + 1, 0, 1, 1, 1, 0);
            }
            for(int j = 0; j < m.cols - 1; j++) {
                bool d, u, r, l;
                d = u = r = l = 1;
                if (m.data[i][j] == m.data[i + 1][j]) l = 0;
                if (m.data[i][j] == m.data[i][j + 1]) u = 0;
                if (m.data[i + 1][j] == m.data[i + 1][j + 1]) d = 0;
                if (m.data[i][j + 1] == m.data[i + 1][j + 1]) r = 0;

                int a = d + u, b = r + l;
                if ((a == 2 && b == 0) || (a == 0 && b == 2) || a + b == 0) continue;

                update(vertexes, c_v, r_v, i + 1, j + 1, d, u, r, l);
            }
            if (m.data[i][m.cols - 1] != m.data[i + 1][m.cols - 1]) {
                update(vertexes, c_v, r_v, i + 1, m.cols, 1, 1, 0, 1);
            }
        }

        update(vertexes, c_v, r_v, m.rows, 0, 0, 1, 1, 0);
        for(int j = 0; j < m.cols - 1; j++) { // firts row vertexes
            if (m.data[m.rows -1][j] != m.data[m.rows-1][j + 1]) {
                update(vertexes, c_v, r_v, m.rows, j + 1, 0, 1, 1, 1);
            }
        }
        update(vertexes, c_v, r_v, m.rows, m.cols, 0, 1, 0, 1);

        // cout <<"hello?" << endl;
        print_rc(c_v, r_v);

        int size = vertexes.size();

        graph.reserve(size);
        for (int i = 0; i < size; i++) graph.push_back(vector<int>());
        //edges.reserve(10);
        // cout <<"hello?" << endl;
        
        generate_graph(vertexes, c_v, r_v);
        
        //cout <<"hello?" << endl;
        print_edges();
        print_graph();
        //cout <<"hello?!" << endl;


        // sweep line
        // define new vertexes and graph and cv rv group
        // std::set
        // for each concave vertex, determine, a b c
        
        set<Edge> active;
        int p = 0;
        for(int r = 0; r <= m.rows; r++) {
            // cout <<"r " << r << ":";
            // first insert all the vertical edge
            for(int i = p; i < size; i++){
                Vertex* curv = vertexes[i];
                if(curv -> x != r) {break;}
                // cout <<" i:" << i;
                int e_size = graph[i].size();
                for (int j = 0; j < e_size; j++){
                    Edge* e = edges[graph[i][j]];
                    if (e -> isDown(curv)){
                        active.insert(*e);
                        // cout<< " ins: <(" << e -> u -> x << ", " << e -> u -> y << ") - (" << e -> v -> x << ", " << e -> v -> y << ")> ";
                    }
                }
                // cout <<endl;
            }

            // iter again for sweep line
            for(int i = p; i < size; i++){
                Vertex* curv = vertexes[i];
                if(curv -> x != r) {p = i; break;}
                // cout <<" i:" << i;
                int e_size = graph[i].size();
                bool hasUpdate = false;

                // if concave
                if (curv -> isConcave() && !(curv -> x == 0 || curv -> x == m.rows || curv -> y == 0 || curv -> y == m.cols)) {

                    for (int j = 0; j < e_size; j++){
                        Edge* e = edges[graph[i][j]];
                        // a down edge
                        if (e -> isLeft(curv)){
                            set<Edge>::iterator up;
                            for (int t = 0; t < e_size; t++){
                                Edge* ee = edges[graph[i][t]];
                                if (ee -> isUp(curv) || ee -> isDown(curv)){
                                    // cout<< " **ud: <(" << ee -> u -> x << ", " << ee -> u -> y << ") - (" << ee -> v -> x << ", " << ee -> v -> y << ")> ";
                                    up = active.upper_bound(*ee);
                                }
                            }
                            Edge upper = *up;
                            int vid = upper.getEnd(curv);
                            curv -> r = 1;
                            update(true_vertexes, true_c_v, true_r_v, curv -> x, curv -> y, curv-> d, curv -> u, curv -> r, curv -> l);
                            if (vid != -1){
                                // cout <<"vid: " << vid << endl;
                                vertexes[vid] -> l = 1; 
                            } else {
                                update(true_vertexes, true_c_v, true_r_v, curv -> x, upper.u -> y, 1, 1, 0, 1);
                            }

                            // cout<< " left: <(" << e -> u -> x << ", " << e -> u -> y << ") - (" << e -> v -> x << ", " << e -> v -> y << ")> ";
                            // cout<< " __ upper: <(" << (upper.u) -> x << ", " << (upper.u) -> y  << ") - (" << (upper.v) -> x << ", " << (upper.v) -> y << ")> ";
                            break;
                        }
                        if (e -> isRight(curv)){
                            set<Edge>::iterator low;
                            for (int t = 0; t < e_size; t++){
                                Edge* ee = edges[graph[i][t]];
                                if (ee -> isUp(curv) || ee -> isDown(curv)){
                                    // cout<< " **ud: <(" << ee -> u -> x << ", " << ee -> u -> y << ") - (" << ee -> v -> x << ", " << ee -> v -> y << ")> ";
                                    low = active.lower_bound(*ee);
                                    break;
                                }
                            }
                            Edge lower = *--low;
                            int vid = lower.getEnd(curv);
                            curv -> l = 1;
                            if (vid != -1){
                                true_vertexes[true_vertexes.size() - 1] -> r = 1; 
                                // cout <<"lower vid: " << vid << endl;

                            } else {
                                int x = curv -> x;
                                int y = lower.u -> y;
                                Vertex* last_v = true_vertexes.back();
                                if (last_v -> x == x && last_v -> y == y){
                                    last_v -> r = 1;
                                } else {
                                    update(true_vertexes, true_c_v, true_r_v, curv -> x, lower.u -> y, 1, 1, 1, 0); 
                                }
                            }   
                            update(true_vertexes, true_c_v, true_r_v, curv -> x, curv -> y, curv-> d, curv -> u, curv -> r, curv -> l);
                            // cout <<"true: " << true_vertexes.size() << curv -> x << curv -> y << endl;

                            // cout<< " right: <(" << e -> u -> x << ", " << e -> u -> y << ") - (" << e -> v -> x << ", " << e -> v -> y << ")> ";
                            // cout<< " __ lower: <(" << (lower.u) -> x << ", " << (lower.u) -> y  << ") - (" << (lower.v) -> x << ", " << (lower.v) -> y << ")> ";
                            break;
                        }
                        
                    }
                } else {
                    update(true_vertexes, true_c_v, true_r_v, curv -> x, curv -> y, curv-> d, curv -> u, curv -> r, curv -> l);
                }
                
                for (int j = 0; j < e_size; j++){
                    Edge* e = edges[graph[i][j]];

                    // a down edge
                    if (e -> isUp(curv)){
                        active.erase(*e);
                        // active.insert(*e);
                        // cout<< " erase: <(" << e -> u -> x << ", " << e -> u -> y << ") - (" << e -> v -> x << ", " << e -> v -> y << ")> ";
                        break;
                    }
                }
                // cout <<endl;
            }
            // cout <<"set: "<< endl;
            for (auto it = active.begin(); it != active.end(); it++) {
                Edge e = *it;
                e.print();
            }

            // cout <<endl;
            
            //if (curv -> x == 0 || curv -> x == m.rows || curv -> y == 0 || curv -> y == m.cols ) continue; // on the outside
        
        }
        // cout <<"set len:" << active.size() << endl;
        
        size = true_vertexes.size();
        edges.clear();
        graph.clear();
        graph.reserve(size);
        for (int i = 0; i < size; i++) graph.push_back(vector<int>());
        
        // cout <<" =========== " << endl;
        generate_graph(true_vertexes, true_c_v, true_r_v);
        vertexes = true_vertexes;
        c_v = true_c_v;
        r_v = true_r_v;
        // print_edges();
        // print_graph();

        // find rects
        
        for(int i = 0; i < size; i++) {
            Vertex* v = vertexes[i];
            if (v -> d == 1 && v -> r == 1){
                // cout <<"Rv: " << v -> x << ", " << v -> y << endl;
                int i1 = v -> x;
                int i2 = v -> y;
                int k1 = 0, k2 = 0;

                // cout <<"!!" << endl;


                Edge* down_edge = find_dir(i, 0);
                Vertex* down_v = down_edge -> other(v);
                int down_id = down_edge -> other_id(v);
                // cout <<"de: " << down_v -> x << " " << down_v -> y << " " << down_id << endl;

                Edge* right_edge = find_dir(down_id, 2);
                while(right_edge == NULL) {
                    Edge* nd_edge = find_dir(down_id, 0);
                    Vertex* nd_v = nd_edge -> other(down_v);
                    int nd_id = nd_edge -> other_id(down_v);
                    right_edge = find_dir(nd_id, 2);
                    down_id = nd_id;
                    down_v = nd_v;
                }


                // cout <<(right_edge == NULL) << endl;
                Vertex* r_v = right_edge -> other(down_v);
                int r_id = right_edge -> other_id(down_v);
                // cout <<"rv: " << r_v -> x << " " << r_v -> y << endl;

                
                Edge* up_edge = find_dir(r_id, 1);
                while(up_edge == NULL) {
                    Edge* nr_edge = find_dir(r_id, 2);
                    Vertex* nr_v = nr_edge -> other(r_v);
                    int nr_id = nr_edge -> other_id(r_v);
                    up_edge = find_dir(nr_id, 1);
                    r_id = nr_id;
                    r_v = nr_v;
                }
                // cout <<"!!" << endl;

                k1 = r_v -> x;
                k2 = r_v -> y;

                // cout <<k1 << k2 << endl;
                int h1 = 0;
                int h2 = m.data[i1][i2];
                rects.push_back(new Rect(i1, k1 - 1, i2, k2 - 1, h1, h2));
            }
        }
        for (int i = 0; i < rects.size(); i++){
            Rect* r = rects[i];
            // cout <<"Rect " << i << ": (" << r -> i1 << "," << r -> i2 << "," << r -> k1 << "," << r -> k2 << "), " << r -> h1 << "-" << r -> h2 << endl;
        }
        /**/
        // cout << "yey" << endl;

        

    };


    ~Graph(){
        edges.clear();
        vertexes.clear();
        graph.clear();
        r_v.clear();
        c_v.clear();
    };


    

    void print_graph() {
        for (int i = 0; i < graph.size(); i++) {
            // cout <<"v_" << i << ". ("<<  vertexes[i] -> x << ",  " << vertexes[i] -> y <<"): ";
            for(int j = 0; j < graph[i].size(); j++) {
                Vertex* u = edges[graph[i][j]] -> u;
                Vertex* v = edges[graph[i][j]] -> v;
                // cout <<"(" << u -> x << "," << u -> y << ")-(" << v -> x << "," << v -> y << ") ";
            }
            // cout <<" | " << vertexes[i] -> u << vertexes[i] -> d << vertexes[i] -> r << vertexes[i] -> l;
            // cout <<endl;
        }
    }

    void print_rc(vector<vector<int>>& c_v, vector<vector<int>>& r_v)
    {
        for (int i = 0; i < r_v.size(); i++) {
            // cout <<"r_v " << i << ": ";
            for(int t = 0; t < r_v[i].size(); t++){
                // cout <<vertexes[r_v[i][t]] -> x << ',' << vertexes[r_v[i][t]] -> y << " ";
            }
            // cout <<endl;
        }
        for (int j = 0; j < c_v.size(); j++) {
            // cout <<"c_v " << j << ": ";
            for(int t = 0; t < c_v[j].size(); t++){
                // cout <<vertexes[c_v[j][t]] -> x << ',' << vertexes[c_v[j][t]] -> y << " ";
            }
            // cout <<endl;
        }
    }

    void print_edges(){
        // cout<<"edge num:" << edges.size() << endl;
        for (int i = 0; i < edges.size(); i++){
            Edge* e = edges[i];
            // cout <<"<(" << e -> u -> x << ", " << e -> u -> y << ") - (" << e -> v -> x << ", " << e -> v -> y << ")>" << endl;
        }
    }

};



#endif