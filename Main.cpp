#include <iostream>
#include <cstdio>
#include <fstream>
#include <ctime>

#include "Matrix.h"
#include "SegementTree.h"
#include "Graph.h"

void gen_rect_strassen(int n, int m, int p, int range)
{
    clock_t time_start;
    clock_t time_end;
    Matrix mat1(n, m, range);
    Matrix mat2(m, p, range);
    Matrix mat3(n, p);
    double t = 0.0;
    int times = 10;

    for(int i = 0; i < times; i++)
    {
        time_start = clock();
        if(!mat1.strassen(mat2, mat3)){ std::cout << "cannot multiply" << std::endl; }
        time_end = clock();
        t += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    }
    printf("(%d,%.3f)", m, t / times);

}

void gen_rect_naive(int n, int m, int p, int range)
{
    clock_t time_start;
    clock_t time_end;
    Matrix mat1(n, m, range);
    Matrix mat2(m, p, range);
    Matrix mat3(n, p);
    double t = 0.0;
    int times = 8;

    for(int i = 0; i < times; i++)
    {
        time_start = clock();
        if(!mat1.naive_multi(mat2, mat3)){ std::cout << "cannot multiply" << std::endl; }
        time_end = clock();
        t += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    }
    printf("(%d,%.3f)", m, t / times);

}

void gen_strassen_thr(int n, int range, int thr)
{
    clock_t time_start;
    clock_t time_end;
    Matrix mat1(n, n, range);
    Matrix mat2(n, n, range);
    mat1.change_threshold(thr);
    mat2.change_threshold(thr);
    Matrix mat3(n, n);
    double t = 0.0;
    int times = 8;

    for(int i = 0; i < times; i++)
    {
        time_start = clock();
        if(!mat1.strassen(mat2, mat3)){ std::cout << "cannot multiply" << std::endl; }

        time_end = clock();
        t += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    }
    printf("(%d,%.3f)", thr, t / times);
}

void gen_sparse_strassen(int n, int m, int p, int range, int rate)
{
    clock_t time_start;
    clock_t time_end;
    SparseMatrix mat1(n, m, range, rate);
    SparseMatrix mat2(m, p, range, rate);
    Matrix mat3(n, p);
    double t = 0.0;
    int times = 8;

    for(int i = 0; i < times; i++)
    {
        time_start = clock();
        if(!mat1.strassen(mat2, mat3)){ std::cout << "cannot multiply" << std::endl; }
        time_end = clock();
        t += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    }
    printf("(%d,%.3f)", rate, t / times);
}

void gen_sparse_naive(int n, int m, int p, int range, int rate)
{
    clock_t time_start;
    clock_t time_end;

    double t = 0.0;
    int times = 8;

    for(int i = 0; i < times; i++)
    {
        SparseMatrix mat1(n, m, range, rate);
        SparseMatrix mat2(m, p, range, rate);
        Matrix mat3(n, p);
        time_start = clock();
        if(!mat1.naive_sparse_multi(mat2, mat3)){ std::cout << "cannot multiply" << std::endl; }
        time_end = clock();
        t += (double)(time_end - time_start) / CLOCKS_PER_SEC;
    }
    printf("(%d,%.3f)", rate, t / times);
}

void time_strassen()
{
    int n = 1000, m_st = 50, gap = 50, m_ed = 1000;
    int range = 50;
    for(int m = m_st; m <= m_ed; m += gap)
    {
        gen_rect_naive(n, m, n, range);
    }
    for(int m = m_st; m <= m_ed; m += gap)
    {
        gen_rect_strassen(n, m, n, range);
    }


}

void time_square_compare()
{
    int n = 1000, m_st = 50, gap = 50, m_ed = 1000;
    int range = 50;

    for(int m = m_st; m <= m_ed; m += gap)
    {
        gen_rect_strassen(m, m, m, range);
    }
    cout << endl;
    for(int m = m_st; m <= m_ed; m += gap)
    {
        gen_rect_naive(m, m, m, range);

    }
}

void time_strassen_thr()
{
    int m_st = 16, gap = 16, m_ed = 256;
    int range = 50;

    for(int m = m_st; m <= m_ed; m += 16)
    {
        gen_strassen_thr(1000, range, m);
    }
}

void time_sparse(int rate)
{
    int n = 1000, m_st = 50, gap = 50, m_ed = 1000;
    int range = 50;

    for(int m = m_st; m <= m_ed; m += gap)
    {
        gen_sparse_strassen(m, m, m, range, rate);
    }
    cout << endl;
    for(int m = m_st; m <= m_ed; m += gap)
    {
        gen_sparse_naive(m, m, m, range, rate);

    }
}

void time_sparse_rate()
{
    int n = 1000, r_st = 0, gap = 5, r_ed = 100;
    int range = 50;

    for(int r = r_st; r <= r_ed; r += gap)
    {
        //gen_sparse_strassen(n, n, n, range, r);
    }
    cout << endl;
    for(int r = r_st; r <= r_ed; r += gap)
    {
        gen_sparse_naive(n, n, n, range, r);
    }
}



void test_strassen()
{

    Matrix mat1(8, 8, 8);
    Matrix mat2(8, 8, 8);
    Matrix mat3; // Result Matrix
    Matrix mat4; // Result Matrix
    mat1.print();
    mat2.print();
    printf("Mat1 equal Mat2: %d\n", mat1.equals(mat2));

    clock_t time_start = clock();
    clock_t time_end = clock();
    // naive multiply
    if(!mat1.naive_multi(mat2, mat3)){ std::cout << "cannot multiply 1" << std::endl; }
    time_end = clock();
    mat3.print();

    printf("Time taken: %.3f\n", (double)(time_end - time_start) / CLOCKS_PER_SEC);

    // strassen multiply
    time_start = clock();
    if(!mat1.strassen(mat2, mat4)){ std::cout << "cannot multiply 2" << std::endl; }
    time_end = clock();
    mat4.print();

    printf("Mat3 equals Mat4 %d\n", mat3.equals(mat4));
    printf("Time taken: %.3f\n", (double)(time_end - time_start) / CLOCKS_PER_SEC);
}

void test_sparse()
{

    SparseMatrix mat1(1000, 1000, 10, 95);
    SparseMatrix mat2(1000, 1000, 10, 95);
    Matrix mat3(1000, 1000);
    Matrix mat4(1000, 1000);
    Matrix mat5(1000, 1000);

    clock_t time_start = clock();
    clock_t time_end = clock();

    if(!mat1.strassen(mat2, mat3)){ std::cout << "cannot multiply 1" << std::endl; }

    time_end = clock();
    mat3.print_hash();

    printf("Time taken: %.3f\n", (double)(time_end - time_start) / CLOCKS_PER_SEC);

    // naive sparse
    time_start = clock();
    if(!mat1.naive_sparse_multi(mat2, mat4)){ std::cout << "cannot multiply 2" << std::endl; }
    time_end = clock();
    mat4.print_hash();

    printf("Mat3 equals Mat4 %d\n", mat3.equals(mat4));
    printf("Time taken: %.3f\n", (double)(time_end - time_start) / CLOCKS_PER_SEC);

    // fast sparse
    time_start = clock();
    if(!mat1.fast_sparse_multi(mat2, mat5)){ std::cout << "cannot multiply 3" << std::endl; }
    time_end = clock();
    mat5.print_hash();

    printf("Mat3 equals Mat5 %d\n", mat3.equals(mat5));
    printf("Time taken: %.3f\n", (double)(time_end - time_start) / CLOCKS_PER_SEC);
}

bool threeDMultiply(Matrix& mat1, Matrix& mat2, Matrix& res){
    if(!mat1.valid_multi(mat2)) return false;

    Graph ga(mat1);
    Graph gb(mat2);
    vector<Rect*> rectas = ga.rects;
    vector<Rect*> rectbs = gb.rects;
    // delete ga;
    // cout << rectas.size() << endl;
    // cout << rectbs.size() << endl;
    // delete gb;

    // cout << "START!" << endl;
    // mat1.print();

    vector<Interval<int, Rect*>> intervals;
    for (int i = 0; i < rectas.size(); i++){
        intervals.push_back(Interval<int, Rect*>(rectas[i] -> k1, rectas[i] -> k2, rectas[i]));
    }

    // cout << "INIT INTERVAL FINISHED!" << endl;
    // mat1.print();

    IntervalTree<int, Rect*> tree = IntervalTree<int, Rect*>(std::move(intervals), 30);
    // cout << "INTERVAL TREE CONSTRUCTED!" << endl;
    // mat1.print();

    vector<vector<Interval<int, int>>> start;
    for(int i = 0; i < mat1.rows; i++) start.push_back(vector<Interval<int, int>>());

    vector<vector<Interval<int, int>>> end;
    for(int i = 0; i < mat1.rows; i++) end.push_back(vector<Interval<int, int>>());
    // cout << "START END CONSTRUCTED!" << endl;
    // mat1.print();

    vector<Interval<int, Rect*>> result;
    for (int i = 0; i < rectbs.size(); i++){
        // cout << "i:" << i << endl;
        int k1_ = rectbs[i] -> i1;
        int k2_ = rectbs[i] -> i2;
        int j1 = rectbs[i] -> k1;
        int j2 = rectbs[i] -> k2;
        int h1_ = rectbs[i] -> h1;
        int h2_ = rectbs[i] -> h2;
        result = tree.findOverlapping(k1_, k2_);
        for(int j = 0; j < result.size(); j++) {
            auto itv = result[j];
            int k1 = itv.start;
            int k2 = itv.stop;
            Rect* rect = itv.value;
            int val = ((rect -> h2) - (rect -> h1)) * (h2_ - h1_) * mat1.card(k1, k2, k1_, k2_);
            // cout << "update:" << val << endl;
            auto temp = new Interval<int, int>(rect ->i1, rect -> i2, val);
            // cout << "k1k2: "<< k1 << ", " << k2 << endl;
            // cout << "pushing: "<< "j1_j2: " << j1 << "_" << j2 << ": " << rect -> i1 << " " << rect -> i2 << " " << val << endl;
            start[j1].push_back(*temp);
            end[j2].push_back(*temp);
        }
    }
    // cout << "START SEGEMENT TREE!" << endl;
    // mat1.print();
    SegementTree sTree = SegementTree(mat1.rows + 1);
    sTree.build(1, 0, mat1.rows);
    for(int j = 0; j < mat1.rows; j++){
        // cout << "j:" << j << endl;
        if (j != 0){
            for(int k = 0; k < end[j - 1].size(); k++){
                sTree.update(1, end[j - 1][k].start, end[j - 1][k].stop, -end[j - 1][k].value);
                if(end[j - 1][k].value!= 0){
                    // cout << "delete: " << end[j - 1][k].start << "_" << end[j - 1][k].stop << ": " << -end[j - 1][k].value << endl;
                }
            }

        }
        // cout << "j:" << j << endl;
        for(int k = 0; k < start[j].size(); k++){
            sTree.update(1, start[j][k].start, start[j][k].stop, start[j][k].value);
            if(start[j][k].value != 0){
                // cout << "update: " << start[j][k].start << "_" << start[j][k].stop << ": " << start[j][k].value << endl;
            }

        }

        // cout << "j:" << j << endl;
        for(int k = 0; k < mat1.rows; k++){

            int t = sTree.query(1, k, k);
            // cout << "t: " << t << endl;
            res.data[k][j]= t;
        }
        //sTree.update(start[j])
    }

    // cout << "THREED FINISHED!" << endl;
    // mat1.print();
    // mat2.print();


    //intervals.push_back(Interval<int, Rect*>(rectbs[i] -> k1, rectbs[i] -> k2, rectbs[i]));

}

void gen_3D(int n, int range, int num)
{
    clock_t time_start;
    clock_t time_end;
    Matrix mat1(n, n, range, num);
    Matrix mat2(n, n, range, num);

    Matrix mat3(n, n);
    double t = 0.0;
    int times = 8;

    for(int i = 0; i < times; i++)
    {
        time_start = clock();
        threeDMultiply(mat1, mat2, mat3);

        time_end = clock();
        t += (double)(time_end - time_start) / CLOCKS_PER_SEC;
        // cout << t/(i+1) << endl;
    }
    printf("(%d,%.3f)", n, t / times);
}

void test_3D(){
    int m_st = 100, gap = 100, m_ed = 800;
    int range = 100;

    cout << "3d" << endl;
    for(int m = m_st; m <= m_ed; m += 100)
    {
        gen_3D(m, range, 32);
    }
    cout << "\n3dm" << endl;
    for(int m = m_st; m <= m_ed; m += 100)
    {
        gen_3D(m, range, (int)(pow((float)m/2.0, 3.0/4.0)));
    }

    cout << "\nnaive" << endl;
    for(int m = m_st; m <= m_ed; m += 100)
    {
        gen_rect_naive(m, m, m, range);
    }
    cout << "\nstrassen" << endl;
    for(int m = m_st; m <= m_ed; m += 100)
    {
        gen_strassen_thr(m, range, 64);
    }
}

int main() {
    srand(time(NULL));
    /*cout << "Time strassen thr" << endl;
    time_strassen_thr();
    cout << "\n\nTime sparse size compare rate 1%" << endl;
    time_sparse(1);*/
    //cout << "\nTime sparse rate compare" << endl;
    //time_sparse_rate();
    //test_sparse();
    //time_strassen();

    Matrix m(1000, 1000, 10, 10);
    Matrix m1(1000, 1000, 10, 10);
    // m.read_from_file("../matrix1");
    // Matrix m1(4, 4);
    // m1.read_from_file("../matrix1");

    Matrix r(1000, 1000);
    Matrix r1(100, 100);
    // m.print();
    // m1.print();
    // gg(m);
    // m.print();
    // vector<Rect*> rects = g.rects;
    test_3D();
    //threeDMultiply(m, m1, r);

    cout << "Finished!" << endl;



    return 0;

} // main


