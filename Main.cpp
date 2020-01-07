#include <iostream>
#include <cstdio>
#include <fstream>
#include <ctime>
#include "Matrix.h"

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


int main() {
    srand(time(NULL));
    /*cout << "Time strassen thr" << endl;
    time_strassen_thr();
    cout << "\n\nTime sparse size compare rate 1%" << endl;
    time_sparse(1);*/
    cout << "\nTime sparse rate compare" << endl;
    time_sparse_rate();
    //test_sparse();
    //time_strassen();

    return 0;

} // main


