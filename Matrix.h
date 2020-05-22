#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <random>
#include <ctime>
#include <fstream>

#include "IntervalTree.h"
#include "SegementTree.h"


using namespace std;

#ifndef MATRIX_H
#define MATRIX_H

class Matrix {
private:


public:
    int rows;
    int cols;
    int threshold = 64;
    int** data;
    Matrix() // Default Constructor - Null or Empty Matrix
    {rows = 0; cols = 0;}
    Matrix(const int rows, const int cols) // Empty Matrix With Defined Size
    {
        this -> rows = rows;
        this -> cols = cols;
        allocate_space();
    }
    // square matrix
    Matrix(const int size)
    {
        this -> rows = size;
        this -> cols = size;
        allocate_space();
    }
    Matrix(const int rows, const int cols, const int range) // Matrix With Defined Size Random Value
    {
        this -> rows = rows;
        this -> cols = cols;
        allocate_space();
        for (int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                data[i][j] = (rand() % (range + 1));

    }

    Matrix(const int rows, const int cols, const int range, const int n) // Matrix With Defined Size Random Value with fewer rectangler 
    {
        this -> rows = rows;
        this -> cols = cols;
        allocate_space();
        int gapr = rows / n;
        int gapc = cols / n;
        for (int i = 0; i < rows; i+=gapr)
            for(int j = 0; j < cols; j+=gapc){
                int t = (rand() % (range + 1));
                for (int k1 = 0; k1 < gapr; k1++)
                    for (int k2 = 0; k2 < gapc; k2++)
                        if(i + k1 < rows && j + k2 < rows) data[i + k1][j + k2] = t;
            }

    }

    ~Matrix()
    {
        for (int i = 0; i < rows; i++) free(data[i]);
        free(data);
    }
    void change_threshold(int thr)
    {
        threshold = thr;
    }

    void fill(const int rows, const int cols, const int range)
    {
        this -> rows = rows;
        this -> cols = cols;
        for (int i = 0; i < rows; i++)
            for(int j = 0; j < cols; j++)
                data[i][j] = (rand() % (range + 1));
    }

    void allocate_space()
    {
        data = (int **)calloc(rows + 2, sizeof(int*));
        for(int i = 0; i < rows + 2; i++) data[i] = (int *)calloc(cols + 2, sizeof(int));
    }

    void print()
    {
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                // printf("%d ", data[i][j]);
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    int print_hash()
    {
        int sum = 0;
        for(int i = 0; i < rows; i+=2)
        {
            for(int j = 0; j < cols; j+=2)
            {
                sum += data[i][j];
            }
        }
        printf("%d\n", sum);
        return sum;
    }

    bool equals(Matrix& mat2){
        if(cols != mat2.cols || rows != mat2.rows) return false;
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                if(data[i][j] != mat2.data[i][j])
                    return false;
            }
        }
        return true;
    }

    void duplicate(Matrix& mat2, int st_row, int ed_row, int st_col, int ed_col)
    {
        mat2 = Matrix(ed_row - st_row + 1, ed_col - st_col + 1);
        for(int i = st_row; i <= ed_row; i++)
            for (int j = st_col; j <= ed_col; j++)
                mat2.data[i - st_row][j - st_col] = data[i][j];
    }

    int card(int k1, int k2, int j1, int j2){
        if (k2 < j1 || j2 < k1) return 0;
        //if (k2 >= j1) return 
        int l = k1 > j1 ? k1 : j1;
        int r = k2 < j2 ? k2 : j2;
        //cout << "card: " << k1 << "_" << k2 << ", " << j1 << "_" << j2 << ": " << r-l+1 << endl;
        return r - l + 1;
    }

    bool naive_multi(Matrix& mat2, Matrix& res);
    bool strassen(Matrix& mat2, Matrix& res);
    void strassen_helper(Matrix& mat1, Matrix& mat2, Matrix& res);
    bool threeDMultiply(Matrix& mat2, Matrix& res);

    bool add(Matrix& mat2, Matrix& res);
    bool subtract(Matrix& mat2, Matrix& res);
    bool valid_multi(Matrix& mat2);
    bool same_shape(Matrix& mat2);
    bool is_square() { return cols == rows;}
    // Matrix from file
    void read_from_file(const char * filename) {
        std::ifstream file;
        file.open(filename);
        file >> rows >> cols;
        allocate_space();
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++){
                file >> data[i][j];
                // std::cout << data[i][j] << std::endl;
            }
        file.close();
    }

}; // Matrix

bool Matrix::naive_multi(Matrix& mat2, Matrix& res)
{
    if (!valid_multi(mat2)) return false;
    //res = Matrix(this -> rows, mat2.cols);

    for(int i = 0; i < res.rows; i++)
        for (int j = 0; j < res.cols; j++){
            res.data[i][j] = 0;
            for(int k = 0; k < this -> cols; k++)
                res.data[i][j] += data[i][k] * mat2.data[k][j];
        }

    return true;
}

bool Matrix::add(Matrix& mat2, Matrix& res)
{
    if (!same_shape(mat2)) return false;

    for(int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            res.data[i][j] = data[i][j] + mat2.data[i][j];

    return true;
}

bool Matrix::subtract(Matrix& mat2, Matrix& res)
{
    if (!same_shape(mat2)) return false;

    for(int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            res.data[i][j] = data[i][j] - mat2.data[i][j];

    return true;
}

bool Matrix::valid_multi(Matrix& mat2)
{
    return ( this -> cols == mat2.rows );
}

bool Matrix::same_shape(Matrix& mat2)
{
    return ( this -> cols == mat2.cols && this -> rows == mat2.rows );
}

bool Matrix::strassen(Matrix& mat2, Matrix& res)
{
    // ensure both are the same size square matrix
    if(!valid_multi(mat2)) return false;
    // res = Matrix(rows, mat2.cols);


    // create padding matrix fills with 0
    Matrix pA(rows, cols);
    Matrix pB(mat2.rows, mat2.cols);
    Matrix pC(rows, mat2.cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            pA.data[i][j] = data[i][j];

    for (int i = 0; i < mat2.rows; i++)
        for (int j = 0; j < mat2.cols; j++)
            pB.data[i][j] = mat2.data[i][j];

    strassen_helper(pA, pB, pC);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < mat2.cols; j++)
        {
            res.data[i][j] = pC.data[i][j];
        }

    return true;
}

void Matrix::strassen_helper(Matrix& mat1, Matrix& mat2, Matrix& res)
{
    // int threshold = 64;
    // when size of matrix is small, use the naive method to avoid high computation overhead
    if (mat1.cols <= threshold || mat1.rows <= threshold || mat2.cols <= threshold){
        mat1.naive_multi(mat2, res);
        return;
    }

    int new_cols = mat2.cols / 2;
    if(mat2.cols % 2) new_cols++;

    int new_rows = mat1.rows / 2;
    if(mat1.rows % 2) new_rows++;

    int new_inter = mat1.cols / 2;
    if(mat1.cols % 2) new_inter++;

    Matrix a11(new_rows, new_inter), a12(new_rows, new_inter), a21(new_rows, new_inter), a22(new_rows, new_inter);
    Matrix b11(new_inter, new_cols), b12(new_inter, new_cols), b21(new_inter, new_cols), b22(new_inter, new_cols);


    Matrix ar(new_rows, new_inter), br(new_inter, new_cols), pr(new_rows, new_cols), pr1(new_rows, new_cols), pr2(new_rows, new_cols);

    for (int i = 0; i < new_rows; i++) {
        for (int j = 0; j < new_inter; j++) {
            a11.data[i][j] = mat1.data[i][j];
            a12.data[i][j] = mat1.data[i][j + new_inter];
            a21.data[i][j] = mat1.data[i + new_rows][j];
            a22.data[i][j] = mat1.data[i + new_rows][j + new_inter];
        }
    }

    for (int i = 0; i < new_inter; i++) {
        for (int j = 0; j < new_cols; j++) {
            b11.data[i][j] = mat2.data[i][j];
            b12.data[i][j] = mat2.data[i][j + new_cols];
            b21.data[i][j] = mat2.data[i + new_inter][j];
            b22.data[i][j] = mat2.data[i + new_inter][j + new_cols];
        }
    }

    // calcuate middle matrixs
    Matrix p1(new_rows, new_cols);
    a11.add(a22, ar);
    b11.add(b22, br);
    strassen_helper(ar, br, p1);  // p1 = (a11 + a22) * (b11 + b22)

    Matrix p2(new_rows, new_cols);
    a21.add(a22, ar);
    strassen_helper(ar, b11, p2); // p2 = (a21 + a22) * b11

    Matrix p3(new_rows, new_cols);
    b12.subtract(b22, br);
    strassen_helper(a11, br, p3); // p3 = (b12 - b22) * a11

    Matrix p4(new_rows, new_cols);
    b21.subtract(b11, br);
    strassen_helper(a22, br, p4); // p4 = (b21 - b11) * a22

    Matrix p5(new_rows, new_cols);
    a11.add(a12, ar);
    strassen_helper(ar, b22, p5); // p5 = (a11 + a12) * b22

    Matrix p6(new_rows, new_cols);
    a21.subtract(a11, ar);
    b11.add(b12, br);
    strassen_helper(ar, br, p6);  // p6 = (a21 - a11) * (b11 + b12)

    Matrix p7(new_rows, new_cols);
    a12.subtract(a22, ar);
    b21.add(b22, br);
    strassen_helper(ar, br, p7);  // p7 = (a12 - a22) * (b21 + b22)

    // calculate the result
    Matrix c11(new_rows, new_cols), c12(new_rows, new_cols), c21(new_rows, new_cols), c22(new_rows, new_cols);
    p1.add(p4, pr);
    pr.add(p7, pr1);
    pr1.subtract(p5, c11); // c11 = p1 + p4 + p7 - p5
    p3.add(p5, c12); // c12 = p3 + p5
    p2.add(p4, c21); // c21 = p2 + p4
    p1.subtract(p2, pr);
    pr.add(p3, pr1);
    pr1.add(p6, c22); // c22 = p1 - p2 + p3 + p6

    for (int i = 0; i < new_rows ; i++) {
        for (int j = 0 ; j < new_cols ; j++) {
            res.data[i][j] = c11.data[i][j];
            res.data[i][j + new_cols] = c12.data[i][j];
            res.data[i + new_rows][j] = c21.data[i][j];
            res.data[i + new_rows][j + new_cols] = c22.data[i][j];
        }
    }

}

bool Matrix::threeDMultiply(Matrix& mat2, Matrix& res)
{
    
}

//#include "Matrix.inl"

struct Node
{
    int data;
    int id;
    Node(){}
    Node(int data, int id)
    {
        this -> data = data;
        this -> id = id;
    }
    bool operator < (const Node& b) const
    {
        return (data < b.data);
    }
    bool operator > (const Node& b) const
    {
        return (data > b.data);
    }
};


class SparseMatrix: public Matrix
{
public:
    std::vector<Node>* row_num;
    std::vector<Node>* col_num;
    int num;
    SparseMatrix(){}
    SparseMatrix(const int rows, const int cols): Matrix(rows, cols)
    {
        row_num = new vector<Node> [rows];
        col_num = new vector<Node> [cols];
        num = 0;
    }
    SparseMatrix(const int rows, const int cols, const int range, const int rate): SparseMatrix(rows, cols)
    {
        for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            int t = 0;
            if (rand() % 100 < rate)
                data[i][j] = (rand() % (range + 1));
        }
        count_num();
    }
    void read_from_file(char * filename)
    {
        Matrix::read_from_file(filename);
        row_num = new vector<Node> [rows];
        col_num = new vector<Node> [cols];
        num = 0;
        count_num();
    }

    void duplicate(SparseMatrix& mat2, int st_row, int ed_row, int st_col, int ed_col)
    {
        mat2 = SparseMatrix(ed_row - st_row + 1, ed_col - st_col + 1);
        for(int i = st_row; i <= ed_row; i++)
            for (int j = st_col; j <= ed_col; j++)
                mat2.data[i - st_row][j - st_col] = data[i][j];
        count_num();
    }

    void count_num()
    {
        for(int i = 0; i < rows; i++)
        {
            for(int j = 0; j < cols; j++)
            {
                int d = data[i][j];
                if(d != 0){
                    row_num[i].push_back(Node(d, j));
                    col_num[j].push_back(Node(d, i));
                    num++;
                }
            }
        }
    }

    void print_sparse()
    {
        cout << "Total Element: " << num << endl;
        for(int i = 0; i < rows; i++)
        {
            auto s = row_num[i].size();
            cout << "Row " << i << " count " << s << ": ";
            for(int j = 0; j < s; j++)
                cout << row_num[i][j].data << "," << row_num[i][j].id << " ";
            cout << endl;
        }

        for(int i = 0; i < cols; i++)
        {
            auto s = col_num[i].size();
            cout << "Col " << i << " count " << s << ": ";
            for(int j = 0; j < s; j++)
                cout << col_num[i][j].data << "," << col_num[i][j].id << " ";
            cout << endl;
        }
    }

    bool naive_sparse_multi(Matrix& mat2, Matrix& res);
    bool fast_sparse_multi(SparseMatrix& mat2, Matrix& res);


};

bool SparseMatrix::naive_sparse_multi(Matrix& mat2, Matrix& res)
{
    if (!valid_multi(mat2)) return false;
    // res = Matrix(rows, mat2.cols);

    for(int i = 0; i < res.rows; i++) {
        int s = row_num[i].size();
        auto n = row_num[i];
        for (int j = 0; j < res.cols; j++){
            res.data[i][j] = 0;
            for(int k = 0; k < s; k++) {
                res.data[i][j] += n[k].data * mat2.data[n[k].id][j];
            }
        }
    }
    return true;
}

bool SparseMatrix::fast_sparse_multi(SparseMatrix& mat2, Matrix& res)
{
    if (!valid_multi(mat2)) return false;
    // res = Matrix(rows, mat2.cols);
    Node val[cols];
    int sum = 0;
    for (int i = 0; i < cols; i++)
    {
        val[i].data = col_num[i].size() * row_num[i].size();
        sum += val[i].data;
        val[i].id = i;
    }
    sort(val, val + cols, greater<Node>());

    int split = cols * 9 / 10;
    double min_com = DBL_MAX;
    cout << "Sum: " << sum << endl;


    for (int i = 0; i < cols; i++)
    {
        double com = 0.0;
        com += (pow((double)rows, 2.8074)) * ((double)i / (double)rows); // strassen time
        cout << com << endl;
        com += (double)sum;
        cout << com << " " << sum << endl;
        if(com < min_com){
            min_com = com;
            split = i;
            cout << i << " - " << min_com << endl;
        }
        sum -= val[i].data;
    }
    split = cols/2;

    Matrix AI(rows, split); // fast rect
    Matrix BI(split, rows); // fast rect
    Matrix CI(rows, rows); // fast rect
    SparseMatrix AJ(rows, rows - split); // naive sparse
    SparseMatrix BJ(rows - split, rows); // naive sparse
    Matrix CJ(rows, rows); // fast rect


    for (int i = 0; i < split; i++) { // for fast rect multi
        // cout << i << ": " << val[i].data << " " << val[i].id << endl;
        for (int j = 0; j < rows; j++) {
            AI.data[j][i] = data[j][val[i].id];
            BI.data[i][j] = mat2.data[val[i].id][j];
        }
    }

    for (int i = split; i < cols; i++) { // for fast rect multi
        // cout << i << ": " << val[i].data << " " << val[i].id << endl;
        for (int j = 0; j < rows; j++) {
            AJ.data[j][i - split] = data[j][val[i].id];
            BJ.data[i - split][j] = mat2.data[val[i].id][j];
        }
    }
    AJ.count_num();
    BJ.count_num();

    clock_t time_start = clock();
    AI.strassen(BI, CI);
    AJ.naive_sparse_multi(BJ, CJ);
    clock_t time_end = clock();
    printf("Time taken using two methods: %.3f\n", (double)(time_end - time_start) / CLOCKS_PER_SEC);

    CI.add(CJ, res);

    return true;
}



#endif // MATRIX_H
