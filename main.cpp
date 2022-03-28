#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char **argv)
{
    string data_image_path = argv[1];
    string query_image_path = argv[2];

    // reading a matrix from a file with first line m*n
    ifstream data_image_file(data_image_path);
    ifstream query_image_file(query_image_path);
    int data_image_m, data_image_n;
    int query_image_m, query_image_n;
    data_image_file >> data_image_m >> data_image_n;
    query_image_file >> query_image_m >> query_image_n;

    // reading a matrix from a file
    // with first line m*n

    int data_image_matrix[data_image_m][data_image_n][3];
    int query_image_matrix[query_image_m][query_image_n][3];

    for (int i = 0; i < data_image_m; i++)
    {
        for (int j = 0; j < data_image_n; j++)
        {
            data_image_file >> data_image_matrix[i][j][0] >> data_image_matrix[i][j][1] >> data_image_matrix[i][j][2];
        }
    }

    for (int i = 0; i < query_image_m; i++)
    {
        for (int j = 0; j < query_image_n; j++)
        {
            query_image_file >> query_image_matrix[i][j][0] >> query_image_matrix[i][j][1] >> query_image_matrix[i][j][2];
        }
    }

    cout << "data_image_m: " << data_image_m << endl;
    cout << "data_image_n: " << data_image_n << endl;
    cout << "query_image_m: " << query_image_m << endl;
    cout << "query_image_n: " << query_image_n << endl;

    cout << "data_image_matrix: " << data_image_matrix[0][0][0] << endl;
    cout << "query_image_matrix: " << query_image_matrix[0][0][0] << endl;

    return 0;
}