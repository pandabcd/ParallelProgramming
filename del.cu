#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<stdio.h>
#include<cuda.h>
#include<math.h>
#include<fstream>
#include<cuda_runtime.h>
#include<cooperative_groups.h>
#include<cuda_runtime_api.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

using namespace std;

// thrust::host_vector<thrust::device_vector<int>> b;

__global__
void func_check(int** b){
    printf("size of vector is %d", b[0][0]);
}


int main(){
    thrust::host_vector<int> a;
    a.push_back(100);
    thrust::device_vector<thrust::device_vector<int>> b;
    b.push_back(a)

    cout << a.size();
    // thrust::device_ptr<int> dp = b.data();
    int **dp = thrust::raw_pointer_cast(b.data());
    func_check<<<1,1>>>(dp);
}





// __global__
// void func_check(){
// 	printf("size of vector is %d", b.size());
// }


// int main(){
// 	thrust::host_vector<thrust::host_vector<int>> a;


// 	thrust::host_vector<int> temp;
// 	temp.push_back(1);
// 	temp.push_back(2);
// 	a.push_back(temp);
// 	temp.clear();
// 	temp.push_back(3);
// 	a.push_back(temp);

// 	thrust::host_vector<thrust::device_vector<int>> b = a;

// 	for(int i=0;i<a.size();i++){
// 		for(int j=0;j<a[i].size();j++){
// 			cout<<b[i][j]<<endl;
// 		}
// 		cout<<endl;
// 	}

// 	func_check<<<1,1>>>();

// }