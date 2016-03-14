//Jianan Wang 800898865
//This is a code to read a txt file that contains a sequence of numbers separated with semicolon, 
//sort it and output it to a txt file named answer.txt at the same directory.


#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;

void exchange (int& a, int& b){
	int c=a;
	a=b;
	b=c;
}
int partition (int a[], int p, int r)
{
  int x=a[r];
  int i, temp;
  i=p-1;
  for(int j=p; j<=(r-1); j++){
	  if a[j]<=x{
		  i++;
		  exchange(a[i],a[j])
	  }
	  exchange(a[i+1],a[r])
  return r;
}
void quicksort(int& a[],int p,int r)}{
	if(p<r){
		q=partition(a,p,r);
		quicksort(a,p,(q-1));
		quicksort(a,(q+1),r);
	}
}
int main() {
	fstream openfile;
	string filename;
	cout<< "Enter the file name that you want to sort?"<<endl;
	cin>>filename;
	openfile.open(filename);
	vector<int> numbers;	
// read numbers from file to a vector of number.
	if (openfile.is_open()) {
		string line;
		while (getline(openfile, line, ';')) {
			istringstream iss(line);
			int value;
			while (iss>>value) {
				numbers.push_back(value);
			};
		}
	}
//conduct the sorting.
quicksort(numbers,0,numbers.size())
//output to the answer.txt file
	ofstream outfile;
	outfile.open("answer.txt");
	for(int i=0;i<numbers.size();i++){
		outfile<<numbers[i]<<";";
	}
	outfile.close();
	cout<<" your output is stored in the file named answer.txt"<< endl;
return 0;
}

